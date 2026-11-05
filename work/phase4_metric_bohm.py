"""
Phase 4, Schritt M (korrigiert): Kinetische Metrik, Massdichte,
Laplace-Beltrami-Operator und Signature des getrimmten gemeinsamen Modells.

Metrikbestimmung per Impuls-Hessian von Hbar (S19.7); die per Hand abgelesenen
Komponenten werden als Gegenprobe gefuehrt und aus dem Hessian korrigiert.

Ausgabe: work/phase4_metric_bohm.json
"""

import json
import os

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 40
OUT = os.path.join("work", "phase4_metric_bohm.json")
res = {}


def dump():
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False, sort_keys=False)


# ----------------------------------------------------------------------------
# 0. Koordinaten, Impulse, Hbar aus S19.7
# ----------------------------------------------------------------------------
a0, a1, ac, c, q = sp.symbols("alpha_0 alpha_1 alpha_c c q", real=True)
P0, P1, Pc, p0, pq = sp.symbols("P_0 P_1 P_c p_0 p_q", real=True)
beta, m, hbar = sp.symbols("beta m hbar", positive=True)

Qs = [a0, a1, ac, c, q]
Pis = [P0, P1, Pc, p0, pq]

z3 = 3 * a1
I0 = sp.besseli(0, z3)
I1 = sp.besseli(1, z3)
I2 = sp.besseli(2, z3)
Cz1 = (sp.besseli(0, a1) - sp.besseli(2, a1)) / 2  # C(alpha_1), Argument alpha_1

Hkin_geo = -(P0**2 * I0 - 4 * P0 * P1 * I1 + 4 * P1**2 * (I0 + I2) / 2 + 4 * Pc**2 * (I0 - I2) / 2)
Hkin_mat = p0**2 * I0 - 4 * p0 * pq * I1 + 4 * pq**2 * (I0 + I2) / 2
Vconf = (
    -sp.exp(a0) / (6 * beta) * a1**2 * Cz1
    + sp.exp(a0) / 2 * q**2 * Cz1
    + m**2 / 2 * sp.exp(3 * a0) * (c**2 * I0 + 2 * c * q * I1 + q**2 * (I0 + I2) / 2)
)
Hbar = beta / 2 * sp.exp(-3 * a0) * (Hkin_geo + Hkin_mat) + Vconf

Hom = beta / 2 * sp.exp(-3 * a0) * (-P0**2 + p0**2) + m**2 / 2 * sp.exp(3 * a0) * c**2
hom_diff = sp.simplify(Hbar.subs({a1: 0, q: 0, Pc: 0, P1: 0, pq: 0}) - Hom)
res["anchor_homogendifferenz"] = {"differenz": str(hom_diff), "gleich_null": hom_diff == 0}

# ----------------------------------------------------------------------------
# 1. Inverse Konfigurationsmetrik G^{AB} = d^2 Hbar / dPi_A dPi_B
# ----------------------------------------------------------------------------
G = sp.Matrix(5, 5, lambda i, j: sp.simplify(sp.diff(Hbar, Pis[i], Pis[j])))
e3 = sp.exp(-3 * a0)
Ghand_falsch = sp.zeros(5, 5)
Ghand_falsch[0, 0] = -beta * e3 * I0
Ghand_falsch[0, 1] = Ghand_falsch[1, 0] = 2 * beta * e3 * I1
Ghand_falsch[1, 1] = -beta * e3 * (I0 + I2)  # Faktor-2-Fehler der Vorversion
Ghand_falsch[2, 2] = -beta * e3 * (I0 - I2)  # Faktor-2-Fehler der Vorversion
Ghand_falsch[3, 3] = beta * e3 * I0
Ghand_falsch[3, 4] = Ghand_falsch[4, 3] = -2 * beta * e3 * I1
Ghand_falsch[4, 4] = beta * e3 * (I0 + I2)  # Faktor-2-Fehler der Vorversion
res["G"] = {
    "komponenten_hessian": {
        f"{i + 1}{j + 1}": str(G[i, j]) for i in range(5) for j in range(5) if G[i, j] != 0
    },
    "vorversion_handfehler_entdeckt": bool((sp.expand(G - Ghand_falsch) != sp.zeros(5, 5))),
    "blockstruktur_keine_geometrie_materie_mischung": all(
        sp.simplify(G[i, j]) == 0 for i in range(3) for j in range(3, 5)
    ),
    "matter_ist_minus_geometrie_block01": bool(
        sp.simplify(G[3, 3] + G[0, 0]) == 0
        and sp.simplify(G[3, 4] + G[0, 1]) == 0
        and sp.simplify(G[4, 4] + G[1, 1]) == 0
    ),
    "kanonische_geschwindigkeit_identity": [
        bool(sp.simplify(sp.diff(Hbar, Pis[i]) - sum(G[i, j] * Pis[j] for j in range(5))) == 0)
        for i in range(5)
    ],
}


def hom_limit(expr):
    """alpha_1 -> 0 unter Beachtung der hebbaren Singulaeritaet I_1(3 a_1)/(3 a_1) -> 1/2."""
    try:
        val = sp.simplify(sp.limit(sp.simplify(expr), a1, 0))
    except Exception:
        val = sp.nan
    if val == sp.nan:
        num = sp.N(expr.subs({a0: sp.Float(1, 30), a1: sp.Float(1, 60) * sp.Rational(1, 10 ** 12), beta: sp.Float(1, 30)}))
        val = sp.Float(num)
    return val


res["G"]["homogener_block_alpha1_0"] = {
    f"{i + 1}{i + 1}": str(hom_limit(G[i, i])) for i in range(5)
}
res["G"]["homogener_block_exakt_gegen_betae3"] = [
    bool(sp.simplify(sp.nsimplify(hom_limit(G[i, i]) - [1, 2, 2, 1, 2][i] * beta * sp.exp(-3 * a0) * (-1 if i < 3 else 1))) == 0)
    for i in range(5)
]
dump()

# ----------------------------------------------------------------------------
# 2. Determinante, Massdichte mu, homogener Anschluss
# ----------------------------------------------------------------------------
detGinv = sp.factor(G.det())
Dp = sp.simplify(I0 * (I0 + I2) - 2 * I1**2)  # Delta'(3 alpha_1)
pred_det = sp.simplify(-8 * beta**5 * sp.exp(-15 * a0) * (I0 - I2) * Dp**2)
res["determinante"] = {
    "det_Ginv_faktoriert": str(detGinv),
    "Delta_strich": str(Dp),
    "vorhersage": str(pred_det),
    "differenz_null": sp.expand(detGinv - pred_det) == 0,
}

mu = beta ** sp.Rational(-5, 2) * sp.exp(sp.Rational(15, 2) * a0) / (
    sp.sqrt(8) * sp.sqrt(I0 - I2) * Dp
)
mu_check = sp.simplify(mu ** (-2) - (-detGinv))
res["massdichte"] = {
    "mu_geschlossen": str(sp.simplify(mu)),
    "check_mu_quadrat_reziprok_det_null": mu_check == 0,
    "homogener_limes": str(hom_limit(mu)),
    "homogener_limes_geteilt_durch_betaregel": str(
        sp.simplify(hom_limit(mu) / (beta * sp.exp(-3 * a0)) ** sp.Rational(-5, 2))
    ),
    "betaregel_mu5": str((beta * sp.exp(-3 * a0)) ** sp.Rational(-5, 2) / sp.sqrt(8)),
}

# b^A = (1/mu) d_B (mu G^{AB}); homogene Grenze direkt bei alpha_1 = 0
Ghom = sp.diag(*[hom_limit(G[i, i]) for i in range(5)])
muh = hom_limit(mu)
bA_hom = [sp.simplify(sum(sp.diff(muh * Ghom[i, j], Qs[j]) for j in range(5)) / muh) for i in range(5)]
bA = [sp.simplify(sum(sp.diff(mu * G[i, j], Qs[j]) for j in range(5)) / mu) for i in range(5)]
res["laplace_koeffizienten"] = {
    "b_alpha_0_allgemein": str(bA[0]),
    "b_alpha_0_homogen": str(bA_hom[0]),
    "b_uebrige_homogen": [str(b) for b in bA_hom[1:]],
    "gamma_normalform_homogen": str(sp.simplify(bA_hom[0] / Ghom[0, 0])),
    "gamma_gleich_9_halb": bool(sp.simplify(bA_hom[0] / Ghom[0, 0] - sp.Rational(9, 2)) == 0),
    "b_alpha_0_homogen_ist_minus_9_halb_beta_e3": bool(
        sp.simplify(bA_hom[0] + sp.Rational(9, 2) * beta * sp.exp(-3 * a0)) == 0
    ),
    "gamma_regel_drei_halb_N_minus_2": bool(
        sp.simplify(bA_hom[0] / Ghom[0, 0] - sp.Rational(3, 2) * (5 - 2)) == 0
    ),
}
dump()

# ----------------------------------------------------------------------------
# 3. Delta'-Positivitaet, Signature, Entartung
# ----------------------------------------------------------------------------


def bessel_triple(z):
    z = mp.mpf(z)
    return (mp.besseli(0, z), mp.besseli(1, z), mp.besseli(2, z))


def Dp_direct(z):
    i0, i1, i2 = bessel_triple(z)
    return i0 * (i0 + i2) - 2 * i1**2


def Ginv_numeric(a0v, a1v, betav):
    i0, i1, i2 = bessel_triple(3 * mp.mpf(a1v))
    pre = mp.mpf(betav) * mp.e ** (-3 * mp.mpf(a0v))
    M = np.zeros((5, 5), dtype=float)
    M[0, 0] = float(-pre * i0)
    M[0, 1] = M[1, 0] = float(2 * pre * i1)
    M[1, 1] = float(-2 * pre * (i0 + i2))
    M[2, 2] = float(-2 * pre * (i0 - i2))
    M[3, 3] = float(pre * i0)
    M[3, 4] = M[4, 3] = float(-2 * pre * i1)
    M[4, 4] = float(2 * pre * (i0 + i2))
    return M


scan = []
for z in ["0", "0.001", "0.1", "0.5", "1", "2", "4", "8", "16", "32", "64"]:
    scan.append({"z": z, "Delta_strich": float(Dp_direct(mp.mpf(z)))})
res["signature_und_entartung"] = {
    "Delta_strich_scan": scan,
    "Delta_strich_minimum_scan": float(min(s["Delta_strich"] for s in scan)),
    "eigenzahlen": {},
}
for a1v in ["0", "0.05", "0.15", "0.3", "0.45", "0.6", "1.0"]:
    M = Ginv_numeric("0.3", a1v, "1.7")
    ev = sorted(np.linalg.eigvalsh(M))
    res["signature_und_entartung"]["eigenzahlen"][f"alpha1={a1v}"] = {
        "eigenzahlen": [float(x) for x in ev],
        "negativ": int(sum(1 for x in ev if x < 0)),
        "positiv": int(sum(1 for x in ev if x > 0)),
        "betrag_kleinste": float(min(abs(x) for x in ev)),
        "det_direct": float(np.linalg.det(M)),
    }

mu_num_checks = []
for (a0v, a1v, betav) in [("0.3", "0.25", "1.7"), ("-0.4", "0.1", "1.0"), ("0.8", "-0.35", "2.5"), ("0.0", "0.6", "1.3")]:
    M = Ginv_numeric(a0v, a1v, betav)
    mu_direct = abs(float(np.linalg.det(M))) ** (-0.5)
    mu_sym = float(mu.subs({a0: a0v, a1: a1v, beta: betav}))
    mu_num_checks.append(
        {
            "alpha_0": a0v,
            "alpha_1": a1v,
            "beta": betav,
            "mu_direkt": mu_direct,
            "mu_formel": mu_sym,
            "relativer_rest": abs(mu_direct - mu_sym) / abs(mu_sym),
        }
    )
res["massdichte"]["numerische_kontrolle"] = mu_num_checks
dump()

print(json.dumps(res["anchor_homogendifferenz"], ensure_ascii=False))
print(json.dumps(res["G"], ensure_ascii=False)[:1600])
print(json.dumps(res["determinante"], ensure_ascii=False)[:900])
print(json.dumps(res["massdichte"], ensure_ascii=False)[:1500])
print(json.dumps(res["laplace_koeffizienten"], ensure_ascii=False)[:900])
print(json.dumps(res["signature_und_entartung"], ensure_ascii=False)[:1800])
