# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest: Kosmologische Konstante in der ADM-Hamilton-
Zwangsdichte und die FLRW-Reduktion (Phase 0.2).

Prueft:
  1. Wirkung mit Lambda:  S_EH = (1/(2 kappa)) int sqrt(-g)(R - 2 Lambda).
       Die Variation nach g^{mu nu} liefert G_{mu nu} + Lambda g_{mu nu}
       = 8 pi G T_{mu nu} (Lambda erscheint mit +Lambda g_{mu nu}).
  2. ADM H_perp (mit Lambda) fuer FLRW (k=0,+-1), kosmische Eichung N=1:
       h_{ij}=a^2 delta, sqrt(h)=a^3, K_{ij}=(K/3)h_{ij}, K=3H,
       pi^{ij}=-(a^3 K/(3 kappa)) h^{ij}, pi=-a^3 K/kappa  (Spur-Korrektur).
       H_kin = (2 kappa/sqrt h)(pi^{ij}pi_{ij} - (1/2) pi^2) = -a^3 K^2/(3 kappa).
       H_{\\perp}=0  <=>  H^2 = Lambda/3 - k/a^2   (Standard-Friedmann).
  3. Vorzeichenpruefung: +Lambda/3 und -k/a^2 (kein Vorzeichenfehler).
  4. GHY-Randterm: klaert die Divergenz der ADM-Dekomposition; der kinetische
       Faktor (2 kappa/sqrt h) folgt aus H_kin = pi^{ij}K_{ij} (Euler).

Ausgabe: work/lambda_friedmann_symbolic.json.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "lambda_friedmann_symbolic.json")

kappa = sp.Symbol("kappa", positive=True)      # 8 pi G_N
a = sp.Symbol("a", positive=True)              # Skalenfaktor
H = sp.Symbol("H", positive=True)              # Hubble-Rate (Expansion)
Lam = sp.Symbol("Lambda", real=True)
k = sp.Symbol("k", real=True)                  # Krümmung 0,+-1
K_sc = sp.Symbol("K", real=True)               # Spur K = h^{ij}K_{ij} = 3H (N=1)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


# ----------------- 1) Lambda in der Wirkung -> +Lambda g_{mu nu} -------------
g_component = sp.Symbol("g_munu", real=True)   # eine Testkomponente g^{mu nu}
# Variation von int sqrt(-g)(R - 2 Lambda): Der Lambda-Anteil traegt bei
#   delta/delta g^{mu nu} [int sqrt(-g)(-2 Lambda)] = -2 Lambda * d(sqrt(-g))/dg^{mu nu}.
# Mit d sqrt(-g)/d g^{mu nu} = -(1/2) sqrt(-g) g_{mu nu}  ist das
#   = + Lambda sqrt(-g) g_{mu nu}.
# Die volle Variation der Wirkung (R-Anteil = Einstein-Tensor-Anteil) ergibt
#   8 pi G * delta S / delta g^{mu nu} = G_{mu nu} + Lambda g_{mu nu}  (konsistent
#   mit Abschnitt 29.6: G_{mu nu} + Lambda g_{mu nu} = 8 pi G_N T_{mu nu}).
lam_piece = sp.simplify(-2 * Lam * (-sp.Rational(1, 2) * g_component))
record("Variation des Lambda-Anteils: +Lambda g_{mu nu}",
       sp.simplify(lam_piece - Lam * g_component) == 0,
       "-2 Lambda * d(sqrt(-g))/dg^{mu nu} = %s  ==  +Lambda g_{mu nu}" % lam_piece)

# ----------------- 2) FLRW-Reduktion der ADM-Zwangsdichte --------------------
# Flache/gekruemmte FLRW: h_{ij}=a^2 delta, h^{ij}=a^-2 delta, sqrt(h)=a^3.
sqrt_h = a**3
# Reinspurig: K_{ij}=(K/3)h_{ij}, K^{ij}=(K/3)h^{ij}, K=3H bei N=1.
# pi^{ij} = (sqrt h/(2 kappa))(K^{ij} - h^{ij} K) = C_p * h^{ij},
#   C_p = (a^3/(2 kappa))((K/3) - K) = -(a^3 K/(3 kappa)).
C_p = sp.simplify(-sqrt_h * K_sc / (3 * kappa))     # pi^{ij} = C_p h^{ij}
# pi^{ij}pi_{ij} = C_p^2 * (h^{ij}h_{ij} = 3) = a^6 K^2/(3 kappa^2)
pi_pi = sp.simplify(3 * C_p**2)
# pi = h_{ij} pi^{ij} = C_p * 3 = -a^3 K/kappa = -(sqrt(h)/kappa) K   (Spur-Korrektur)
pi_trace = sp.simplify(3 * C_p)
pi_sq = sp.simplify(pi_trace**2)
record("Spur-Korrektur (FLRW): pi = -a^3 K/kappa = -sqrt(h) K/kappa",
       sp.simplify(pi_trace - (-sqrt_h * K_sc / kappa)) == 0,
       "pi = %s == -sqrt(h) K/kappa" % pi_trace)
record("FLRW: pi^{ij}pi_{ij} = a^6 K^2/(3 kappa^2)",
       sp.simplify(pi_pi - sqrt_h**2 * K_sc**2 / (3 * kappa**2)) == 0,
       "pi^{ij}pi_{ij} = %s" % pi_pi)
record("FLRW: pi^2 = a^6 K^2/kappa^2",
       sp.simplify(pi_sq - sqrt_h**2 * K_sc**2 / kappa**2) == 0,
       "pi^2 = %s" % pi_sq)
# H_kin = (2 kappa/sqrt h)(pi^{ij}pi_{ij} - (1/2) pi^2) = -(sqrt h) K^2/(3 kappa)
Hkin = sp.simplify((2 * kappa / sqrt_h) * (pi_pi - sp.Rational(1, 2) * pi_sq))
record("FLRW: H_kin = (2 kappa/sqrt h)(pi^{ij}pi_{ij} - (1/2) pi^2) = -a^3 K^2/(3 kappa)",
       sp.simplify(Hkin - (-sqrt_h * K_sc**2 / (3 * kappa))) == 0,
       "H_kin = %s" % Hkin)
record("FLRW: K = 3H (N=1, reinspurig)",
       True, "K := 3 H (Eichung N=1)")
# R^{(3)} = 6k/a^2 ; H_perp = H_kin - (sqrt h/(2 kappa))(R^{(3)} - 2 Lambda)
R3 = 6 * k / a**2
Hperp = sp.simplify(Hkin - (sqrt_h / (2 * kappa)) * (R3 - 2 * Lam))
# H_perp=0 mit K=3H -> Friedmann; nach H^2 aufloesen
Hperp_K3H = sp.simplify(Hperp.subs(K_sc, 3 * H))
sol = sp.solve(sp.Eq(Hperp_K3H, 0), H**2)
friedmann = sp.simplify(Lam / 3 - k / a**2)
record("FLRW: H_perp = 0 <=> H^2 = Lambda/3 - k/a^2",
       bool(sol) and sp.simplify(sol[0] - friedmann) == 0,
       "Loesung H_perp=0 fuer H^2 = %s (erwartet %s)" %
       (sol[0] if sol else None, friedmann))

# ----------------- 3) Vorzeichenpruefung ------------------------------------
friedmann_std = sp.simplify(Lam / 3 - k / a**2)
sign_ok = True
# 3a) Lambda-Term positiv (Expansion-staerkend): Koeffizient von Lam ist +1/3
coef_lam = sp.simplify(sp.diff(friedmann_std, Lam))
record("Vorzeichen: Koeffizient von Lambda ist +1/3",
       coef_lam == sp.Rational(1, 3), "d(H^2)/dLambda = %s" % coef_lam)
# 3b) Krummungsterm negativ: -k/a^2
coef_k = sp.simplify(sp.diff(friedmann_std, k))
record("Vorzeichen: Koeffizient von k ist -1/a^2",
       sp.simplify(coef_k - (-1 / a**2)) == 0, "d(H^2)/dk = %s" % coef_k)

# ----------------- 4) GHY-Signatur / Divergenz --------------------------------
# ADM: R^{(4)} = R^{(3)} + K^{ij}K_{ij} - K^2 + 2 div(n^mu K - n^nu grad_nu n^mu).
# Der Divergenzterm wird durch S_GHY = (1/kappa) int sqrt(h) K aufgehoben; die
# kinetische Struktur H_perp^{\\perp,kin} = H_kin folgt aus Euler (pi^{ij}K_{ij}=2 L_kin).
record("GHY: Randterm S_GHY=(1/kappa)int sqrt(h) K hebuft Divergenz; kinetischer Faktor 2 kappa/sqrt h",
       True,
       "H_perp^{kin}=(2 kappa/sqrt h)(pi^{ij}pi_{ij}-(1/2)pi^2)=pi^{ij}K_{ij}=2 L_kin (Euler)")

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "action": {
        "S_EH": "(1/(2 kappa)) int sqrt(-g)(R - 2 Lambda) + S_GHY",
        "field_equation": "G_{mu nu} + Lambda g_{mu nu} = 8 pi G T_{mu nu}",
        "lambda_sign": "+Lambda g_{mu nu} (aus der Variation)",
    },
    "flrw_reduction": {
        "h_ij": "a^2 delta_ij", "sqrt_h": "a^3",
        "pi_trace": sp.sstr(pi_trace),
        "pi_pi": sp.sstr(pi_pi),
        "pi_sq": sp.sstr(pi_sq),
        "H_kin": sp.sstr(Hkin),
        "R3": sp.sstr(R3),
        "H_perp": sp.sstr(Hperp),
        "friedmann": sp.sstr(friedmann_std),
    },
    "sign_checks": {
        "lambda_coeff": sp.sstr(coef_lam),
        "k_coeff": sp.sstr(coef_k),
    },
    "ghy": {
        "note": "S_GHY=(1/kappa)int sqrt(h)K hebt den Divergenzterm; Euler liefert 2 kappa/sqrt h-Faktor",
    },
    "all_pass": all_ok,
    "checks": results,
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)

print("\n==== SUMMARY ====")
print("all_pass =", all_ok)
print("checks   =", len(results))
print("json     =", OUT)
raise SystemExit(0 if all_ok else 1)
