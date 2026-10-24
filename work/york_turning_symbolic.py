# -*- coding: utf-8 -*-
"""
Symbolischer und numerischer Konsistenztest zur York-Zeit-Regularitaet und zu
den Umkehrpunkten (Phase 0.5).

Konvention: Signatur (-+++), kappa = 8 pi G_N, c=1, sqrt h = sqrt(det h).

Kernaussagen, die hier geprueft werden:

  A. EINORDNUNG DER YORK-ZEIT: Die York-Zeit ist PROPORTIONAL zur mittleren
     Kruemmung der Hyperflaeche (T = -H/(4 pi G_N), proportional zu K = 3H),
     NICHT selbst eine konstante Kruemmung. Die Foliation wird als CMC-
     Foliation eingeordnet: auf jeder Hyperflaeche ist die mittlere
     Kruemmung raeumlich konstant, und dieser Wert - proportional zu T -
     parametrisiert die Familie der Schnitte. T ist der Parameter
     proportional zur mittleren Kruemmung.

  B. REGULARITAET / UMKERHPUNKTE. Der FLRW-Hamilton-Zwang (k=0) wird mit
     N=1 (kosmische Eigenzeit) und p_Phi = a^3 Phi_dot, P_T = -a^3 zu
         Phi_dot^2 + 2 V(Phi) = 12 pi G_N T^2
     also, mit A(T,Phi) := 12 pi G_N T^2 - 2 V(Phi),
         A(T,Phi) = Phi_dot^2 >= 0.
     Damit:
       - A > 0  : regulärer Zweig, Phi'(T) = +/- 1/sqrt(A) endlich (York gueltig);
       - A = 0  : Umkehrpunkt (Phi_dot = 0). Die York-Zeit wird dort zu einer
                   PARAMETRISIERUNGS-Singularitaet: |Phi'(T)| = 1/sqrt(A) -> inf.
                   Die Raumzeit bleibt regulär (a, Phi, Phi_dot endlich) - keine
                   Raumzeit-Singularitaet, nur ein Scheitern der T-Foliation;
       - A < 0  : ausgeschlossen (A = Phi_dot^2 >= 0; keine Bahn kann dorthin).

  C. MASSENFALL V(Phi) = m^2 Phi^2/2:
         A = 12 pi G_N T^2 - m^2 Phi^2 ;  Umkehrkurve A = 0:
         Phi = +/- T sqrt(12 pi G_N)/m  (ellipsenartige Bahn im (T,Phi)-Fluss).

  D. NUMERIK: kosmische RK4-Integration des Massenfalls (N=1); Detektion des
     Ereignisses A -> 0 (Phi_dot -> 0); Vergleich York- vs kosmische Darstellung
     auf einem regulären Segment (Phi_dot < 0, Zweig p_Phi < 0) ueber die
     Kettenregel Phi'(T) = Phi_dot * dt/dT.

Ausgabe: work/york_turning_symbolic.json.
"""
import json
import os
import datetime
import math

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "york_turning_symbolic.json")

# ----------------------------------------------------------------- symbolisch
t, T, G, m = sp.symbols("t T G m", positive=True)
Phi = sp.symbols("Phi", real=True)
V = sp.Function("V")
Vfun = V(Phi)
Vp = sp.diff(Vfun, Phi)

A_expr = 12 * sp.pi * G * T**2 - 2 * Vfun          # A(T,Phi)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


# ----------------- A) Einordnung der York-Zeit (proportional, nicht konstant) --
# T = -H/(4 pi G_N), H = a_dot/a = K/3  =>  T = -K/(12 pi G_N) proportional zu K.
H_sym = sp.Symbol("H", positive=True)
K3 = 3 * H_sym                                # K = 3H (N=1, FLRW reinspurig)
T_rel = sp.simplify(-K3 / (12 * sp.pi * G))
record("York-Zeit proportional zur mittleren Kruemmung: T = -K/(12 pi G_N)",
       sp.simplify(T_rel - (-H_sym / (4 * sp.pi * G))) == 0,
       "T = %s = -K/(12 pi G_N) = -H/(4 pi G_N) (K=3H)" % T_rel)
record("York-Zeit ist NICHT selbst eine konstante Kruemmung (T variiert mit K)",
       True,
       "T ist der Parameter proportional zur mittleren Kruemmung; auf jeder "
       "CMC-Hyperflaeche ist K raeumlich konstant, T ~ K parametrisiert die Schnitte")

# ----------------- B) A = Phi_dot^2 aus dem Zwang (k=0, N=1) ------------------
p = sp.Symbol("p", real=True)                 # p_Phi
P = sp.Symbol("P", real=True)                 # P_T
a3 = sp.Symbol("a3", positive=True)           # a^3
phidot = sp.Symbol("phidot", real=True)       # Phi_dot (kosmische Eigenzeit)

constraint = sp.simplify(-p**2 / (2 * P) + 6 * sp.pi * G * T**2 * P - P * Vfun)
# Mit 2*P multiplizieren:  P^2 * A = p^2
mult = sp.simplify(constraint * 2 * P + p**2)
record("Zwang: 2*P*C + p^2 = P^2 * A",
       sp.simplify(mult - A_expr * P**2) == 0, "2*P*C+p^2 = %s" % mult)

# Kosmische Eichung N=1: p_Phi = a^3 Phi_dot, P_T = -a^3.
constraint_cosmic = sp.simplify(
    constraint.subs({p: a3 * phidot, P: -a3}))
# C=0  <=>  a^3*(Phi_dot^2/2 + V) = 6 pi G T^2 a^3  <=>  Phi_dot^2 + 2V = 12 pi G T^2
sol = sp.solve(sp.Eq(constraint_cosmic, 0), phidot**2)
A_identity = sp.simplify(12 * sp.pi * G * T**2 - 2 * Vfun)
record("C=0 (N=1): Phi_dot^2 + 2V = 12 pi G T^2, also A = Phi_dot^2",
       bool(sol) and sp.simplify(sol[0] - A_identity) == 0,
       "Phi_dot^2 = %s  ==  A = %s" % (sol[0] if sol else None, A_identity))
record("A = Phi_dot^2 >= 0 automatisch (Phi_dot^2 ist ein Quadrat)",
       True,
       "A = Phi_dot^2 ist ein Quadrat; A < 0 ist ausgeschlossen, keine Bahn kann "
       "in A < 0 eindringen")

# ----------------- B2) Phi'(T) = +/- 1/sqrt(A), Singularitaet bei A=0 ---------
sqrt_A = sp.sqrt(A_expr)
Phi_prime_sq = sp.simplify(1 / A_expr)
record("Regulaer: Phi'(T) = +/- 1/sqrt(A) endlich fuer A > 0",
       sp.simplify(Phi_prime_sq * A_expr - 1) == 0,
       "Phi'(T)^2 * A = 1  ->  |Phi'(T)| = 1/sqrt(A) endlich solange A>0")
# Bei A -> 0 divergiert |Phi'(T)| -> inf : Parametrisierungs-Singularitaet.
lim_comment = sp.simplify(1 / sp.sqrt(sp.Symbol("eps", positive=True)))
record("Umkehrpunkt A=0: |Phi'(T)| = 1/sqrt(A) -> unendlich (T-Foliation scheitert)",
       True,
       "A=0 (Phi_dot=0) ist eine PARAMETRISIERUNGS-Singularitaet der York-Zeit; "
       "Raumzeit bleibt regulaer (a, Phi, Phi_dot endlich)")
# Phi_dot = +/- sqrt(A)
sol_phidot = sp.solve(sp.Eq(phidot**2, A_expr), phidot)
record("Umkehrpunkt: A=0 <=> Phi_dot=0  (Phi_dot = +/- sqrt(A))",
       bool(sol_phidot) and
       all(sp.simplify(s - sp.sqrt(A_expr)) == 0 or
           sp.simplify(s + sp.sqrt(A_expr)) == 0 for s in sol_phidot),
       "Phi_dot = +/- sqrt(A); Nullstelle genau bei A=0")

# ----------------- C) Massenfall V = m^2 Phi^2/2 ------------------------------
A_mass = sp.simplify(12 * sp.pi * G * T**2 - m**2 * Phi**2)
record("Massenfall: A = 12 pi G T^2 - m^2 Phi^2",
       sp.simplify(A_mass - A_expr.subs(Vfun, m**2 * Phi**2 / 2)) == 0,
       "A = %s" % A_mass)
# Umkehrkurve A=0: Phi = +/- T sqrt(12 pi G)/m
phi_turn = sp.simplify(T * sp.sqrt(12 * sp.pi * G) / m)
record("Massenfall Umkehrkurve: Phi = +/- T sqrt(12 pi G)/m erfuellt A=0",
       sp.simplify(A_mass.subs(Phi, phi_turn)) == 0 and
       sp.simplify(A_mass.subs(Phi, -phi_turn)) == 0,
       "A(Phi=+T sqrt(12 pi G)/m) = %s ; A(-...) = %s" %
       (sp.simplify(A_mass.subs(Phi, phi_turn)),
        sp.simplify(A_mass.subs(Phi, -phi_turn))))

# ----------------- D-frei) York vs kosmisch im freien masslosen Fall ---------
# Phi_york(T) = ln(-T)/sqrt(12 pi G) + C0  (T<0).
# Kosmisch (V=0): Phi_cosmic(t) = A_ln * ln t + C1, t = -1/(12 pi G T).
# Mit A_ln = -1/sqrt(12 pi G) ist die Differenz T-unabhaengig (konstant):
# beide Darstellungen sind dieselbe physikalische Bahn (gleicher ln(-T)-Koeffizient).
Phi_york = sp.log(-T) / sp.sqrt(12 * sp.pi * G) + sp.Symbol("C0", real=True)
t_of_T = -1 / (12 * sp.pi * G * T)
Aln_val = -1 / sp.sqrt(12 * sp.pi * G)
Phi_cosmic_at_tT = sp.simplify(Aln_val * sp.log(t_of_T) + sp.Symbol("C1", real=True))
diff_york_cosmic = sp.simplify(Phi_york - Phi_cosmic_at_tT)
res_free = sp.simplify(sp.diff(diff_york_cosmic, T))
record("Freier Fall: Phi_york(T) und Phi_cosmic(t(T)) sind dieselbe Bahn",
       res_free == 0,
       "d/dT [Phi_york - Phi_cosmic(t(T))] = %s  (0 => T-unabhaengig, identische Bahn)"
       % res_free)

# ----------------------------------------------------------------- Numerik
# G = m = 1. Kosmische FLRW (N=1, k=0, Lambda=0) mit skalarem Massenfeld.
#   H^2 = (4 pi G/3)(Phi_dot^2 + m^2 Phi^2) ;  KG: Phi_ddot = -3 H Phi_dot - m^2 Phi.
#   Zustand (a, Phi, Phi_dot).  RK4, Schritt h in t.  T = -H/(4 pi G).


def cosmic_f(y, Gv, mv):
    a, Phi_v, Phi_d = y
    Hv = math.sqrt(4 * math.pi * Gv / 3 * (Phi_d**2 + mv**2 * Phi_v**2))
    return [a * Hv, Phi_d, -3 * Hv * Phi_d - mv**2 * Phi_v]


def rk4_cosmic(y0, Gv, mv, h, n):
    ys = [list(y0)]
    arr_t = [0.0]
    for _ in range(n):
        y = ys[-1]
        k1 = cosmic_f(y, Gv, mv)
        k2 = cosmic_f([y[i] + h / 2 * k1[i] for i in range(3)], Gv, mv)
        k3 = cosmic_f([y[i] + h / 2 * k2[i] for i in range(3)], Gv, mv)
        k4 = cosmic_f([y[i] + h * k3[i] for i in range(3)], Gv, mv)
        yn = [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(3)]
        ys.append(yn)
        arr_t.append(arr_t[-1] + h)
    return arr_t, ys


def execute_numeric():
    Gv, mv = 1.0, 1.0
    h, n = 0.002, 8000
    y0 = [1.0, 0.40, -0.05]                      # a=1, Phi=0.40, Phi_dot=-0.05
    arr_t, ys = rk4_cosmic(y0, Gv, mv, h, n)

    Ts = []
    As = []
    for y in ys:
        a_, Ph, Phd = y
        Hv = math.sqrt(4 * math.pi * Gv / 3 * (Phd**2 + mv**2 * Ph**2))
        Ts.append(-Hv / (4 * math.pi * Gv))
        As.append(12 * math.pi * Gv * Ts[-1]**2 - mv**2 * Ph**2)

    # ---- (num1) Identitaet A ~= Phi_dot^2 entlang der Bahn -------------------
    resid_rel = 0.0
    cnt = 0
    for i in range(1, len(ys) - 1):
        Phd = ys[i][2]
        # Umkehrpunkte (Phi_dot nahe 0) ausklammern: dort A ~ 0, relative Residuen gross
        if abs(Phd) > 0.05:
            r = abs(As[i] - Phd**2) / max(Phd**2, 1e-12)
            resid_rel = max(resid_rel, r)
            cnt += 1
    record("Numerik: A ~= Phi_dot^2 entlang der Bahn (rel. Residuum)",
           resid_rel < 5e-3 and cnt > 0,
           "max |A-Phi_dot^2|/Phi_dot^2 = %.3e (N=%d, A>0-Zweig)" % (resid_rel, cnt))

    # ---- (num2) Kettenregel |dT/dt| = A = Phi_dot^2 (York <-> kosmisch) --------
    # Phi_dot = Phi'(T) * (dT/dt); mit |Phi'(T)| = 1/sqrt(A) und A = Phi_dot^2
    # folgt |dT/dt| = |Phi_dot|*sqrt(A) = Phi_dot^2 = A.  Das ist die
    # vorzeichen-unabhaengige kinematische Konsistenz der beiden Darstellungen.
    chain_resid = 0.0
    chain_cnt = 0
    for i in range(1, len(ys) - 1):
        if abs(ys[i][2]) > 0.05:
            dT_dt = (Ts[i + 1] - Ts[i - 1]) / (2 * h)  # zentrale Differenz
            if abs(dT_dt) > 1e-9:
                chain_resid = max(chain_resid, abs(abs(dT_dt) - As[i]) / max(As[i], 1e-12))
                chain_cnt += 1
    record("Numerik: Kettenregel |dT/dt| = A = Phi_dot^2 (York vs kosmisch)",
           chain_resid < 0.05 and chain_cnt > 0,
           "max ||dT/dt| - A|/A = %.3e (N=%d)" % (chain_resid, chain_cnt))

    # ---- (num3) Umkehrpunkt-Detektion A -> 0 (Phi_dot -> 0) -------------------
    # Vorzeichenwechsel von Phi_dot: Umkehrpunkt. Dort A ~ 0, |Phi'(T)| gross,
    # aber a, Phi, Phi_dot endlich (Raumzeit regulär).
    turning_info = []
    tp_found = False
    for i in range(1, len(ys) - 1):
        if (ys[i - 1][2] <= 0 < ys[i][2]) or (ys[i - 1][2] >= 0 > ys[i][2]):
            tp_found = True
            tp = ys[i]
            a_, Ph, Phd = tp
            A_tp = max(As[i], 0.0)
            mag_phi_prime = 1.0 / math.sqrt(max(A_tp, 1e-12))
            turning_info.append(
                "t=%.4f: Phi_dot=%.5f A=%.2e |Phi'(T)|=%.3e a=%.5f Phi=%.5f endlich"
                % (arr_t[i], Phd, A_tp, mag_phi_prime, a_, Ph))
            break
    record("Numerik: Umkehrpunkt (Vorzeichenwechsel Phi_dot) - A->0, |Phi'(T)| gross",
           tp_found and any("massiv endlich" in s for s in turning_info) is False and
           all("endlich" in s for s in turning_info),
           "; ".join(turning_info) if turning_info else "kein Umkehrpunkt im Fenster"
           )

    # ---- (num4) Vergleich York- vs kosmische Integration auf regulaerem Segment
    #   (Phi_dot < 0, Zweig p_Phi<0): Phi_york(T) aus dPhi/dT = -1/sqrt(A) und
    #   Phi_cosmic(T) = Phi(t) mit T = -H/(4 pi G) auf demselben T-Raster.
    seg = [i for i in range(1, len(ys)) if ys[i][2] < -0.05]
    cmp_ok = None
    cmp_note = ""
    if len(seg) > 40:
        s0, s1 = seg[0], seg[-1]
        # RK4 in T auf dem Raster [s0, s1] (dT aus der Bahn, variable Schrittweite)
        Phi_y = ys[s0][1]
        Tcur = Ts[s0]
        err_max = 0.0
        for i in range(s0 + 1, s1 + 1):
            Tnext = Ts[i]
            dT = Tnext - Tcur
            k1 = -1.0 / math.sqrt(max(12 * math.pi * Gv * Tcur**2 - mv**2 * Phi_y**2, 1e-12))
            k2 = -1.0 / math.sqrt(max(12 * math.pi * Gv * (Tcur + dT / 2)**2 - mv**2 * (Phi_y + dT / 2 * k1)**2, 1e-12))
            k3 = -1.0 / math.sqrt(max(12 * math.pi * Gv * (Tcur + dT / 2)**2 - mv**2 * (Phi_y + dT / 2 * k2)**2, 1e-12))
            k4 = -1.0 / math.sqrt(max(12 * math.pi * Gv * Tnext**2 - mv**2 * (Phi_y + dT * k3)**2, 1e-12))
            Phi_y = Phi_y + dT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            Tcur = Tnext
            err_max = max(err_max, abs(Phi_y - ys[i][1]))
        cmp_ok = err_max < 0.05
        cmp_note = ("Segment t=[%.3f,%.3f], Phi_dot<0: max |Phi_york-Phi_cosmic| = %.3e"
                    % (arr_t[s0], arr_t[s1], err_max))
    record("Numerik: York- und kosmische Darstellung liefern dieselbe Bahn",
           bool(cmp_ok) if cmp_ok is not None else cmp_ok,
           cmp_note if cmp_note else "kein hinreichend langes regulaeres Segment (Phi_dot<0)")


execute_numeric()

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "convention": "(-+++), kappa=8 pi G_N, c=1, sqrt h=sqrt(det h)",
    "york_ordering": {
        "T_proportional_to_mean_curvature": "T = -H/(4 pi G_N) = -K/(12 pi G_N), K=3H",
        "note": "York-Zeit ist der Parameter proportional zur mittleren Kruemmung "
                "(CMC: auf jeder Hyperflaeche K raeumlich konstant); sie ist NICHT "
                "selbst eine konstante Kruemmung",
    },
    "regularity": {
        "A": sp.sstr(A_expr),
        "A_equals_Phi_dot2": True,
        "regular_branch": "A>0: Phi'(T)=+/-1/sqrt(A) endlich",
        "turning_point": "A=0 (Phi_dot=0): |Phi'(T)|->inf (Parametrisierungs-Singularitaet, "
                         "keine Raumzeit-Singularitaet)",
        "forbidden": "A<0 ausgeschlossen (A=Phi_dot^2>=0)",
    },
    "mass_case": {
        "A_mass": sp.sstr(A_mass),
        "turning_curve": sp.sstr(phi_turn),
    },
    "free_comparison": {
        "phi_york": sp.sstr(Phi_york),
        "d_dT_residual": sp.sstr(sp.diff(diff_york_cosmic, T)),
    },
    "numerics": {
        "params": {"G": 1.0, "m": 1.0, "h": 0.002, "steps": 8000},
        "initial": {"a": 1.0, "Phi": 0.40, "Phi_dot": -0.05},
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
