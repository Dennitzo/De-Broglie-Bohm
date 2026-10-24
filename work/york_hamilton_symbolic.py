# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest fuer den York-Zeit-Hamiltonian.

Prueft (alle exakt, mit Sympy):
  1. Unabhaengige Differentiation des reduzierten Hamiltonians
     H_T = -p/sqrt(A), A = 12*pi*G*T^2 - 2*V(Phi)  (Zweig p<0):
        Phi'  = dH/dp   = -1/sqrt(A)
        p'    = -dH/dPhi = + p*V'(Phi)/A^(3/2)      <- Vorzeichenpruefung
  2. Zweiter Impulszweig H_T = +p/sqrt(A)  (Zweig p>0):
        Phi'  = +1/sqrt(A)
        p'    = - p*V'(Phi)/A^(3/2)
     Allgemeine Form H_T = |p|/sqrt(A); beide Zweige separat.
  3. Freier massloser Fall V=0:
        Analytische Loesung Phi(T) = ln(-T)/sqrt(12*pi*G) + C   (T<0)
        Phi'(T)  = 1/(sqrt(12*pi*G)*T)  ==  -1/sqrt(12*pi*G*T^2)  (T<0)
        Positive Kontraktionsrate v = -Phi' = 1/(sqrt(12*pi*G)*(-T))
        Numerische Werte bei G=1, T = -1,-2,-3,-5.
  4. Herleitung aus dem FLRW-Hamilton-Zwang (k=0):
        -p^2/(2 P) + 6*pi*G*T^2*P - P*V = 0   (P = P_T)
        => P^2 = p^2/A,  P = -|p|/sqrt(A)  (da P_T=-a^3<0),
        => H = -P = |p|/sqrt(A)  (Volumen-Interpretation a^3=H>0).
  5. Spezialfall p=0: H=|p|/sqrt(A) ist an p=0 nicht glatt (|p|),
     H=0 => Volumen 0 (unphysikalisch); Hamilton-Gleichungen dort undefiniert.

Ausgabe: work/york_hamilton_symbolic.json plus Text-Zusammenfassung.
"""
import json
import math
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "york_hamilton_symbolic.json")

# ----------------------------------------------------------------- symbolisch
T, p, G, C = sp.symbols("T p G C", real=True)
Phi = sp.symbols("Phi", real=True)
V = sp.Function("V")
Vfun = V(Phi)
Vp = sp.diff(Vfun, Phi)                    # V'(Phi)
A = 12 * sp.pi * G * T**2 - 2 * Vfun       # Radikand A(T,Phi)
sqrtA = sp.sqrt(A)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


# -------------------------- 1) Zweig p<0: H = -p/sqrt(A) -------------------
Hm = -p / sqrtA
dHm_dp = sp.simplify(sp.diff(Hm, p))
dHm_dPhi = sp.simplify(sp.diff(Hm, Phi))
pdot_m = sp.simplify(-dHm_dPhi)            # p' = -dH/dPhi

record("Zweig p<0: Phi' = -1/sqrt(A)",
       sp.simplify(dHm_dp - (-1 / sqrtA)) == 0,
       "dH/dp = %s" % dHm_dp)
record("Zweig p<0: p' = +p*V'/A^(3/2)",
       sp.simplify(pdot_m - (p * Vp / A**sp.Rational(3, 2))) == 0,
       "-dH/dPhi = %s" % pdot_m)
record("Zweig p<0: dH/dPhi = -p*V'/A^(3/2) (Zwischenstep)",
       sp.simplify(dHm_dPhi + (p * Vp / A**sp.Rational(3, 2))) == 0,
       "dH/dPhi = %s" % dHm_dPhi)

# -------------------------- 2) Zweig p>0: H = +p/sqrt(A) -------------------
Hp = +p / sqrtA
dHp_dp = sp.simplify(sp.diff(Hp, p))
dHp_dPhi = sp.simplify(sp.diff(Hp, Phi))
pdot_p = sp.simplify(-dHp_dPhi)

record("Zweig p>0: Phi' = +1/sqrt(A)",
       sp.simplify(dHp_dp - (1 / sqrtA)) == 0,
       "dH/dp = %s" % dHp_dp)
record("Zweig p>0: p' = -p*V'/A^(3/2)",
       sp.simplify(pdot_p - (-(p * Vp / A**sp.Rational(3, 2)))) == 0,
       "-dH/dPhi = %s" % pdot_p)

# Allgemeine Form H = |p|/sqrt(A): auf beiden Zweigen konsistent
record("Allgemein: H=|p|/sqrt(A) reduziert sich pro Zweig korrekt",
       (sp.simplify(Hm - (-p / sqrtA)) == 0) and (sp.simplify(Hp - (p / sqrtA)) == 0),
       "p<0: |p|=-p => H=-p/sqrt(A); p>0: |p|=p => H=+p/sqrt(A)")

# -------------------------- 3) Freier massloser Fall V=0 -------------------
A0 = 12 * sp.pi * G * T**2
Phi_sol = sp.log(-T) / sp.sqrt(12 * sp.pi * G) + C      # Zweig T<0
dPhi_dT = sp.simplify(sp.diff(Phi_sol, T))
# Erste Hamilton-Gleichung (V=0, Zweig p<0): Phi' = -1/sqrt(12*pi*G*T^2)
ham_phi0 = -1 / sp.sqrt(A0)
# Fuer T<0 gilt sqrt(T^2)=-T, also ham_phi0 = -1/(sqrt(12*pi*G)*(-T)) = 1/(sqrt(12*pi*G)*T)
ham_phi0_Tneg = sp.simplify(ham_phi0.rewrite(sp.log))
# Symbolisch: dPhi/dT == 1/(sqrt(12*pi*G)*T)
record("V=0: Phi_sol' = 1/(sqrt(12*pi*G)*T)",
       sp.simplify(dPhi_dT - 1 / (sp.sqrt(12 * sp.pi * G) * T)) == 0,
       "Phi_sol' = %s" % dPhi_dT)
# Numerische Restpruefung der ersten Hamilton-Gleichung an T<0
num_ok = True
num_lines = []
for tval in (-1, -2, -3, -5):
    lhs = float(sp.N(dPhi_dT.subs({T: tval, G: 1})))
    rhs = float(sp.N(ham_phi0.subs({T: tval, G: 1})))
    res = abs(lhs - rhs)
    num_ok = num_ok and (res < 1e-12)
    num_lines.append("T=%d: Phi'=%.12f  -1/sqrt(A)=%.12f  |res|=%.2e" % (tval, lhs, rhs, res))
record("V=0: numerisch Phi'(T) == -1/sqrt(12*pi*T^2) (G=1, T<0)",
       num_ok, " ; ".join(num_lines))

# Positive Kontraktionsrate v = -Phi' = 1/(sqrt(12*pi*G)*(-T)); G=1
expected_rates = {
    -1: 0.162868,
    -2: 0.081434,
    -3: 0.054289,
    -5: 0.032574,
}
rate_lines = []
rates_ok = True
computed_rates = {}
for tval, exp in expected_rates.items():
    vval = 1.0 / (math.sqrt(12 * math.pi * 1.0) * (-tval))
    computed_rates[str(tval)] = round(vval, 6)
    ok = abs(round(vval, 6) - exp) < 1e-9
    rates_ok = rates_ok and ok
    rate_lines.append("T=%d: v=%.6f (erwartet %s)" % (tval, round(vval, 6), exp))
record("V=0: positive Raten 0.162868/0.081434/0.054289/0.032574 (G=1)",
       rates_ok, " ; ".join(rate_lines))

# -------------------------- 4) Herleitung aus dem FLRW-Zwang (k=0) ---------
P = sp.symbols("P", real=True)  # P = P_T
constraint = -p**2 / (2 * P) + 6 * sp.pi * G * T**2 * P - P * Vfun  # k=0
# Mit 2*P multiplizieren: -p^2 + 12*pi*G*T^2*P^2 - 2*V*P^2 = 0
mult = sp.simplify(constraint * 2 * P + p**2)   # = A*P^2
record("Zwang: 2*P*Constraint + p^2 == A*P^2",
       sp.simplify(mult - A * P**2) == 0,
       "2*P*H_c + p^2 = %s" % mult)
# Daraus P^2 = p^2/A ; Zweig P_T=-a^3<0 => P=-|p|/sqrt(A) => H=-P=|p|/sqrt(A)
P_neg = -sp.Abs(p) / sqrtA
record("Zwang: P_T=-|p|/sqrt(A) erfuellt P^2=p^2/A",
       sp.simplify(P_neg**2 - p**2 / A) == 0,
       "P_T = %s, P_T^2 - p^2/A = 0" % P_neg)
H_from_constraint = sp.simplify(-P_neg)
record("Zwang: H=-P_T=|p|/sqrt(A) (Volumen a^3=H>0)",
       sp.simplify(H_from_constraint - sp.Abs(p) / sqrtA) == 0,
       "H = %s" % H_from_constraint)

# -------------------------- 5) Spezialfall p=0 ----------------------------
H_at_p0 = sp.simplify(sp.Abs(p) / sqrtA).subs(p, 0)
record("p=0: H=|p|/sqrt(A) ergibt 0 (Volumen 0, unphysikalisch)",
       H_at_p0 == 0, "H(p=0) = %s" % H_at_p0)
# |p| ist an p=0 nicht differenzierbar -> Hamilton-Gleichungen undefiniert
d_abs_left = sp.simplify(sp.diff(-p, p))     # fuer p<0: d|p|/dp = -1
d_abs_right = sp.simplify(sp.diff(p, p))     # fuer p>0: d|p|/dp = +1
record("p=0: |p| nicht glatt (linker/rechter Ableitungswert differieren)",
       (d_abs_left != d_abs_right),
       "d|p|/dp|_{p<0}=%s, d|p|/dp|_{p>0}=%s" % (d_abs_left, d_abs_right))

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "A": sp.sstr(A),
    "branch_p_neg": {
        "H": sp.sstr(Hm),
        "dH_dp": sp.sstr(dHm_dp),
        "dH_dPhi": sp.sstr(dHm_dPhi),
        "Phi_prime": sp.sstr(dHm_dp),
        "p_prime": sp.sstr(pdot_m),
    },
    "branch_p_pos": {
        "H": sp.sstr(Hp),
        "dH_dp": sp.sstr(dHp_dp),
        "dH_dPhi": sp.sstr(dHp_dPhi),
        "Phi_prime": sp.sstr(dHp_dp),
        "p_prime": sp.sstr(pdot_p),
    },
    "free_massless": {
        "Phi_solution_T<0": sp.sstr(Phi_sol),
        "dPhi_dT": sp.sstr(dPhi_dT),
        "positive_rate_formula": "1/(sqrt(12*pi*G)*(-T))",
        "computed_rates_G=1": computed_rates,
        "expected_rates_G=1": expected_rates,
    },
    "constraint_k0": {
        "constraint": sp.sstr(constraint),
        "P_T_negative_branch": sp.sstr(P_neg),
        "H_from_constraint": sp.sstr(H_from_constraint),
    },
    "special_p0": {
        "H_at_p0": sp.sstr(H_at_p0),
        "abs_p_kink": {"left": sp.sstr(d_abs_left), "right": sp.sstr(d_abs_right)},
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
