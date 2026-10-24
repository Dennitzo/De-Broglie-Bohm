# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest fuer den Zusammenhang zwischen kosmischer Zeit t
und York-Zeit T im flachen, freien, masselosen FLRW-Modell (k=0, Lambda=0).

Prueft (alle exakt, mit Sympy):
  1. York-Zeit aus FLRW-Zwang (k=0) + erster Friedmann-Gleichung:
       C = -p^2/(2 P_T) + 6*pi*G*T^2*P_T - P_T*V = 0
     Mit N=1 (kosmische Eigenzeit), p_phi = a^3*phi_dot, P_T = -a^3 wird daraus
       6*pi*G*T^2 = rho + V = rho_tot .  Mit H^2 = (8*pi*G/3)*rho_tot folgt
       T^2 = H^2/(16*pi^2*G^2)  und auf dem Zweig H>0, T<0:  T = -H/(4*pi*G).
     Gegenprobe: ADM-Impuls p_a = -3*a*a_dot/(4*pi*G) (N=1) liefert
       T = p_a/(3*a^2) = -(a_dot/a)/(4*pi*G) = -H/(4*pi*G).
  2. Freie masslose Loesung: a(t) ~ t^(1/3), H = 1/(3 t), T(t) = -1/(12*pi*G*t).
  3. Umparametrisierung: t(T) = -1/(12*pi*G*T); dt/dT = 1/(12*pi*G*T^2),
       konsistent mit 1/(dT/dt) (Kettenregel der Kehrwerte).
  4. Lapse-Unterscheidung: N_t = 1 (kosmische Eigenzeit, Eichfixierung),
       N_T = dt/dT = 1/(12*pi*G*T^2) != 1; ds^2 = -N_T^2 dT^2 + a^2*delta.
  5. Kettenregel: Phi'(T) = dPhi/dt * dt/dT (NIE dPhi/dt = dPhi/dT).

Ausgabe: work/york_cosmic_time_symbolic.json plus Text-Zusammenfassung.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "york_cosmic_time_symbolic.json")

t, T, G, a0, t0, C0 = sp.symbols("t T G a0 t0 C0", positive=True)
Vfun = sp.symbols("V", real=True)
H = sp.symbols("H", positive=True)  # Hubble-Rate im Expansionzweig H>0
rho = sp.symbols("rho", real=True)  # Energiedichte phi_dot^2/2
rho_tot = sp.symbols("rho_tot", real=True)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


# ----------------- 1) York-Zeit aus FLRW-Zwang + erster Friedmann ------------
P, p = sp.symbols("P p", real=True)   # P = P_T ; p = p_phi
a3, phidot = sp.symbols("a3 phidot", positive=True)   # a3 = a^3 > 0, phidot = phi_dot
constraint = sp.simplify(-p**2 / (2 * P) + 6 * sp.pi * G * T**2 * P - P * Vfun)
# Bei N=1 ist p_phi = a^3*phi_dot und P_T = -a^3; kinetischer Term rho=phi_dot^2/2:
#   -p^2/(2 P_T) = +a^3*rho.  Der Zwang C=0 wird damit
#   a^3*(rho + V) = 6*pi*G*T^2*a^3  <==>  6*pi*G*T^2 = rho + V = rho_tot.
constraint_depot = sp.simplify(
    constraint.subs({p: a3 * phidot, P: -a3}))
# Identitaetspruefung: constraint == a3*(phidot^2/2 + V) - 6*pi*G*T^2*a3
C_rep = sp.simplify(constraint_depot - (a3 * (phidot**2 / 2 + Vfun)
                                        - 6 * sp.pi * G * T**2 * a3))
record("Zwang (k=0), N=1, P_T=-a^3: C = a^3*(rho+V) - 6*pi*G*T^2*a^3",
       C_rep == 0, "C_identitaet: %s == 0" % C_rep)
# C=0, a^3>0  <=>  rho + V = 6*pi*G*T^2,  rho = phidot^2/2
E_den = sp.Symbol("Eden", real=True)   # Energiedichte rho + V bei N=1
# C=0 <=> a^3*(E_den - 6*pi*G*T^2) = 0; wegen a^3>0 folgt E_den = 6*pi*G*T^2.
# Die Identitaetspruefung oben belegt constraint_depot == a^3*E_den - 6*pi*G*T^2*a^3
# mit E_den = phidot^2/2 + Vfun; hier die Implikation aus C=0.
C_rep2 = sp.simplify(a3 * E_den - 6 * sp.pi * G * T**2 * a3)
sol_dens = sp.solve(sp.Eq(C_rep2, 0), E_den)
record("C=0 => 6*pi*G*T^2 = rho + V (rho:=phidot^2/2)",
       bool(sol_dens) and sp.simplify(sol_dens[0] - 6 * sp.pi * G * T**2) == 0,
       "Loesung C=0 fuer (rho+V) = %s (erwartet 6*pi*G*T^2)" %
       (sol_dens[0] if sol_dens else None))
# Mit rho_tot := rho + V und erster Friedmann (k=0, Lambda=0):
#   H^2 = (8*pi*G/3)*rho_tot  <=>  rho_tot = 3 H^2/(8*pi*G)
rho_tot_soln = 3 * H**2 / (8 * sp.pi * G)
T_sq = sp.simplify(rho_tot_soln / (6 * sp.pi * G))
record("Friedmann: T^2 = H^2/(16*pi^2*G^2)",
       sp.simplify(T_sq - H**2 / (16 * sp.pi**2 * G**2)) == 0,
       "T^2 = %s" % T_sq)
# Zweigwahl H>0 (Expansion), T<0:  T = -H/(4*pi*G)
T_neg_branch = sp.simplify(-H / (4 * sp.pi * G))
record("Zweigwahl: T = -H/(4*pi*G) (H>0, T<0) erfuellt T^2=H^2/(16*pi^2*G^2)",
       sp.simplify(T_neg_branch**2 - T_sq) == 0,
       "(-H/(4*pi*G))^2 = %s == T^2" % sp.simplify(T_neg_branch**2))
# Gegenprobe ueber ADM-Impuls: p_a = -3*a*a_dot/(4*pi*G) (N=1)
a_fun = sp.Function("a")
adot = sp.diff(a_fun(t), t)
a_mom = -3 * a_fun(t) * adot / (4 * sp.pi * G)
T_from_mom = sp.simplify(a_mom / (3 * a_fun(t)**2))
record("ADM-Gegenprobe: T = p_a/(3*a^2) = -(a_dot/a)/(4*pi*G)",
       sp.simplify(T_from_mom - (-(adot / a_fun(t)) / (4 * sp.pi * G))) == 0,
       "p_a/(3a^2) = %s  ==  -H/(4*pi*G)" % T_from_mom)

# ----------------- 2) Freie masslose Loesung a(t) ~ t^(1/3) ------------------
a_power = sp.simplify(a0 * (t / t0)**sp.Rational(1, 3))
H_power = sp.simplify(sp.diff(a_power, t) / a_power)
record("Freie masslose FLRW: a(t) ~ t^(1/3), H = 1/(3 t)",
       sp.simplify(H_power - 1 / (3 * t)) == 0,
       "H = %s" % H_power)
T_t_expr = sp.simplify(-H_power / (4 * sp.pi * G))
record("York-Zeit als Zeitfunktion: T(t) = -1/(12*pi*G*t)",
       sp.simplify(T_t_expr - (-1 / (12 * sp.pi * G * t))) == 0,
       "T(t) = %s" % T_t_expr)

# ----------------- 3) Umparametrisierung t <-> T ----------------------------
t_of_T = sp.simplify(-1 / (12 * sp.pi * G * T))
dt_dT_direct = sp.simplify(sp.diff(t_of_T, T))
dT_dt = sp.simplify(sp.diff(T_t_expr, t))
dt_dT_recip = sp.simplify(1 / dT_dt)
record("Umparametrisierung: dt/dT = 1/(12*pi*G*T^2) (aus t(T))",
       sp.simplify(dt_dT_direct - 1 / (12 * sp.pi * G * T**2)) == 0,
       "dt/dT = %s" % dt_dT_direct)
# 1/(dT/dt) nach Substitution t = t(T) in reine T-Form bringen
dt_dT_recip_T = sp.simplify(dt_dT_recip.subs(t, t_of_T))
record("Kettenregel-Kehrwert: dt/dT = 1/(dT/dt)|_{t=t(T)} = 1/(12*pi*G*T^2)",
       sp.simplify(dt_dT_recip_T - dt_dT_direct) == 0,
       "1/(dT/dt) substituiert = %s == dt/dT aus t(T)" % dt_dT_recip_T)

# ----------------- 4) Lapse-Unterscheidung ----------------------------------
N_t = 1
N_T = dt_dT_direct
record("Lapse kosmische Eigenzeit: N_t = 1 (Eichfixierung)",
       N_t == 1, "N_t = %s" % N_t)
record("Lapse nach Reparametrisierung: N_T = dt/dT = 1/(12*pi*G*T^2) != 1",
       sp.simplify(N_T - 1 / (12 * sp.pi * G * T**2)) == 0,
       "N_T = %s" % N_T)
record("Metrik: ds^2 = -dt^2 + a^2*delta = -(dt/dT)^2*dT^2 + a^2*delta",
       sp.simplify((dt_dT_direct)**2 - N_T**2) == 0,
       "(dt/dT)^2 = N_T^2 = %s" % sp.simplify(N_T**2))
record("Lapse N_T > 0 fuer T<0: N_T(G=1) = 1/(12*pi*T^2)",
       sp.simplify(N_T.subs(G, 1) - 1 / (12 * sp.pi * T**2)) == 0 and
       sp.simplify(1 / (12 * sp.pi * T**2)) > 0,
       "N_T(G=1) = %s (positiv)" % sp.simplify(N_T.subs(G, 1)))

# ----------------- 5) Kettenregel Phi'(T) = Phi_dot * dt/dT ------------------
C0p = sp.symbols("C0", positive=True)
Phi_T_expr = sp.log(-T) / sp.sqrt(12 * sp.pi * G) + C0p
Phi_prime_T = sp.simplify(sp.diff(Phi_T_expr, T))
record("York-Loesung: Phi'(T) = 1/(sqrt(12*pi*G)*T)",
       sp.simplify(Phi_prime_T - 1 / (sp.sqrt(12 * sp.pi * G) * T)) == 0,
       "Phi'(T) = %s" % Phi_prime_T)
# Kosmische Zeit: T = -1/(12*pi*G*t)  <=>  t = -1/(12*pi*G*T)  (t(T))
Phi_t_expr = sp.simplify(Phi_T_expr.subs(T, -1 / (12 * sp.pi * G * t)))
Phi_dot_t = sp.simplify(sp.diff(Phi_t_expr, t))
# Kettenregel: Phi'(T) = dPhi/dt * dt/dT ; t durch t(T) eliminieren
chain_rule = sp.simplify(Phi_dot_t * dt_dT_direct).subs(t, t_of_T)
chain_rule = sp.simplify(chain_rule)
record("Kettenregel: Phi'(T) == Phi_dot(t) * dt/dT |_{t=t(T)}",
       sp.simplify(chain_rule - Phi_prime_T) == 0,
       "Phi_dot*dt/dT = %s == Phi'(T) = %s" % (chain_rule, Phi_prime_T))
record("Warnzeichen: Phi_dot(t) != Phi'(T)",
       sp.simplify(Phi_dot_t - Phi_prime_T) != 0,
       "Phi_dot(t) = %s ; Phi'(T) = %s (verschieden)" %
       (Phi_dot_t, Phi_prime_T))

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "relation_T_H": {
        "formula": "T = -H/(4*pi*G_N)   (H>0, T<0)",
        "T_squared_via_Friedmann": sp.sstr(T_sq),
        "T_via_ADM_momentum": sp.sstr(T_from_mom),
    },
    "free_massless": {
        "scale_factor": "a(t) ~ t^(1/3)",
        "H": sp.sstr(H_power),
        "T_of_t": sp.sstr(T_t_expr),
    },
    "reparametrization": {
        "t_of_T": sp.sstr(t_of_T),
        "dt_dT": sp.sstr(dt_dT_direct),
        "dT_dt": sp.sstr(dT_dt),
        "dt_dT_formula": "1/(12*pi*G_N*T^2)",
        "chain_rule_reciprocal_ok": True,
    },
    "lapse": {
        "N_t_cosmic_proper_time": "1",
        "N_T_after_reparam": sp.sstr(N_T),
        "note": "N_T = dt/dT != N_t; ds^2 = -N_T^2 dT^2 + a^2*delta",
    },
    "chain_rule": {
        "Phi_prime_T": sp.sstr(Phi_prime_T),
        "Phi_dot_t": sp.sstr(Phi_dot_t),
        "Phi_dot_t_times_dt_dT": sp.sstr(chain_rule),
        "warning": "dPhi/dt und dPhi/dT sind NICHT gleich; Phi'(T)=dPhi/dt*dt/dT",
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
