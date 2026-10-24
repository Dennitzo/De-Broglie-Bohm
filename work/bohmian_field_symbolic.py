# -*- coding: utf-8 -*-
"""
Phase 2 (Teil 1) - Herleitung der Bohmschen Bewegung aus der Schroedinger-
Gleichung fuer ein endlich dimensionales (modenzerlegtes) Quantenfeld.

Das modenzersierte freie Skalarfeld ist ein System von Oszillatoren; die
Feldkonfiguration ist q = (q_1..q_n). Behandelte Schr odinger-Gleichung:

  i hbar d_t psi = H psi,     H = -sum_i (hbar^2/(2 m_i)) d_i^2 + V(q).

Ohne Fuehrungsannahme wird psi = R exp(i S/hbar) eingesetzt; dann ist die
Schr odinger-Gleichung exakt aequivalent zu

  Kontinuitaet:   d_t(R^2) + sum_i d_i( R^2 (d_i S)/m_i ) = 0,
  Quanten-HJ:     d_t S + sum_i (d_i S)^2/(2 m_i) + V + Q = 0,
                  Q := -sum_i (hbar^2/(2 m_i)) (d_i^2 R)/R.

Die Fuehrungsgeschwindigkeit wird strombasiert ausgewertet,
  v_i = j_i/rho,   j_i = hbar Im(psi* d_i psi)/m_i,   rho = |psi|^2,
was mit d_i S/m_i identisch ist (P1).

Pruefungen:
  S1  SE  <=>  (Kontinuitaet + Quanten-HJ) fuer n=2 mit allgemeinen R,S,V.
  P1  j_i/rho = d_i S/m_i  (Stromform == Phasenform).
  S2  Kontinuitaet in rho-Form ist dieselbe Erhaltungsaussage wie in R-Form.
  S3  Newton-Form m(d_t + v.grad)v_i = -d_i(V+Q) exakt aus d_i(Quanten-HJ).
  S4  HO-Grundzustand: Q = E0 - V, also V+Q = E0 konstant, v = 0 (stationaere Bahn).
  S5  Superposition psi_0+psi_1: Dichteknoten, Q divergiert am Knoten,
      Fuehrungsgeschwindigkeit nicht identisch Null und zeitabhaengig.
  S6  Reelle Wellenfunktion => v identisch 0 (relative Phaasen noetig).
  S7  Produktzustand: d v_1/d q_2 = 0; verschaenkter Zustand (|00>+|11>)/sqrt2:
      d v_1/d q_2 != 0 (Nichtseparabilitaet der Fuhrung).

Ausgabe: work/bohmian_field_symbolic.json.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "bohmian_field_symbolic.json")

hbar = sp.Symbol("hbar", positive=True)
results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def rho_of(psf):
    return sp.simplify(sp.expand(psf * sp.conjugate(psf)))


def v_of(psf, coord, mass):
    """Strombasierte Fuehrungsgeschwindigkeit j/rho (reell)."""
    cur = sp.simplify((sp.conjugate(psf) * sp.diff(psf, coord)
                       - psf * sp.diff(sp.conjugate(psf), coord)) / (2 * sp.I))
    return sp.simplify(hbar * cur / (mass * rho_of(psf)))


def R_of(psf):
    return sp.sqrt(rho_of(psf))


def is_real(expr):
    """Kriterium: Ausdruck ist (nach Zusammenfassung) reell."""
    return sp.simplify(sp.im(sp.expand_complex(sp.together(expr)))) == 0


# =====================================================================  S1 ====
q1, q2, tt = sp.symbols("q1 q2 t", real=True)
m1, m2 = sp.Symbol("m1", positive=True), sp.Symbol("m2", positive=True)
R = sp.Function("R", real=True)(q1, q2, tt)
S = sp.Function("S", real=True)(q1, q2, tt)
Vg = sp.Function("V", real=True)(q1, q2, tt)

masses = {q1: m1, q2: m2}
coords = [q1, q2]
psi = R * sp.exp(sp.I * S / hbar)

Hpsi = sum(-hbar**2 / (2 * masses[c]) * sp.diff(psi, c, 2) for c in coords) + Vg * psi
resid = sp.simplify((sp.I * hbar * sp.diff(psi, tt) - Hpsi) / psi)

Qexpr = sp.simplify(-sum(hbar**2 / (2 * masses[c]) * sp.diff(R, c, 2) / R
                         for c in coords))
QHJ_resid = sp.simplify(sp.diff(S, tt)
                        + sum(sp.diff(S, c)**2 / (2 * masses[c]) for c in coords)
                        + Vg + Qexpr)
Bexpr = sp.simplify(sp.diff(R, tt)
                    + sum(sp.diff(R**2 * sp.diff(S, c) / masses[c], c)
                          for c in coords) / (2 * R))
# Exakte Identitaet (nicht vorausgesetzt, sondern als Behauptung formuliert):
#   R * resid  =  -R * QHJ_resid  +  i hbar * Bexpr,
# d.h. (i hbar d_t Psi - H Psi)/Psi = -QHJ + i hbar B/R. Da beide Seiten reell
# bzw. imaginaerteilig sind, folgt fuer R > 0:
#   SE  <=>  (QHJ_resid = 0 und Bexpr = 0).
diff_S1 = sp.simplify(R * resid + R * QHJ_resid - sp.I * hbar * Bexpr)
record("S1: Schroedinger-Gleichung  <=>  Kontinuitaet + Quanten-HJ,   Q = %s" % Qexpr,
       diff_S1 == 0, "Residual R*(SE/Psi) + R*QHJ - i hbar B = %s" % diff_S1)

# ====================================================================== P1 ====
cur_p1 = sp.simplify(hbar * (sp.conjugate(psi) * sp.diff(psi, q1)
                             - psi * sp.diff(sp.conjugate(psi), q1)) / (2 * sp.I)
                     / rho_of(psi) / m1)
diff_P1 = sp.simplify(cur_p1 - sp.diff(S, q1) / m1)
record("P1: Stromform j_i/rho  ==  Phasenform d_i S/m_i", diff_P1 == 0,
       "Differenz = %s" % diff_P1)

# ====================================================================== S2 ====
rho_g = R**2
cont_rho = sp.simplify(sp.diff(rho_g, tt)
                       + sum(sp.diff(rho_g * sp.diff(S, c) / masses[c], c)
                             for c in coords))
diff_S2 = sp.simplify(cont_rho - 2 * R * Bexpr)
record("S2: d_t rho + sum d_i(rho d_iS/m_i) = 0  ist dieselbe Aussage wie die R-Form",
       diff_S2 == 0, "rho-Form = %s ;  2R*B = %s" % (cont_rho, sp.simplify(2 * R * Bexpr)))

# ====================================================================== S3 ====
res_S3 = []
for c in coords:
    mi = masses[c]
    vc = sp.diff(S, c) / mi
    mat = sp.simplify(sp.diff(vc, tt)
                      + sum(sp.diff(S, k) / masses[k] * sp.diff(vc, k)
                            for k in coords))
    lhs = sp.simplify(mi * mat + sp.diff(Vg, c) + sp.diff(Qexpr, c))
    res_S3.append(sp.simplify(lhs - sp.diff(QHJ_resid, c)))
record("S3: Newton-Form m(d_t+v.grad)v_i = -d_i(V+Q)  aus d_i(Quanten-HJ)",
       all(r == 0 for r in res_S3), "Residuale (n=2) = %s" % res_S3)

# ============================================== explizite Oszillator-Zustaende =
x, t = sp.symbols("x t", real=True)
m, om = sp.Symbol("m", positive=True), sp.Symbol("omega", positive=True)
nrm0 = (m * om / (sp.pi * hbar))**sp.Rational(1, 4)
alpha = sp.sqrt(m * om / hbar)

psi0 = nrm0 * sp.exp(-m * om * x**2 / (2 * hbar)) * sp.exp(-sp.I * om * t / 2)
psi1 = nrm0 * sp.sqrt(2) * alpha * x \
    * sp.exp(-m * om * x**2 / (2 * hbar)) * sp.exp(-sp.I * 3 * om * t / 2)
Vho = sp.Rational(1, 2) * m * om**2 * x**2


def Q_of(psf, coord, mass):
    Rv = sp.sqrt(rho_of(psf))
    return sp.simplify(-(hbar**2 / (2 * mass)) * sp.diff(Rv, coord, 2) / Rv)


# ---------------------------------------------------------- S4: Grundzustand --
Q0 = Q_of(psi0, x, m)
E0 = hbar * om / 2
v0 = v_of(psi0, x, m)
# Erwartung (abgeleitet aus d^2 R/dx^2 / R = -m omega/hbar + m^2 omega^2 x^2/hbar^2):
#   Q = E0 - V   =>   V + Q = E0 konstant  =>  Kraft -d(V+Q)/dx = 0.
ok_S4 = (sp.simplify(Q0 - (E0 - Vho)) == 0
         and sp.simplify(sp.diff(Vho + Q0, x)) == 0
         and sp.simplify(v0) == 0)
record("S4: HO-Grundzustand  Q = E0 - V,  V+Q = E0 = %s konstant,  v = 0 (stationaere Bahn)" % E0,
       ok_S4, "Q = %s ; erwartet E0-V = %s ; d(V+Q)/dx = %s ; v = %s"
       % (Q0, sp.simplify(E0 - Vho), sp.simplify(sp.diff(Vho + Q0, x)), v0))

# -------------------------------------- S5: Superposition psi_0 + psi_1 -------
psisum = sp.simplify((psi0 + psi1) / sp.sqrt(2))
rhos = rho_of(psisum)
Qs = Q_of(psisum, x, m)
vs = v_of(psisum, x, m)
x_node = -1 / (sp.sqrt(2) * alpha)                    # Nullstelle von psi_0+psi_1 bei t=0
node_ok = sp.simplify(rhos.subs(t, 0).subs(x, x_node)) == 0 \
    and sp.simplify(rhos.subs(t, 0).subs(x, 0)) != 0
vs_late = sp.simplify(vs.subs(t, sp.pi / (2 * om)))
ok_S5 = node_ok and vs_late != 0 and sp.simplify(sp.diff(vs, t)) != 0 and is_real(vs)
record("S5: Superposition psi_0+psi_1: Dichteknoten bei x = %s, Fuehrung nichttrivial und zeitabhaengig"
       % x_node, ok_S5,
       "rho(Knoten,t=0) = %s ; v(t=pi/2omega) = %s ; dv/dt != 0: %s"
       % (sp.simplify(rhos.subs(t, 0).subs(x, x_node)), vs_late,
          sp.simplify(sp.diff(vs, t)) != 0))

# ------------------------------------------------------ S6: reeller psi -------
psi_real = sp.simplify(psi0.subs(t, 0) + psi1.subs(t, 0))   # reell bei t=0
v_real = v_of(psi_real, x, m)
record("S6: reelle Wellenfunktion (relative Phase 0 oder pi)  ->  v identisch 0",
       sp.simplify(v_real) == 0, "v = %s" % v_real)

# ------------------------------------------------- S7: Produkt vs. verschränkt -
# Zeitabhaengige Phaasen dazu, sonst ware der Zustand reell und v = 0 (S6).
psi_a0 = sp.exp(-m * om * q1**2 / (2 * hbar)) * sp.exp(-sp.I * om * tt / 2)
psi_a1 = q1 * sp.exp(-m * om * q1**2 / (2 * hbar)) * sp.exp(-sp.I * 3 * om * tt / 2)
psi_b0 = sp.exp(-m * om * q2**2 / (2 * hbar)) * sp.exp(-sp.I * om * tt / 2)
psi_b1 = q2 * sp.exp(-m * om * q2**2 / (2 * hbar)) * sp.exp(-sp.I * 3 * om * tt / 2)
prod_state = sp.simplify((psi_a0 + psi_a1) * psi_b0)
v1_prod = v_of(prod_state, q1, m)
dv1_dq2_prod = sp.simplify(sp.diff(v1_prod, q2))
ent_state = sp.simplify((psi_a0 * psi_b0 + psi_a1 * psi_b1) / sp.sqrt(2))
v1_ent = v_of(ent_state, q1, m)
dv1_dq2_ent = sp.simplify(sp.diff(v1_ent, q2))
ok_S7 = (dv1_dq2_prod == 0 and dv1_dq2_ent != 0
         and is_real(v1_ent) and is_real(dv1_dq2_ent))
record("S7: Produktzustand d v_1/d q_2 = 0;  verschaenkter (|00>+|11>)/sqrt2: d v_1/d q_2 != 0 "
       "(beide v reell)",
       ok_S7, "Produkt: %s ;  verschaenkt: %s ;  v_reell: %s"
       % (dv1_dq2_prod, dv1_dq2_ent, is_real(v1_ent)))

# ----------------------------------------------------------------- Ergebnis ---
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "setup": {
        "hamiltonian": "H = -sum_i hbar^2/(2 m_i) d_i^2 + V(q), n=2, allgemeine R,S,V",
        "decomposition": "psi = R exp(i S/hbar)",
        "quantum_potential": sp.sstr(Qexpr),
        "guidance_equivalence": "v_i = j_i/rho = d_i S/m_i",
    },
    "equivalence": {
        "continuity_rho_form": "d_t rho + sum_i d_i(rho d_iS/m_i) = 0",
        "quantum_hj": "d_t S + sum_i (d_iS)^2/(2m_i) + V + Q = 0",
        "newton_form": "m(d_t + v.grad)v_i = -d_i(V+Q)",
    },
    "explicit_states": {
        "ho_ground_Q": sp.sstr(Q0),
        "ho_superposition_v": sp.sstr(vs),
        "density_node": sp.sstr(x_node),
        "entangled_dv1_dq2": sp.sstr(dv1_dq2_ent),
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
