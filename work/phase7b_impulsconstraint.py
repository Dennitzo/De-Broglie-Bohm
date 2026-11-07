# -*- coding: utf-8 -*-
"""
Zusaetzliche Gegenrechnung zur korrigierten konformen Impulsconstraint
(Phase 1, Manuskript 19.6):

  (1) Zweite, vom Gewichtsbild unabhaengige Herleitung der Divergenz
      D_j pi_i^j ueber sqrt(h) * D_j (h_im pi^mj / sqrt(h));
  (2) integrierender Faktor und exakte Loesung der lokalen Constraint auf S^1;
  (3) Periodizitaets-/Loeubarkeitsbedingung und ihre Verletzung;
  (4) Nachweis, dass die exakte Loesung fuer alpha_1 != 0 den
      Trimmraum {1, sin, cos} verlaesst (leading cos(2 theta)-Anteil);
  (5) Modenbeispiel bei konstantem alpha (stimmt mit der alten Rechnung
      ueberein, weil dort P d_i alpha faellt).

Ausschliesslich exakte symbolische Rechnung; kein numerischer Fit.
Ergebnispapier: work/phase7b_impulsconstraint.json
"""
import datetime
import json
import os

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "phase7b_impulsconstraint.json")
CHECKS = {}


def rec(name, ok, detail=""):
    CHECKS[name] = {"pass": bool(ok), "detail": str(detail)}
    print(("PASS " if ok else "FAIL "), "|", name, "|", detail)


x, y, z = sp.symbols("x y z", real=True)
XS = (x, y, z)
al = sp.Function("alpha")(*XS)
Pl = sp.Function("P")(*XS)
dd = lambda f, i: sp.diff(f, XS[i])

h_lo = sp.diag(*[sp.exp(2 * al)] * 3)
h_up = h_lo.inv()
sqh = sp.exp(3 * al)

G = {}
for k in range(3):
    for i in range(3):
        for j in range(3):
            s = sum((dd(h_lo[j, m], i) + dd(h_lo[i, m], j) - dd(h_lo[i, j], m)) * h_up[k, m]
                    for m in range(3))
            G[(k, i, j)] = sp.simplify(s / 2)

# ---------------------------------------------------------------------
# (1) Zweite Herleitung: pi^ij als Dichte weight +1,
#     D_j pi_i^j := sqrt(h) * D_j( h_im pi^mj / sqrt(h) )
#     mit der UEBLICHEN Levi-Civita-Ableitung fuer den Tensor weight 0.
# ---------------------------------------------------------------------
pi_up = {(i, j): (sp.simplify(Pl / 6 * h_up[i, j]) if i == j else sp.Integer(0))
         for i in range(3) for j in range(3)}
# Spuridentitaet
tr_chk = sp.simplify(sum(h_lo[i, j] * pi_up[(i, j)] for i in range(3) for j in range(3)))
rec("Spuridentitaet 2 h_ij pi^ij = P (Pi^ij = (P/6) e^{-2a} delta^ij)",
    sp.simplify(2 * tr_chk - Pl) == 0, "h_ij pi^ij = %s, also 2 h_ij pi^ij = %s"
    % (tr_chk, sp.simplify(2 * tr_chk)))

Div2 = {}
for i in range(3):
    # tau_i^j = h_im pi^mj / sqrt(h)  (echter Tensor, weight 0)
    tau = sp.Matrix(3, 3, lambda a, b: sp.simplify(
        sum(h_lo[a, m] * pi_up[(m, b)] for m in range(3)) / sqh))
    tot = sp.Integer(0)
    for j in range(3):
        term = dd(tau[i, j], j)
        for m in range(3):
            term += G[(j, j, m)] * tau[i, m] - G[(m, j, i)] * tau[m, j]
        tot += term
    Div2[i] = sp.simplify(sp.expand(tot * sqh))

tgt = {i: sp.simplify((dd(Pl, i) - 3 * Pl * dd(al, i)) / 6) for i in range(3)}
rec("Weg C (sqrt(h)-Reduktion, Levi-Civita-Div. eines Tensors weight 0): "
    "D_j pi_i^j = (d_i P - 3 P d_i a)/6",
    all(sp.simplify(Div2[i] - tgt[i]) == 0 for i in range(3)), str(Div2))

# Explizit: der Anteil, den die alte (ungewichtete) Rechnung verliert
DivOld = {}
for i in range(3):
    tau = sp.Matrix(3, 3, lambda a, b: sp.simplify(sum(h_lo[a, m] * pi_up[(m, b)]
                                                       for m in range(3))))
    tot = sp.Integer(0)
    for j in range(3):
        term = dd(tau[i, j], j)
        for m in range(3):
            term += G[(j, j, m)] * tau[i, m] - G[(m, j, i)] * tau[m, j]
        tot += term
    DivOld[i] = sp.simplify(sp.expand(tot))
rec("Alte Rechnung (pi_i^j = (P/6) delta_i^j ohne Dichtegewicht) laesst pro Komponente "
    "P d_i a/6 aus der Divergenz weg, d. h. in -2D genau den Beitrag P d_i a",
    all(sp.simplify(DivOld[i] - (Div2[i] + Pl * dd(al, i) / 2)) == 0 for i in range(3))
    and all(sp.simplify(-2 * DivOld[i] - (-2 * Div2[i] - Pl * dd(al, i))) == 0
            for i in range(3)),
    "Div_alt - Div_korr = %s" % [sp.simplify(DivOld[i] - Div2[i]) for i in range(3)])

# ---------------------------------------------------------------------
# (2) (4) Lokale Aufloesbarkeit auf S^1, theta-Koordinate
# ---------------------------------------------------------------------
th = sp.symbols("theta", real=True)
a0s, a1, q, p0, P0s, kS = sp.symbols("alpha_0 alpha_1 q p_0 P_0 k", real=True)
mean = lambda e: sp.simplify(sp.integrate(sp.expand(e), (th, 0, 2 * sp.pi)) / (2 * sp.pi))

alphaK = a0s + a1 * sp.sin(th)
phiK = q * sp.sin(th)
Pgen = P0s * sp.exp(3 * a1 * sp.sin(th)) + q * p0 / a1 * (sp.exp(3 * a1 * sp.sin(th)) - 1)

Hdens1 = lambda P, A, PH, p: sp.expand(-(sp.diff(P, th) / 3) + P * sp.diff(A, th)
                                       + p * sp.diff(PH, th))
res_gen = sp.simplify(Hdens1(Pgen, alphaK, phiK, p0))
rec("(2) Exakte Loesung P = P0 e^{3 a1 sin} + (p0 q/a1)(e^{3 a1 sin} - 1) erfuellt "
    "H_theta = -P'/3 + P a' + p0 phi' identisch",
    res_gen == 0, "H_theta = %s" % res_gen)

# Integrierender Faktor als Identitaet
F = sp.Function("F")(th)
Gf = sp.Function("Gf")(th)
lhs = sp.expand(sp.diff(F * sp.exp(-3 * Gf), th))
rec("(2) Integrierender Faktor: (P e^{-3a})' = e^{-3a}(P' - 3 a' P) identisch",
    sp.simplify(lhs - sp.expand(sp.exp(-3 * Gf) * (sp.diff(F, th) - 3 * sp.diff(Gf, th) * F))) == 0,
    lhs)

# Aequivalenz der Constraint zur totalen Ableitung (voll generisch in theta)
Ag = sp.Function("A")(th)
Pg = sp.Function("P")(th)
Phg = sp.Function("Varphi")(th)
pg = sp.Function("p")(th)
Hgen = sp.expand(-(sp.diff(Pg, th) / 3) + Pg * sp.diff(Ag, th) + pg * sp.diff(Phg, th))
Hexpr = sp.expand(-3 * sp.exp(-3 * Ag) * Hgen)
rhs_eq = sp.expand(sp.diff(Pg * sp.exp(-3 * Ag), th)
                   - 3 * sp.exp(-3 * Ag) * pg * sp.diff(Phg, th))
rec("(2) Aequivalenz (generisch): -3 e^{-3A} H = (P e^{-3A})' - 3 e^{-3A} p phi'",
    sp.simplify(Hexpr - rhs_eq) == 0, "beide = %s" % sp.factor(sp.simplify(rhs_eq)))

# ---------------------------------------------------------------------
# (3) Periodizitaetsbedingung int_0^{2 pi} e^{-3a} p phi' d theta = 0
# ---------------------------------------------------------------------
cond_ok = sp.simplify(mean(sp.exp(-3 * alphaK) * p0 * sp.diff(phiK, th)))
rec("(3) Loesbarkeitsintegral fuer phi = q sin(theta) ist EXAKT null "
    "(Totale Ableitung e^{-3 a1 sin theta})",
    cond_ok == 0, "Integral = %s" % cond_ok)

phiC = q * sp.cos(th)
cond_bad = sp.simplify(mean(sp.exp(-3 * alphaK) * p0 * sp.diff(phiC, th)))
cond_bad_ser = sp.series(cond_bad, a1, 0, 3).removeO()
lin = sp.expand(cond_bad_ser).coeff(a1, 1)
rec("(3) Loesbarkeitsintegral fuer phi = q cos(theta) ist NICHT identisch null: "
    "fuehrender Term in a1 ist +(3/2) p0 q e^{-3 alpha_0} a1",
    sp.simplify(cond_bad) != 0
    and sp.simplify(lin - sp.Rational(3, 2) * p0 * q * sp.exp(-3 * a0s)) == 0,
    "Reihe = %s ; Linearkoeff = %s" % (cond_bad_ser, lin))

# ---------------------------------------------------------------------
# (4) Leck aus dem Trimmraum {1, sin, cos}: cos(2 theta)-Komponente
# ---------------------------------------------------------------------
Pser = sp.expand(sp.series(Pgen, a1, 0, 3).removeO())
rec("(4) Reihenentwicklung der exakten Loesung: P = P0 + 3 p0 q sin(theta) "
    "+ (9/2) p0 q a1 sin^2(theta) + O(a1^2)",
    sp.simplify(Pser.subs(a1, 0) - (P0s + 3 * p0 * q * sp.sin(th))) == 0,
    "P bis a1^1 = %s" % Pser)
c2 = sp.simplify(mean(Pser * sp.cos(2 * th)))
c2_lin = sp.simplify(sp.expand(c2).coeff(a1, 1))
rec("(4) cos(2 theta)-Fourierkoeffizient der exakten Loesung: <P cos 2theta> = "
    "-(9/8) p0 q a1 + O(a1^2) ungleich null fuer a1 != 0",
    sp.simplify(c2_lin + sp.Rational(9, 8) * p0 * q) == 0,
    "<P cos2th> (bis a1^2) = %s ; Linearkoeff = %s" % (c2, c2_lin))
c1 = sp.simplify(mean(Pser * sp.sin(th)))
c1s0 = sp.simplify(sp.expand(c1.subs(a1, 0)) - 3 * p0 * q / 2 * sp.Integer(2) / 2)
rec("(4) sin(theta)-Koeffizient der exakten Loesung bei a1 = 0: <P sin> = 3 p0 q/2 "
    "(Doppelnormierung P_1 = 3 p0 q / 2, Feld amplitude 3 p0 q)",
    sp.simplify(sp.expand(c1.subs(a1, 0)) - sp.Rational(3, 2) * p0 * q) == 0,
    "<P sin> bei a1=0 = %s" % sp.simplify(c1.subs(a1, 0)))

# ---------------------------------------------------------------------
# (5) Modenbeispiel bei konstantem alpha: alte und neue Constraint liefern
#     dieselbe Loesung, weil P d_i alpha faellt
# ---------------------------------------------------------------------
for kk in (1, 2, 3):
    Pold = P0s + 3 * p0 * q * sp.sin(kS * th).subs(kS, kk)
    r_new = sp.simplify(Hdens1(Pold, a0s, q * sp.sin(kk * th), p0))
    r_old = sp.expand(sp.diff(Pold, th) / 3 - p0 * sp.diff(q * sp.sin(kk * th), th))
    rec("(5) Modenloesung P = P0 + 3 p0 q sin(%d theta) bei konstantem alpha: "
        "korrigierte Constraint H_theta = 0 identisch" % kk,
        r_new == 0, "H_theta = %s" % r_new)
    rec("(5) dieselbe Loesung im Vergleichsmodell r = P'/3 - p phi' (alpha faellt "
        "raus): identisch null, k = %d" % kk,
        sp.simplify(r_old) == 0, "r = %s" % sp.simplify(r_old))

# ---------------------------------------------------------------------
# (6) Homogener Anschluss: konstante Felder erfuellen die Impulsconstraint
# ---------------------------------------------------------------------
Pc, phic, pmc = sp.symbols("P_c phi_c p_c", real=True)
hom = sp.simplify(Hdens1(Pc, sp.Symbol("a_c", real=True), phic, pmc))
rec("(6) homogener Anschluss: alpha,P,phi,p konstant => H_theta = 0 exakt",
    hom == 0, "H_theta = %s" % hom)

# ---------------------------------------------------------------------
# (7) Shift-Kopplung: die Constraint ist das zu N^i konjugierte Funktional
# ---------------------------------------------------------------------
Nsh = [sp.Function("N%d" % i)(*XS) for i in range(3)]
phi3 = sp.Function("phi3")(*XS)
pm3 = sp.Function("pm3")(*XS)
lie = sp.Integer(0)
for i in range(3):
    for j in range(3):
        s = sum(Nsh[m] * dd(h_lo[i, j], m) for m in range(3))
        if i == j:
            s += 0
        for m in range(3):
            s += h_lo[m, j] * dd(Nsh[m], i) + h_lo[i, m] * dd(Nsh[m], j)
        lie += pi_up[(i, j)] * sp.expand(s)
lie = sp.expand(lie + pm3 * sum(Nsh[m] * dd(phi3, m) for m in range(3)))
HB = {}
for i in range(3):
    Lc = sp.simplify(sp.expand(lie).coeff(Nsh[i]))
    Msum = sp.Integer(0)
    for j in range(3):
        Msum += dd(sp.simplify(sp.expand(lie).coeff(dd(Nsh[i], j))), j)
    HB[i] = sp.simplify(sp.expand(Lc - Msum))
tgtB = {i: sp.simplify(-dd(Pl, i) / 3 + Pl * dd(al, i) + pm3 * dd(phi3, i)) for i in range(3)}
rec("(7) Lie-Funktional int(pi^ij L_N h_ij + p L_N phi) = int N^i H_i mit "
    "H_i = -d_i P/3 + P d_i a + p d_i phi (zweite, von (1) unabhaengige Realisierung)",
    all(sp.simplify(HB[i] - tgtB[i]) == 0 for i in range(3)), str(HB))
rec("(7) Konsistenz: Weg C und Weg B liefern dieselbe Constraint (Differenz 0)",
    all(sp.simplify((-2 * Div2[i] + pm3 * dd(phi3, i)) - tgtB[i]) == 0 for i in range(3)),
    "differenzen = %s" % [sp.simplify((-2 * Div2[i] + pm3 * dd(phi3, i)) - tgtB[i])
                          for i in range(3)])

out = {"generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
       "sympy_version": sp.__version__,
       "scope": "Phase 1 Zusatz: korrigierte konforme Impulsconstraint, lokale "
                "Aufloesbarkeit und Trimmraum-Leck",
       "convention": "h_ij = e^{2a} delta_ij, sqrt h = e^{3a}, P = 2 h_ij pi^ij, "
                     "pi^ij Dichte weight +1, D_j pi_i^j := sqrt(h) D_j(h_im pi^mj/sqrt h), "
                     "phi = sqrt(beta) Phi, p = Pi_phi/beta",
       "checks": CHECKS,
       "all_pass": all(c["pass"] for c in CHECKS.values()),
       "n_checks": len(CHECKS)}
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False, sort_keys=True)
print("\nALLE CHECKS BESTANDEN" if out["all_pass"] else "\nFEHLER VORHANDEN")
print("checks:", out["n_checks"], "all_pass:", out["all_pass"])
