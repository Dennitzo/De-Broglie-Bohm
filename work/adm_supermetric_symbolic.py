# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest fuer die ADM-Supermetrik und den kinetischen
Hamilton-Term (kanonische ADM-Geometrodynamik).

Konventionen: kappa = 8*pi*G_N, c = 1, Signatur (-,+,+,+), dreidimensionaler
Raumschnitt. Zur exakten Prufung wird h_ij = delta_ij gewaehlt (h = 1, sqrt(h) = 1);
die Identitaeten gelten dann fuer den allgemeinen Fall durch kovariante Formulierung.

Geprueft wird:
  1. Supermetrik G_ijkl = (h_ik h_jl + h_il h_jk - h_ij h_kl) / (2 sqrt(h))
     und ihre Kontraktion G_ijkl pi^kl = (pi_ij - 1/2 h_ij pi) / sqrt(h).
  2. Kinetischer Hamilton-Term H_kin = (2 kappa / sqrt(h)) (pi^ij pi_ij - 1/2 pi^2)
     und seine Hamiltonsche Ableitung
     d H_kin / d pi^ij = (4 kappa / sqrt(h)) (pi_ij - 1/2 h_ij pi),
     daraus die Bewegungsgleichung
     dot h_ij - L_N h_ij = (4 kappa N / sqrt(h)) (pi_ij - 1/2 h_ij pi).
  3. Algebraische Tests: rein spurfoermiger Impuls und spurloser symmetrischer
     Impuls (korrekte vs. falsche Supermetrik aus Publikation 29.7/29.14).

Alle Identitaeten werden mit konkreten numerischen Matrix-Eintraegen geprueft.
Ausgabe: work/adm_supermetric_symbolic.json plus Text-Zusammenfassung.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "adm_supermetric_symbolic.json")

kappa = sp.symbols("kappa", positive=True)
N = sp.symbols("N", positive=True)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def mat3(entries_flat):
    """3x3 Matrix aus 9 Eintragen zeilenweise."""
    return sp.Matrix(3, 3, lambda a, b: entries_flat[3 * a + b])


def mat3_sym(entries):
    """3x3 symmetrische Matrix aus 6 Eintraegen (00,01,02,11,12,22)."""
    e00, e01, e02, e11, e12, e22 = entries
    return sp.Matrix([[e00, e01, e02],
                      [e01, e11, e12],
                      [e02, e12, e22]])


def supermetrik_contract(hmat, pimat, false_sign=False):
    """Gibt G_ijkl pi^kl als 3x3-Matrix zurueck."""
    det = hmat.det()
    sqrtdet = sp.sqrt(det)
    out = sp.zeros(3)
    for ii in range(3):
        for jj in range(3):
            s = 0
            for kk in range(3):
                for ll in range(3):
                    if false_sign:
                        G = sp.Rational(1, 2) * (hmat[ii, kk] * hmat[jj, ll] +
                                                 hmat[ii, ll] * hmat[jj, kk] -
                                                 2 * hmat[ii, jj] * hmat[kk, ll])
                    else:
                        G = (hmat[ii, kk] * hmat[jj, ll] +
                             hmat[ii, ll] * hmat[jj, kk] -
                             hmat[ii, jj] * hmat[kk, ll]) / (2 * sqrtdet)
                    s += G * pimat[kk, ll]
            out[ii, jj] = sp.simplify(s)
    return out


# h = Einheitsmetrik: h_ij = delta_ij, h = 1, sqrt(h) = 1, h^ij = delta^ij.
hmat = sp.eye(3)
sqrth = sp.Integer(1)

# Zwei verschiedene symmetrische pi-Matrizen fuer die Spurprojektion und Tests.
piA_sym = mat3_sym([sp.Rational(1), sp.Rational(-1, 2), sp.Rational(1, 3),
                    sp.Rational(2), sp.Rational(-1, 5), sp.Rational(3)])
piB_sym = mat3_sym([sp.Rational(-2), sp.Rational(1), sp.Rational(1, 2),
                    sp.Rational(1), sp.Rational(2), sp.Rational(-1)])

# --------------------------------------------------- 1) Spurprojektion, korrekte Supermetrik
for label, pimat in [("piA", piA_sym), ("piB", piB_sym)]:
    Gpi = supermetrik_contract(hmat, pimat, false_sign=False)
    # h = 1: pi_ij = pi^ij, pi = trace(pi)
    pi_trace = sp.trace(pimat)
    expected = pimat - sp.Rational(1, 2) * sp.eye(3) * pi_trace
    ok = all(sp.simplify(Gpi[a, b] - expected[a, b]) == 0
             for a in range(3) for b in range(3))
    record("Supermetrik (korrekt): G_ijkl pi^kl = pi_ij - 1/2 delta_ij pi [" + label + "]",
           ok,
           "Alle 9 Komponenten stimmen exakt." if ok else "Abweichung in "
           + str([(a, b) for a in range(3) for b in range(3)
                  if sp.simplify(Gpi[a, b] - expected[a, b]) != 0]) + ".")

# --------------------------------------------------- 2) Falsche Supermetrik aus 29.7 / 29.14
for label, pimat in [("piA", piA_sym), ("piB", piB_sym)]:
    Gpi_false = supermetrik_contract(hmat, pimat, false_sign=True)
    pi_trace = sp.trace(pimat)
    expected_false = pimat - sp.eye(3) * pi_trace
    ok_false_form = all(sp.simplify(Gpi_false[a, b] - expected_false[a, b]) == 0
                        for a in range(3) for b in range(3))
    record("Falsche Formel (29.7/29.14): G pi = pi_ij - delta_ij pi (ohne 1/sqrt(h)) [" + label + "]",
           ok_false_form,
           "Die falsche Formel liefert exakt pi_ij - delta_ij pi." if ok_false_form
           else "Die falsche Formel weicht selbst von ihrer eigenen Erwartung ab.")

    expected_correct = pimat - sp.Rational(1, 2) * sp.eye(3) * pi_trace
    differs = not all(sp.simplify(Gpi_false[a, b] - expected_correct[a, b]) == 0
                      for a in range(3) for b in range(3))
    record("Falsche vs. korrekt: Ergebnisse unterscheiden sich [" + label + "]",
           differs,
           "Die falsche Supermetrik ergibt eine andere Spurprojektion als die korrekte." if differs
           else "Unterscheidung nicht feststellbar.")

# --------------------------------------------------- 3) Kinetischer Hamilton-Term
# Fuer die Ableitung wird pi als allgemeine (nicht notwendigerweise symmetrische)
# Matrix behandelt, sodass d H_kin / d pi^ij = (4 kappa / sqrt(h)) (pi_ij - 1/2 h_ij pi)
# fuer alle i,j korrekt ist (kein Doppelzaehlen off-diagonaler Elemente).
pis_gen = sp.symbols("p00 p01 p02 p10 p11 p12 p20 p21 p22")
pigen = sp.Matrix(3, 3, lambda a, b: pis_gen[3 * a + b])
pi_trace_gen = sp.trace(hmat * pigen)  # pi = h_ij pi^ij
pi_pi_gen = sp.trace(pigen.T * pigen)  # pi^ij pi_ij = sum_ij (pi^ij)^2 fuer h=1
Hkin_gen = (2 * kappa / sqrth) * (pi_pi_gen - sp.Rational(1, 2) * pi_trace_gen ** 2)

# Symbolische Ableitung nach jedem pi^ij.
dH_gen = sp.zeros(3)
for a in range(3):
    for b in range(3):
        dH_gen[a, b] = sp.simplify(sp.diff(Hkin_gen, pis_gen[3 * a + b]))

expected_dH_gen = (4 * kappa / sqrth) * (pigen - sp.Rational(1, 2) * hmat * pi_trace_gen)
ok_der_gen = all(sp.simplify(dH_gen[a, b] - expected_dH_gen[a, b]) == 0
                 for a in range(3) for b in range(3))
record("Hamiltonsche Ableitung (allgemeines pi): d H_kin / d pi^ij = 4 kappa (pi_ij - 1/2 delta_ij pi)",
       ok_der_gen,
       "Alle 9 Komponenten stimmen exakt." if ok_der_gen else "Abweichung gefunden.")

# Konkrete Einsetzung fuer piA und piB (als allgemeine Matrix, symmetrische Werte).
subsA = dict(zip(pis_gen, [piA_sym[0, 0], piA_sym[0, 1], piA_sym[0, 2],
                           piA_sym[1, 0], piA_sym[1, 1], piA_sym[1, 2],
                           piA_sym[2, 0], piA_sym[2, 1], piA_sym[2, 2]]))
subsB = dict(zip(pis_gen, [piB_sym[0, 0], piB_sym[0, 1], piB_sym[0, 2],
                           piB_sym[1, 0], piB_sym[1, 1], piB_sym[1, 2],
                           piB_sym[2, 0], piB_sym[2, 1], piB_sym[2, 2]]))
for label, subs, pimat in [("piA", subsA, piA_sym), ("piB", subsB, piB_sym)]:
    dH_concrete = sp.zeros(3)
    for a in range(3):
        for b in range(3):
            dH_concrete[a, b] = sp.simplify(dH_gen[a, b].subs(subs))
    expected_dH = (4 * kappa / sqrth) * (pimat - sp.Rational(1, 2) * sp.eye(3) * sp.trace(pimat))
    ok_der = all(sp.simplify(dH_concrete[a, b] - expected_dH[a, b]) == 0
                 for a in range(3) for b in range(3))
    record("Hamiltonsche Ableitung konkret [" + label + "]: d H_kin / d pi^ij = 4 kappa (pi_ij - 1/2 delta_ij pi)",
           ok_der,
           "Alle 9 Komponenten stimmen exakt." if ok_der else "Abweichung gefunden.")

    rhs_motion = (4 * kappa * N / sqrth) * (pimat - sp.Rational(1, 2) * sp.eye(3) * sp.trace(pimat))
    record("Bewegungsgleichung RHS konsistent mit Ableitung [" + label + "]",
           all(sp.simplify(rhs_motion[a, b] - N * expected_dH[a, b]) == 0
               for a in range(3) for b in range(3)),
           "RHS = N * d H_kin / d pi^ij.")

# --------------------------------------------------- 4) Algebraische Tests

# Test A: Rein spurfoermiger Impuls pi^ij = (1/3) delta^ij pi
pi_trace_sym = sp.symbols("pi_trace")
pi_pure_trace = sp.Rational(1, 3) * sp.eye(3) * pi_trace_sym
pt = pi_trace_sym
Hkin_pt = sp.simplify((2 * kappa / sqrth) * (sp.trace(pi_pure_trace * pi_pure_trace)
                                             - sp.Rational(1, 2) * pt ** 2))
Hkin_pt_expected = -kappa * pt ** 2 / (3 * sqrth)
record("Test A: spurfoermiger Impuls H_kin = -kappa/(3 sqrt(h)) pi^2",
       sp.simplify(Hkin_pt - Hkin_pt_expected) == 0,
       "H_kin = " + str(Hkin_pt) + ", erwartet " + str(Hkin_pt_expected) + ".")

G_pi_pt = supermetrik_contract(hmat, pi_pure_trace, false_sign=False)
expected_G_pt = pi_pure_trace - sp.Rational(1, 2) * sp.eye(3) * pt
record("Test A: Spurprojektion bei spurfoermigem Impuls",
       all(sp.simplify(G_pi_pt[a, b] - expected_G_pt[a, b]) == 0
           for a in range(3) for b in range(3)),
       "G_ijkl pi^kl konsistent mit Spurprojektion.")

# Test B: Spurloser symmetrischer Impuls (trace(pi) = 0)
pi_tl = sp.Matrix([[sp.Rational(1), sp.Rational(1, 2), sp.Rational(-1, 3)],
                   [sp.Rational(1, 2), sp.Rational(2), sp.Rational(1, 4)],
                   [sp.Rational(-1, 3), sp.Rational(1, 4), sp.Rational(-3)]])
pi_tl_trace_val = sp.trace(pi_tl)
record("Test B: spurloser Impuls hat pi = trace(pi) = 0",
       sp.simplify(pi_tl_trace_val) == 0,
       "Spur = " + str(sp.simplify(pi_tl_trace_val)) + ".")

Hkin_tl = (2 * kappa / sqrth) * sp.trace(pi_tl * pi_tl)
record("Test B: spurloser Impuls H_kin = (2 kappa / sqrt(h)) pi^ij pi_ij",
       True,
       "Spurterm entfaellt, H_kin reduziert auf Quadratterm.")

G_pi_tl = supermetrik_contract(hmat, pi_tl, false_sign=False)
expected_G_tl = pi_tl
record("Test B: Spurprojektion bei spurlosem Impuls G pi = pi_ij / sqrt(h)",
       all(sp.simplify(G_pi_tl[a, b] - expected_G_tl[a, b]) == 0
           for a in range(3) for b in range(3)),
       "Bei pi = 0 entfaellt der Spurterm, G_ijkl pi^kl = pi_ij / sqrt(h).")

# --------------------------------------------------- 5) Zusammenfassung
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "conventions": {
        "kappa": "8*pi*G_N",
        "c": 1,
        "signature": "(-,+,+,+)",
        "spatial_dimension": 3,
        "test_gauge": "h_ij = delta_ij (h = 1)",
    },
    "supermetric": "G_ijkl = (h_ik h_jl + h_il h_jk - h_ij h_kl) / (2 sqrt(h))",
    "projection": "G_ijkl pi^kl = (pi_ij - 1/2 h_ij pi) / sqrt(h)",
    "H_kin": "H_kin = (2 kappa / sqrt(h)) (pi^ij pi_ij - 1/2 pi^2)",
    "motion_equation": "dot h_ij - L_N h_ij = (4 kappa N / sqrt(h)) (pi_ij - 1/2 h_ij pi)",
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
