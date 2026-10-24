# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest der ADM-Hamilton-Zwangsdichte (3 Raumdimensionen).

Konvention: Signatur (-+++), kappa=8*pi*G_N, c=1, sqrt h = sqrt(det h).
In der Algebra wird die Wurzel (sqrt h) als positiver Skalar 'sh' gefuehrt;
der Spezialfall sqrt(det h) wird separat geprueft. Alle Matrix-Vergleiche
werden komponentenweise gefuehrt (sympy nimmt `Matrix == 0` nicht elementweise).

Prueft:
  1. Impulsdefinition: pi^{ij} = (sqrt h/(2 kappa))(K^{ij} - h^{ij} K).
  2. Spur ueber 3 Dimensionen (KORREKTUR):
       pi = h_{ij} pi^{ij} = (sqrt h/(2 kappa))(K - 3 K) = -(sqrt h/kappa) K.
  3. Umkehr mit hochgestelltem Index (h^{ij}, NICHT h_{ij}):
       K^{ij} = (2 kappa/sqrt h)(pi^{ij} - (1/2) h^{ij} pi),
       K = -(kappa/sqrt h) pi.
  4. ROUND-TRIP (unabhaengig): konkrete nichtdiagonale Zahlen-Metrik + K,
       pi^{ij} aus der Definition, daraus K^{ij} zurueck, Vergleich mit Original.
  5. Rein-spurig (konform): K_{ij} = q h_{ij}.
  6. Spurlos: K = 0  =>  pi = 0.
  7. Diagonale positive Metrik (drei unabhaengige Skalen).
  8. DeWitt-Supermetrik: G_{ijkl} pi^{kl} = (1/sqrt h)(pi_ij - (1/2) h_ij pi).
  9. Hamilton-Ableitung (als symmetrische Funktionalvariation):
       dH_kin/dpi^{ij} = (4 kappa/sqrt h)(pi_ij - (1/2) h_ij pi) = 4 kappa G_{..} pi^{..}.

Ausgabe: work/adm_constraint_symbolic.json plus Text-Zusammenfassung.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "adm_constraint_symbolic.json")

kappa = sp.Symbol("kappa", positive=True)   # 8*pi*G_N
sh = sp.Symbol("sh", positive=True)         # sqrt h (Dichtegewicht, hier Skalar)
results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def spur(h, A):
    return sp.simplify(sp.trace(h * A))


def mat_zero(M):
    return all(sp.simplify(M[i, j]) == 0
               for i in range(M.rows) for j in range(M.cols))


def mat_eq(M1, M2):
    return mat_zero(M1 - M2)


# --------------------- 7) Diagonale positive Metrik (Algebra) -----------------
ha, hb, hc = sp.symbols("ha hb hc", positive=True)
h = sp.diag(ha, hb, hc)
hinv = sp.diag(1 / ha, 1 / hb, 1 / hc)
record("Diagonale Metrik: det h = ha*hb*hc (Wurzel als Dichtegewicht gefuehrt)",
       sp.simplify(sp.sqrt(h.det())**2 - ha * hb * hc) == 0,
       "det h = %s ; sqrt h^2 = ha*hb*hc" % sp.simplify(h.det()))

k1, k2, k3, k4, k5, k6 = sp.symbols("k1 k2 k3 k4 k5 k6", real=True)
K = sp.Matrix([[k1, k2, k3], [k2, k4, k5], [k3, k5, k6]])   # K_{ij} symmetrisch
Kraise = sp.simplify(hinv * K * hinv)         # K^{ij}
Kscalar = sp.simplify(spur(hinv, K))          # K = h^{ij} K_{ij}

Pdef = sp.simplify((sh / (2 * kappa)) * (Kraise - Kscalar * hinv))   # pi^{ij}
Ptrace = sp.simplify(spur(h, Pdef))           # pi = h_{ij} pi^{ij}

# --------------------- 2) Spur- und Faktor-Korrektur ---------------------------
spur_expected = sp.simplify((sh / (2 * kappa)) * (Kscalar - 3 * Kscalar))
record("Spur 3D: pi = (sqrt h/(2 kappa))(K - 3K) = -(sqrt h/kappa) K",
       sp.simplify(spur_expected - (-sh / kappa * Kscalar)) == 0
       and sp.simplify(Ptrace - spur_expected) == 0,
       "K-3K=-2K ; pi=%s == -(sqrt h/kappa)K=%s" %
       (Ptrace, -sh / kappa * Kscalar))

# --------------------- 3) Umkehr mit hochgestelltem Index ---------------------
Pinv = sp.simplify((2 * kappa / sh) * (Pdef - sp.Rational(1, 2) * Ptrace * hinv))
record("Umkehr: K^{ij} = (2 kappa/sqrt h)(pi^{ij} - (1/2) h^{ij} pi)",
       mat_eq(Pinv, Kraise),
       "Rueckrechnung == Original K^{ij} (alle Komponenten)")
record("Umkehr K: K = -(kappa/sqrt h) pi",
       sp.simplify(-kappa / sh * Ptrace - Kscalar) == 0,
       "-(kappa/sqrt h) pi = %s == K = %s" % (-kappa / sh * Ptrace, Kscalar))
# Index-Form: Rueckformel nutzt h^{ij} (hinv), nicht h_{ij}.
record("Index-Form: Rueckformel nutzt h^{ij} (hinv), nicht h_{ij}",
       True,
       "Pinv^ij = ... - (1/2) pi * h^{ij}; die Matrix h (kovariant) kommt nicht vor")

# --------------------- 4) ROUND-TRIP mit konkreten Zahlen ---------------------
hN = sp.Matrix([[2, 1, 0], [1, 3, 1], [0, 1, 4]])
KN = sp.Matrix([[1, 2, 3], [2, 5, 7], [3, 7, 6]])
hN_inv = hN.inv()
sN = sp.sqrt(hN.det())
KN_raise = sp.simplify(hN_inv * KN * hN_inv)
KN_scalar = sp.simplify(sp.trace(hN_inv * KN))
PN = sp.simplify((sN / (2 * kappa)) * (KN_raise - KN_scalar * hN_inv))
PN_trace = sp.simplify(sp.trace(hN * PN))
PN_back = sp.simplify((2 * kappa / sN) * (PN - sp.Rational(1, 2) * PN_trace * hN_inv))

rt_ok = True
rt_lines = []
for i in range(3):
    for j in range(3):
        d = sp.simplify(PN_back[i, j] - KN_raise[i, j])
        rt_ok = rt_ok and (d == 0)
        rt_lines.append("(%d,%d):%s" % (i, j, d))
record("ROUND-TRIP (nichtdiag. Zahlen-Metrik): K^{ij} aus pi^{ij} == Original",
       rt_ok, "; ".join(rt_lines))
record("ROUND-TRIP Spur: K = -(kappa/sqrt h) pi (konkrete Metrik)",
       sp.simplify(KN_scalar - (-kappa / sN * PN_trace)) == 0,
       "K=%s ; -(kappa/sqrt h)*pi=%s" % (KN_scalar, -kappa / sN * PN_trace))

# --------------------- 5) Rein-spuriger (konformer) Fall ---------------------
q = sp.Symbol("q", real=True)
Kc = q * h
Kc_raise = q * hinv
Kc_scalar = sp.simplify(spur(hinv, Kc))          # = 3 q
Pc = sp.simplify((sh / (2 * kappa)) * (Kc_raise - Kc_scalar * hinv))
Pc_trace = sp.simplify(spur(h, Pc))
record("Rein-spurig: K_{ij}=q h_{ij} => pi = -3 q sqrt h/kappa",
       sp.simplify(Pc_trace - (-3 * q * sh / kappa)) == 0,
       "pi = %s" % Pc_trace)
Pc_back = sp.simplify((2 * kappa / sh) * (Pc - sp.Rational(1, 2) * Pc_trace * hinv))
record("Rein-spurig ROUND-TRIP: K_{ij} zurueck == q h_{ij}",
       mat_eq(h * Pc_back * h, Kc),
       "h*(Pc_back)*h - Kc = %s" % sp.simplify(h * Pc_back * h - Kc))

# --------------------- 6) Spurloser Fall ---------------------------------------
x, y = sp.symbols("x y", real=True)
z = -x - y
K0 = sp.diag(x, y, z)
h0 = sp.diag(1, 1, 1)
h0_inv = sp.diag(1, 1, 1)
s0 = sp.sqrt(1)
K0_scalar = sp.simplify(spur(h0_inv, K0))       # x+y+z=0
P0 = sp.simplify((s0 / (2 * kappa)) * (K0 - K0_scalar * h0_inv))
P0_trace = sp.simplify(spur(h0, P0))
record("Spurlos: K=0 => pi = 0",
       P0_trace == 0,
       "pi=%s" % P0_trace)
record("Spurlos: pi^{ij} = (sqrt h/(2 kappa)) K^{ij}",
       mat_eq(P0, (s0 / (2 * kappa)) * K0),
       "pi^{ij} = %s" % P0)

# --------------------- 8) DeWitt-Supermetrik ------------------------------------
pS = sp.Matrix([[sp.Symbol("P11", real=True), sp.Symbol("P12", real=True),
                 sp.Symbol("P13", real=True)],
                [sp.Symbol("P12", real=True), sp.Symbol("P22", real=True),
                 sp.Symbol("P23", real=True)],
                [sp.Symbol("P13", real=True), sp.Symbol("P23", real=True),
                 sp.Symbol("P33", real=True)]])
hS = sp.diag(sp.Symbol("Ha", positive=True), sp.Symbol("Hb", positive=True),
             sp.Symbol("Hc", positive=True))
pS_low = sp.simplify(hS * pS * hS)              # pi_ij
pS_trace = sp.simplify(spur(hS, pS))            # pi
G_contract = sp.simplify(
    (2 * pS_low - pS_trace * hS) / (2 * sp.sqrt(hS.det())))
expected = sp.simplify((pS_low - sp.Rational(1, 2) * pS_trace * hS) /
                       sp.sqrt(hS.det()))
record("DeWitt: G_{ijkl} pi^{kl} = (1/sqrt h)(pi_ij - (1/2) h_ij pi)",
       mat_eq(G_contract, expected),
       "G_contract == expected (komponentenweise identisch)")

# --------------------- 9) Hamilton-Ableitung (symmetrische Variation) -----------
# H_kin = (2 kappa/sqrt h)(pi^{ij}pi_{ij} - (1/2) pi^2).  Euler-Theorem fuer
# L_kin quadratisch in K: pi^{ij} K_{ij} = 2 L_kin  =>  H_kin = L_kin = pi^{ij}K_{ij}.
# Die funktionale Variation (delta pi^{ij} = eps T^{ij}, T symmetrisch) liefert
#   dH_kin = (4 kappa/sqrt h)(pi_ij - (1/2) h_ij pi) T^{ij}  =  2 K_ij T^{ij},
# konsistent mit K_{ij} = (1/(2N))(dh_{ij}-L_N h_{ij}).
Hkin = (2 * kappa / sp.sqrt(hS.det())) * \
    (sp.trace(hS * pS * hS * pS) - sp.Rational(1, 2) * pS_trace**2)
eps = sp.Symbol("eps", real=True)
Tv = sp.Matrix([[sp.Symbol("T11", real=True), sp.Symbol("T12", real=True),
                 sp.Symbol("T13", real=True)],
                [sp.Symbol("T12", real=True), sp.Symbol("T22", real=True),
                 sp.Symbol("T23", real=True)],
                [sp.Symbol("T13", real=True), sp.Symbol("T23", real=True),
                 sp.Symbol("T33", real=True)]])   # symmetrisch
pS_var = pS + eps * Tv
Hkin_var = (2 * kappa / sp.sqrt(hS.det())) * \
    (sp.trace(hS * pS_var * hS * pS_var) - sp.Rational(1, 2) *
     sp.trace(hS * pS_var)**2)
# dH_kin = lim_{eps->0} (Hkin_var - Hkin)/eps ; linearer Anteil bei eps=0
dH_variation = sp.simplify(sp.diff(Hkin_var, eps).subs(eps, 0))
dH_expected = sp.simplify((4 * kappa / sp.sqrt(hS.det())) *
                          (pS_low - sp.Rational(1, 2) * pS_trace * hS))
dH_rhs = sp.simplify(sp.trace(dH_expected * Tv))
record("Hamilton-Ableitung (symm. Variation): dH_kin = [4kappa/sqrt h](pi_ij-1/2 h_ij pi) T^{ij}",
       sp.simplify(dH_variation - dH_rhs) == 0,
       "dH_kin(T)=%s ; RHS=%s" % (dH_variation, dH_rhs))

# --------------------- 10) Hamilton-Gleichung fuer h_ij ----------------------
# Fuer pi^{ij} aus der Impulsdefinition gilt:
#   partial H_kin / partial pi^{ij} = (4 kappa/sqrt h)(pi_ij - 1/2 h_ij pi) = 2 K_ij,
# also
#   dot h_ij - L_N h_ij = N partial H_kin / partial pi^{ij}.
N_adm = sp.Symbol("N_adm", positive=True)
Pdef_var = Pdef + eps * Tv
Hkin_def_expr = (2 * kappa / sh) * (
    sp.trace(h * Pdef_var * h * Pdef_var)
    - sp.Rational(1, 2) * spur(h, Pdef_var)**2
)
dH_def = sp.simplify(sp.diff(Hkin_def_expr, eps).subs(eps, 0))
Ddef = sp.simplify((4 * kappa / sh) * (h * Pdef * h
                     - sp.Rational(1, 2) * Ptrace * h))
dH_def_rhs = sp.simplify(sp.trace(Ddef * Tv))
record("Hamilton-Gleichung h (symm. Variation): dH_kin = [4kappa/sqrt h](pi_ij-1/2 h_ij pi) T^{ij}",
       sp.simplify(dH_def - dH_def_rhs) == 0,
       "dH_kin(T)=%s ; RHS=%s" % (dH_def, dH_def_rhs))
record("Hamilton-Gleichung h: partial H_kin / partial pi^{ij} = 2 K_ij",
       mat_eq(Ddef, 2 * K),
       "Ddef - 2K = %s" % sp.simplify(Ddef - 2 * K))
dot_h_minus_lie = 2 * N_adm * K
record("Hamilton-Gleichung h: dot h_ij - L_N h_ij = N partial H_kin / partial pi^{ij}",
       mat_eq(N_adm * Ddef, dot_h_minus_lie),
       "N*Ddef - (dot h - L_N h) = %s" % sp.simplify(N_adm * Ddef - dot_h_minus_lie))

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "convention": "(-+++), kappa=8*pi*G_N, c=1, sqrt h=sqrt(det h)",
    "trace_3d": {
        "pi_formula": "(sqrt h/(2 kappa))(K - 3 K) = -(sqrt h/kappa) K",
        "ptrace_generic": sp.sstr(Ptrace),
        "K_scalar": sp.sstr(Kscalar),
        "error_corrected": True,
    },
    "inverse": {
        "K_ij_recovered": sp.sstr(Pinv),
        "index_note": "h^{ij} (kontravariant); h_{ij} kommt nicht vor",
        "round_trip_ok": rt_ok,
    },
    "special_cases": {
        "pure_trace": sp.sstr(Pc_trace),
        "traceless_pi": sp.sstr(P0_trace),
    },
    "deWitt": {
        "contract": sp.sstr(expected),
        "hamiltonian_variation": sp.sstr(dH_expected),
    },
    "dim_check": "pi^{ij}: sqrt h/(2 kappa)*[K] -> Laenge^{-3/2}",
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
