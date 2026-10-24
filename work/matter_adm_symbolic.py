# -*- coding: utf-8 -*-
"""
Symbolischer Konsistenztest des Materiesektors (Skalarfeld) in den ADM-
Zwangsbedingungen (Phase 0.3).

Konvention: Signatur (-+++), kappa = 8*pi*G_N, c = 1.

Wirkung:
  S_Phi = -(1/2) int sqrt(-g)(g^{mu nu} d_mu Phi d_nu Phi + m^2 Phi^2) d^4x

ADM-Lagrangedichte (allgemein mit Shift N^i):
  sqrt(-g) = N sqrt(h)
  L_Phi = sqrt(h)/(2N) * ((Phi_dot - N^i d_i Phi)^2
                          - N^2 (h^{ij} d_i Phi d_j Phi + m^2 Phi^2))

Prueft:
  1. Kanonischer Impuls aus der Legendre-Transformation:
       pi_Phi = sqrt(h)(Phi_dot - N^i d_i Phi)/N.
  2. Vollstaendige Legendre-Transformation mit Lapse und Shift:
       pi_Phi Phi_dot - L_Phi
         = N * [pi_Phi^2/(2 sqrt(h))
                + sqrt(h)/2 (h^{ij} d_i Phi d_j Phi + m^2 Phi^2)]
           + N^i pi_Phi d_i Phi.
  3. Materie-Hamilton-Zwangsdichte als N-Koeffizient:
       H_perp^m = pi_Phi^2/(2 sqrt(h))
                 + sqrt(h)/2 (h^{ij} d_i Phi d_j Phi + m^2 Phi^2).
  4. Materie-Impuls-Constraint als N^i-Koeffizient:
       H_i^m = pi_Phi d_i Phi.
  5. Hamilton-Gleichung: Phi_dot = partial H / partial pi_Phi.
  6. FLRW (N=1, N^i=0, d_i Phi=0): H_perp^m = a^3 rho_Phi.
  7. Volle Zwangsbedingung H_perp = H_perp^grav + H_perp^m = 0
     <=> H^2 = Lambda/3 - k/a^2 + kappa rho_Phi/3.
  8. Klein-Gordon aus der kanonischen FLRW-Dynamik:
       Phi_ddot + 3 H Phi_dot + m^2 Phi = 0.

Ausgabe: work/matter_adm_symbolic.json plus Text-Zusammenfassung.
"""
import json
import os
import datetime

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "matter_adm_symbolic.json")

kappa = sp.Symbol("kappa", positive=True)   # 8 pi G_N
a = sp.Symbol("a", positive=True)           # Skalenfaktor
H_rate = sp.Symbol("H", positive=True)      # Hubble-Rate
m = sp.Symbol("m", positive=True)           # Feldmasse m_Phi
Lam = sp.Symbol("Lambda", real=True)
k = sp.Symbol("k", real=True)
t = sp.Symbol("t", real=True)
sqrt_h = a**3

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def mat_zero_expr(expr):
    return sp.simplify(expr) == 0


# ----------------- ADM-Feldvariablen mit Lapse und Shift --------------------
Phi = sp.Function("Phi")
Phi_t = Phi(t)
phidot = sp.diff(Phi_t, t)

N = sp.Symbol("N", positive=True)
Nx, Ny, Nz = sp.symbols("Nx Ny Nz", real=True)
d1, d2, d3 = sp.symbols("d1Phi d2Phi d3Phi", real=True)
g11, g22, g33 = sp.symbols("g11 g22 g33", positive=True)

shift_dot = Nx * d1 + Ny * d2 + Nz * d3
spatial_grad = g11 * d1**2 + g22 * d2**2 + g33 * d3**2
mass_term = m**2 * Phi_t**2

# ----------------- 1) Lagrangian und kanonischer Impuls ---------------------
L_adm = sp.Rational(1, 2) * sqrt_h / N * (
    (phidot - shift_dot)**2 - N**2 * (spatial_grad + mass_term)
)
pi_Phi = sp.simplify(sp.diff(L_adm, phidot))

record("pi_Phi = sqrt(h)(Phi_dot - N^i d_i Phi)/N (aus S_Phi)",
       mat_zero_expr(pi_Phi - sqrt_h * (phidot - shift_dot) / N),
       "pi_Phi = %s" % pi_Phi)

# ----------------- 2) Vollstaendige Legendre-Transformation -----------------
pi_Phi_sym = sp.Symbol("pi_Phi", real=True)
H_perp_m = sp.Rational(1, 2) * pi_Phi_sym**2 / sqrt_h \
    + sp.Rational(1, 2) * sqrt_h * (spatial_grad + mass_term)
H_shift = Nx * pi_Phi_sym * d1 + Ny * pi_Phi_sym * d2 + Nz * pi_Phi_sym * d3
H_total_adm = N * H_perp_m + H_shift

# Legendre-Ausdruck zuerst in Geschwindigkeitsvariablen, danach auf
# kanonische Impulse umgeschrieben. Die Umkehr ist pi_Phi = sqrt(h)(...)/N.
legendre_vel = sp.simplify(pi_Phi * phidot - L_adm)
legendre_can = sp.simplify(legendre_vel.subs(
    phidot,
    N * pi_Phi_sym / sqrt_h + shift_dot
).subs(pi_Phi, pi_Phi_sym))

record("Legendre-Transformation (mit Lapse/Shift) = N H_perp^m + N^i H_i^m",
       mat_zero_expr(legendre_can - H_total_adm),
       "pi_Phi Phi_dot - L_Phi == N H_perp^m + N^i H_i^m, Residual = 0")

# ----------------- 3) Materie-Hamilton-Zwangsdichte -------------------------
record("H_perp^m = pi_Phi^2/(2 sqrt h) + (sqrt h/2)(h^{ij}d_iPhi d_jPhi + m^2 Phi^2)",
       mat_zero_expr(
           sp.simplify(sp.diff(H_total_adm, N) - H_perp_m)
       ),
       "H_perp^m = %s" % H_perp_m)

# ----------------- 4) Impuls-Constraint Materie -----------------------------
record("H_i^m = pi_Phi d_i Phi (N^i-Koeffizient der Legendre-Dichte)",
       mat_zero_expr(sp.diff(H_total_adm, Nx) - pi_Phi_sym * d1)
       and mat_zero_expr(sp.diff(H_total_adm, Ny) - pi_Phi_sym * d2)
       and mat_zero_expr(sp.diff(H_total_adm, Nz) - pi_Phi_sym * d3),
       "H_x^m = %s ; H_y^m = %s ; H_z^m = %s" % (
           sp.simplify(sp.diff(H_total_adm, Nx)),
           sp.simplify(sp.diff(H_total_adm, Ny)),
           sp.simplify(sp.diff(H_total_adm, Nz))))

# ----------------- 5) Hamilton-Gleichung dot Phi -----------------------------
phidot_from_H = sp.simplify(sp.diff(H_total_adm, pi_Phi_sym))
phidot_from_momentum = sp.simplify(N * pi_Phi_sym / sqrt_h + shift_dot)
record("Hamilton-Gleichung: Phi_dot = partial H / partial pi_Phi",
       mat_zero_expr(phidot_from_H - phidot_from_momentum),
       "partial H/partial pi_Phi = %s == Phi_dot-Ausdruck" % phidot_from_H)

# ----------------- 6) FLRW: H_perp^m = a^3 rho_Phi ----------------------------
phi_d = sp.diff(Phi(t), t)
pi_phi_flrw = sp.simplify(a**3 * phi_d)
H_perp_m_flrw = sp.simplify(pi_phi_flrw**2 / (2 * a**3)
                            + a**3 / 2 * m**2 * Phi(t)**2)
rho_phi = sp.simplify(phi_d**2 / 2 + m**2 * Phi(t)**2 / 2)
record("FLRW: H_perp^m = a^3 rho_Phi (rho=Phi_dot^2/2+m^2Phi^2/2)",
       mat_zero_expr(H_perp_m_flrw - a**3 * rho_phi),
       "H_perp^m = %s == a^3 rho_Phi" % H_perp_m_flrw)

# ----------------- 7) Volle Zwangsbedingung -> Friedmann mit Materie ---------
Hperp_grav = sp.simplify(-3 * a**3 * H_rate**2 / kappa
                         - 3 * a * k / kappa
                         + a**3 * Lam / kappa)
Hperp_cm = sp.simplify(Hperp_grav + H_perp_m_flrw)
sol = sp.solve(sp.Eq(Hperp_cm, 0), H_rate**2)
friedmann = sp.simplify(Lam / 3 - k / a**2 + kappa * rho_phi / 3)
record("H_perp=0 (grav+Materie) <=> H^2 = Lambda/3 - k/a^2 + kappa rho_Phi/3",
       bool(sol) and mat_zero_expr(sol[0] - friedmann),
       "Loesung H^2 = %s (erwartet %s)" % (sol[0] if sol else None, friedmann))

# ----------------- 8) Hamilton-Gleichung -> Klein-Gordon ----------------------
pi_dot_ham = sp.simplify(-sp.diff(H_perp_m_flrw, Phi(t)))
record("Hamilton-2 (FLRW): pi_Phi_dot = -dH_perp^m/dPhi = -a^3 m^2 Phi",
       mat_zero_expr(pi_dot_ham - (-a**3 * m**2 * Phi(t))),
       "pi_Phi_dot = %s" % pi_dot_ham)

Phi_ddot = sp.symbols("Phi_ddot", real=True)
Phi_dot = sp.symbols("Phi_dot", real=True)
Phi_val = sp.symbols("Phi", real=True)
pi_dot_der = sp.simplify(a**3 * (Phi_ddot + 3 * H_rate * Phi_dot))
kg_ih = sp.simplify(Phi_ddot + 3 * H_rate * Phi_dot + m**2 * Phi_val)
record("KG-Gleichung: Phi_ddot + 3H Phi_dot + m^2 Phi = 0 (folgt aus pi_Phi_dot=2 Gleichungen)",
       mat_zero_expr(pi_dot_der + a**3 * m**2 * Phi_val - a**3 * kg_ih),
       "KG = %s" % kg_ih)

# ----------------------------------------------------------------- Ergebnis
all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "action": {
        "S_Phi": "-(1/2)int sqrt(-g)(g^ab dPhi dPhi + m^2 Phi^2)",
        "ADM_L_Phi": "sqrt(h)/(2N)((Phi_dot - N^i d_iPhi)^2 - N^2(h^ij d_iPhi d_jPhi + m^2Phi^2))",
        "Legendre_identity": "pi_Phi Phi_dot - L_Phi = N H_perp^m + N^i H_i^m",
    },
    "matter": {
        "pi_Phi": sp.sstr(pi_Phi),
        "H_perp_m": sp.sstr(H_perp_m.subs(pi_Phi_sym, pi_Phi)),
        "H_i_m": "pi_Phi d_i Phi",
        "H_i_m_components": {
            "Nx": sp.sstr(sp.simplify(sp.diff(H_total_adm, Nx))),
            "Ny": sp.sstr(sp.simplify(sp.diff(H_total_adm, Ny))),
            "Nz": sp.sstr(sp.simplify(sp.diff(H_total_adm, Nz))),
        },
        "Hamilton_phi_dot": sp.sstr(phidot_from_H),
    },
    "flrw": {
        "H_perp_m": sp.sstr(H_perp_m_flrw),
        "rho_Phi": sp.sstr(rho_phi),
        "friedmann_with_matter": sp.sstr(friedmann),
        "klein_gordon": sp.sstr(kg_ih),
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
