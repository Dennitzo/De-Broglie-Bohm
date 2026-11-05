# Gegenpruefung Phase 3, Iteration B: Herleitung der konformen Koeffizienten.
# Prift: (A1) konforme Legendre-Kette: Symplektik-Impuls P = Sum pi*hdot/alphad,
# Spur tau = h_ij pi^ij, Identitaet P = 2 tau; Hamilton-Kinetik in P.
# (A2) Projektion der ADM-Impulsconstraint auf den reinen Spurimpuls mit
# ausgerechneten Christoffelsymbolen: -2 D_j (tau/3 delta_i^j) = -(2/3) d_i tau.
# (A3) Materieanteil: Pi_Phi d_i Phi = p d_i phi in Modellnormierung (Pi_phi = beta*p),
# und beta^{-1} p d_i phi in dichtarisierter Normierung (p = Pi_phi).
# (N1) numerisch: korrigierte Modenloesung (1/3)P' - p0*phi' = 0 fuer
# P = P0 + 3 p0 q sin(kx). (N2) Negativkontrolle: alte Form mit 1/beta verletzten.
import json
import sympy as sp

out = {}

# ---------- A1: konforme Legendre-Kette ----------
alpha, alphad, Nsym = sp.symbols('alpha alphad N', real=True)
c, be, P = sp.symbols('c beta', positive=True), sp.symbols('beta', positive=True), sp.symbols('P', real=True)
c = c[0] if isinstance(c, tuple) else c
c = sp.Symbol('c', positive=True)
be = sp.Symbol('beta', positive=True)
# Modell-Lagrangedichte (homo-Teil des konformen Ansatzes, Konvention §9.2/§11.1):
# L = -6 c e^{3 alpha} alphad^2 / N ; Kette: alphaad = e^{-2a} trace(hdot)/6
Lkin = -6*c*sp.exp(3*alpha)*alphad**2/Nsym
dalphad_dhdot = sp.exp(-2*alpha)/6  # Faktor pro delta_ij-Anteil
pi_comp = sp.diff(Lkin, alphad)*dalphad_dhdot          # pi^ij = pi_comp * delta^ij
Pv = sp.simplify(pi_comp*2*sp.exp(2*alpha)*3)          # Sum_ij pi^ij hdot_ij / alphad
tau = sp.simplify(3*sp.exp(2*alpha)*pi_comp)           # h_ij pi^ij
sols = sp.solve(sp.Eq(Pv, P), alphad)
Hkin = sp.expand((Pv*alphad - Lkin).subs(alphad, sols[0]))
a1 = {
    "P_minus_2tau": sp.simplify(Pv - 2*tau),
    "Hkin_minus_generic": sp.simplify(Hkin - (-Nsym*P**2*sp.exp(-3*alpha)/(24*c))),
    "Hkin_minus_model": sp.simplify(Hkin.subs(c, 1/(12*be)) - (-Nsym*be*P**2*sp.exp(-3*alpha)/2)),
    "P_hom_eq_p_alpha": sp.simplify(Pv.subs(c, 1/(12*be)) - (-sp.exp(3*alpha)*alphad/(be*Nsym))),
}
out["A1_legendre"] = {k: str(v) for k, v in a1.items()}

# ---------- A2: Impulsconstraint-Projektion, echte Christoffelsymbole ----------
x0, x1, x2 = sp.symbols('x0 x1 x2', real=True)
xs = [x0, x1, x2]
A = sp.Function('alpha')(*xs)
T = sp.Function('tau')(*xs)
KD = sp.KroneckerDelta

def Gamma(k, i, j):
    # Gamma^k_ij = delta^k_i d_j alpha + delta^k_j d_i alpha - delta_ij delta^{kl} d_l alpha
    return (KD(k, i)*sp.diff(A, xs[j]) + KD(k, j)*sp.diff(A, xs[i])
            - KD(i, j)*sp.diff(A, xs[k]))

def DV(i):
    # D_j V_i^j fuer V_i^j = (tau/3) delta_i^j
    term1 = sum(sp.diff(T/3*KD(i, j), xs[j]) for j in range(3))
    term2 = sum(Gamma(j, j, m)*(T/3)*KD(i, m) for j in range(3) for m in range(3))
    term3 = sum(Gamma(m, j, i)*(T/3)*KD(m, j) for j in range(3) for m in range(3))
    return sp.simplify(term1 + term2 - term3)

a2 = [sp.simplify(-2*DV(i) + sp.Rational(2, 3)*sp.diff(T, xs[i])) for i in range(3)]
out["A2_projection"] = [str(v) for v in a2]

# ---------- A3: Materie-Reskalierung ----------
Pi_phi, dphi, p_m, beta_s = sp.symbols('Pi_phi dphi p beta', positive=True)
# Phi = phi/sqrt(beta), Pi_Phi = Pi_phi/sqrt(beta): Pi_Phi dPhi = Pi_phi dphi / beta
expr_model = sp.simplify((Pi_phi/sp.sqrt(beta_s))*(dphi/sp.sqrt(beta_s)) - p_m*dphi)
expr_model = expr_model.subs(Pi_phi, beta_s*p_m)
expr_dense = sp.simplify((Pi_phi/sp.sqrt(beta_s))*(dphi/sp.sqrt(beta_s)) - Pi_phi*dphi/beta_s)
out["A3_matter"] = {"model_norm": str(expr_model), "dense_norm": str(expr_dense)}

# ---------- N1/N2: numerisch mit mpmath ----------
from mpmath import mp, mpf, sin, cos, pi as mp_pi
mp.dps = 50

def resid_correct(p0, q, k, x):
    # P(x) = P0 + 3 p0 q sin(kx)  ->  (1/3) P'(x) - p0 phi'(x), phi = q sin(kx)
    return mpf(3)*p0*q*k*cos(k*x)/mpf(3) - p0*q*k*cos(k*x)

def resid_old(p0, q, k, x, b):
    # alte Form (1/3)P' - beta^{-1} p0 phi' mit korrigiertem P (ohne 1/beta)
    return p0*q*k*cos(k*x) - p0*q*k*cos(k*x)/b

probes = [(mpf('0.7'), mpf('1.3'), mpf('1'), mpf('0.3')),
          (mpf('-2.1'), mpf('0.9'), mpf('2'), mpf('1.7')),
          (mpf('0.5'), mpf('-1.25'), mpf('3'), mpf('0.11'))]
out["N1_corrected_residuals"] = [str(resid_correct(*pr)) for pr in probes]
out["N2_negative_old_form"] = [str(resid_old(pr[0], pr[1], pr[2], pr[3], mpf('2.0'))) for pr in probes]

ok = (str(a1["P_minus_2tau"]) == "0"
      and str(a1["Hkin_minus_generic"]) == "0"
      and str(a1["Hkin_minus_model"]) == "0"
      and str(a1["P_hom_eq_p_alpha"]) == "0"
      and all(v == "0" for v in out["A2_projection"])
      and str(expr_model) == "0" and str(expr_dense) == "0")
out["ALL_SYMBOLIC_ZERO"] = bool(ok)
print(json.dumps(out, ensure_ascii=False, indent=1))
