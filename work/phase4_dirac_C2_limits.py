# Phase 4, Schritt P2: Limit des Sekundaerconstraints auf dem homogenen Zweig + Klammer {C_c, C2}.
import json
import sympy as sp

a0, a1, ac, c, q = sp.symbols('alpha0 alpha1 alpha_c c q', real=True)
P0, P1, Pc, p0, pq = sp.symbols('P0 P1 Pc p0 pq', real=True)
be, m = sp.symbols('beta m', positive=True)
Ps, p1 = 2*P1, 2*pq
I0 = lambda z: sp.besseli(0, z)
I1 = lambda z: sp.besseli(1, z)
I2 = lambda z: sp.besseli(2, z)
M0m = I0(-3*a1); M1m = I1(-3*a1); M2m = (I0(-3*a1)+I2(-3*a1))/2
C2g = (I0(a1)-I2(a1))/2
M0p = I0(3*a1); M1p = I1(3*a1); M2p = (I0(3*a1)+I2(3*a1))/2
kinG = -(be/2)*sp.exp(-3*a0)*(P0**2*M0m + 2*P0*Ps*M1m + Ps**2*M2m)
kinM = (be/2)*sp.exp(-3*a0)*(p0**2*M0m + 2*p0*p1*M1m + p1**2*M2m)
gA   = -(sp.exp(a0)/(6*be))*a1**2*C2g
gP   = (sp.exp(a0)/2)*q**2*C2g
pot  = (m**2/2)*sp.exp(3*a0)*(c**2*M0p + 2*c*q*M1p + q**2*M2p)
Hbar = kinG + kinM + gA + gP + pot
Cc = sp.Rational(2,3)*P1 - 3*p0*q
coords = [a0, a1, ac, c, q]; mom = [P0, P1, Pc, p0, pq]
def pb(f, g):
    s = 0
    for x, p in zip(coords, mom):
        s += sp.diff(f, x)*sp.diff(g, p) - sp.diff(f, p)*sp.diff(g, x)
    return sp.simplify(s)

res = {}
C2 = sp.simplify(pb(Cc, Hbar))
# (a) homogener Zweig: q=0, P1=0, pq=0, Pc=0; Limit a1 -> 0 (zweiwege pruefen)
C2hom_expr = sp.simplify(C2.subs({q: 0, P1: 0, pq: 0, Pc: 0}))
res["C2hom_expr_simpl"] = str(C2hom_expr)
lim = sp.limit(C2hom_expr, a1, 0)
res["C2hom_limit_a1_0"] = str(lim)
# numerischer Limit fuer endliche P0, p0, c
from mpmath import mp, mpf
mp.dps = 30
F = sp.lambdify((a0, a1, c, P0, p0, be, m), C2hom_expr, modules='mpmath')
vals = [(mpf('0.3'), mpf('1e-4'), mpf('0.6'), mpf('-1.1'), mpf('0.8'), mpf('1.7'), mpf('0.9'))]
res["C2hom_numeric_a1_1e-4"] = str(F(*vals[0]))
res["C2hom_numeric_a1_1e-2"] = str(F(0.3, 1e-2, 0.6, -1.1, 0.8, 1.7, 0.9))
# (b) Klammer {C_c, C2}: zweites-klassig?
K = sp.simplify(pb(Cc, C2))
res["Cc_C2_full"] = str(K)
Kpoly = sp.Poly(sp.expand(K), P0, P1, p0, pq)
res["Cc_C2_nonzero_monomials"] = [str(mm) for mm in Kpoly.monoms()]
res["Cc_C2_zero"] = str(sp.simplify(K) == 0)
print(json.dumps(res, ensure_ascii=False, indent=1))
