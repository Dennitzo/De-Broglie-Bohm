# Phase 4, Schritt P-KORR: Dirac-Bilanz des getrimmten konformen Modells in einheitlicher
# orthonormaler (verdoppelter) Impulsnormierung, mit alpha_c als dualer Cos-Geometriemode.
#
# Konfiguration:  alpha = a0 + a1 sin(th) + ac cos(th),  phi = c + q sin(th)
# Impulse:        P = P0 + 2 P1 sin(th) + 2 Pc cos(th),  p = p0 + 2 pq sin(th)
# Kanonische Paare: (a0,P0), (a1,P1), (ac,Pc), (c,p0), (q,pq)   [Paarung in (1) exakt nachgerechnet]
# Modelldichte (Ansatz-Ritz-Projektion, identisch zu work/phase4_Hbar_check.py):
#   H = (beta/2) e^{-3a} (-P^2 + p^2) - e^{a} (a')^2/(6 beta) + e^{a} (phi')^2/2 + (m^2/2) e^{3a} phi^2
#   Hbar = <H>,  <f> = (2pi)^-1 int_0^{2pi} f dth
# Constraints:    Cs = <r sin>, Cc = <r cos> mit r = P'/3 - p phi'
# Schnellausfuehrung: keine globalen simplify-Aufrufe, Poly nur mit domain='EX'.
import json
import sys

import sympy as sp

OUT = "work/phase4_dirac_v2.json"
out = {}


def dump(tag):
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, sort_keys=True)
    print("### done", tag, flush=True)


th = sp.symbols('theta', real=True)
a0, a1, ac, c, q = sp.symbols('alpha0 alpha1 alpha_c c q', real=True)
P0, P1, Pc, p0, pq = sp.symbols('P0 P1 Pc p0 pq', real=True)
be, m = sp.symbols('beta m', positive=True)
a0d, a1d, acd, cd, qd = sp.symbols('adot0 adot1 adotc cdot qdot', real=True)


def mean(f):
    return sp.expand(sp.integrate(sp.expand(f), (th, 0, 2 * sp.pi)) / (2 * sp.pi))


alpha = a0 + a1 * sp.sin(th) + ac * sp.cos(th)
phi = c + q * sp.sin(th)
P = P0 + 2 * P1 * sp.sin(th) + 2 * Pc * sp.cos(th)
p = p0 + 2 * pq * sp.sin(th)

# (1) Paarung
pairG = sp.expand(mean(P * (a0d + a1d * sp.sin(th) + acd * sp.cos(th))))
pairM = sp.expand(mean(p * (cd + qd * sp.sin(th))))
out["pairing_geom_diff"] = str(sp.expand(pairG - (P0 * a0d + P1 * a1d + Pc * acd)))
out["pairing_matter_diff"] = str(sp.expand(pairM - (p0 * cd + pq * qd)))

# (2) Impulsconstraint-Projektionen, (2pi)-Mittel
r = sp.expand(sp.Rational(1, 3) * sp.diff(P, th) - p * sp.diff(phi, th))
Cs = mean(r * sp.sin(th))
Cc = mean(r * sp.cos(th))
out["residual"] = str(r)
out["mean_r"] = str(mean(r))
out["mean_r_sin"] = str(Cs)
out["mean_r_cos"] = str(Cc)
out["sin2_coeff"] = str(sp.simplify(2 * mean(r * sp.sin(2 * th))))
out["residual_minus_proj"] = str(sp.expand(r - (Cs * sp.sin(th) + Cc * sp.cos(th))))
dump("1-2")

# (3) geschlossene Form von Hbar fuer ac = 0 (inkl. Pc^2-Anteil)
I0 = lambda z: sp.besseli(0, z)
I1 = lambda z: sp.besseli(1, z)
I2 = lambda z: sp.besseli(2, z)
M0m, M1m = I0(-3 * a1), I1(-3 * a1)
M2m = (I0(-3 * a1) + I2(-3 * a1)) / 2
C2m = (I0(-3 * a1) - I2(-3 * a1)) / 2
C2a = (I0(a1) - I2(a1)) / 2
M0p, M1p = I0(3 * a1), I1(3 * a1)
M2p = (I0(3 * a1) + I2(3 * a1)) / 2
Hbar5 = sp.expand((be / 2) * sp.exp(-3 * a0) * (
    -(P0 ** 2 * M0m + 4 * P0 * P1 * M1m + 4 * P1 ** 2 * M2m + 4 * Pc ** 2 * C2m)
    + (p0 ** 2 * M0m + 4 * p0 * pq * M1m + 4 * pq ** 2 * M2m))
    - (sp.exp(a0) / (6 * be)) * a1 ** 2 * C2a + (sp.exp(a0) / 2) * q ** 2 * C2a
    + (m ** 2 / 2) * sp.exp(3 * a0) * (c ** 2 * M0p + 2 * c * q * M1p + q ** 2 * M2p))

# On-Flaehe Cs=Cc=0 (Pc=0, P1 = 3 p0 q / 2) vs. gedruckter §19.7-Ausdruck
Hbar_on = sp.expand(Hbar5.subs({Pc: 0, P1: sp.Rational(3, 2) * p0 * q}))
printed = sp.expand((be / 2) * sp.exp(-3 * a0) * (
    -(P0 ** 2 * M0m + 6 * P0 * p0 * q * M1m + 9 * p0 ** 2 * q ** 2 * M2m)
    + p0 ** 2 * M0m + 4 * p0 * pq * M1m + 4 * pq ** 2 * M2m)
    - (sp.exp(a0) / (6 * be)) * a1 ** 2 * C2a + (sp.exp(a0) / 2) * q ** 2 * C2a
    + (m ** 2 / 2) * sp.exp(3 * a0) * (c ** 2 * M0p + 2 * c * q * M1p + q ** 2 * M2p))
out["Hbar_on_vs_printed_diff"] = str(sp.simplify(Hbar_on - printed))
dump("3")

# numerische Gegenprobe: geschlossene Form (ac=0) gegen Direktquadratur
from mpmath import mp, mpf, quad, sin as ms, cos as mc, exp as me, besseli as mI, findroot
mp.dps = 30


def H_direct(av, a1v, acv, cv, qv, P0v, P1v, Pcv, p0v, pqv, bev, mv):
    def dens(t):
        A = av + a1v * ms(t) + acv * mc(t)
        Ap = a1v * mc(t) - acv * ms(t)
        Pph = P0v + 2 * P1v * ms(t) + 2 * Pcv * mc(t)
        pph = p0v + 2 * pqv * ms(t)
        return (bev / 2) * me(-3 * A) * (-Pph ** 2 + pph ** 2) \
            - me(A) * Ap ** 2 / (6 * bev) + me(A) * (qv * mc(t)) ** 2 / 2 \
            + (mv ** 2 / 2) * me(3 * A) * (cv + qv * ms(t)) ** 2
    return quad(dens, [0, 2 * mp.pi]) / (2 * mp.pi)


def H_closed0(av, a1v, cv, qv, P0v, P1v, Pcv, p0v, pqv, bev, mv):
    m0, m1 = mI(0, -3 * a1v), mI(1, -3 * a1v)
    m2 = (mI(0, -3 * a1v) + mI(2, -3 * a1v)) / 2
    c2m = (mI(0, -3 * a1v) - mI(2, -3 * a1v)) / 2
    c2a = (mI(0, a1v) - mI(2, a1v)) / 2
    q0, q1 = mI(0, 3 * a1v), mI(1, 3 * a1v)
    q2 = (mI(0, 3 * a1v) + mI(2, 3 * a1v)) / 2
    kin = (bev / 2) * me(-3 * av) * (
        -(P0v ** 2 * m0 + 4 * P0v * P1v * m1 + 4 * P1v ** 2 * m2 + 4 * Pcv ** 2 * c2m)
        + (p0v ** 2 * m0 + 4 * p0v * pqv * m1 + 4 * pqv ** 2 * m2))
    return kin - me(av) * a1v ** 2 * c2a / (6 * bev) + me(av) * qv ** 2 * c2a / 2 \
        + (mv ** 2 / 2) * me(3 * av) * (cv ** 2 * q0 + 2 * cv * qv * q1 + qv ** 2 * q2)


Wv = [mpf('0.3'), mpf('0.25'), mpf('0.15'), mpf('0.6'), mpf('0.4'),
      mpf('-1.1'), mpf('0.2'), mpf('-0.45'), mpf('0.8'), mpf('-0.35'), mpf('1.7'), mpf('0.9')]
out["Hbar_closed_vs_quadrature_ac0"] = str(
    H_closed0(Wv[0], Wv[1], Wv[3], Wv[4], Wv[5], Wv[6], Wv[7], Wv[8], Wv[9], Wv[10], Wv[11])
    - H_direct(Wv[0], Wv[1], mpf('0'), Wv[3], Wv[4], Wv[5], Wv[6], Wv[7], Wv[8], Wv[9], Wv[10], Wv[11]))
dump("3b")

# (4) Dirac-Ketten (expand statt simplify)
coords = [a0, a1, ac, c, q]
moms = [P0, P1, Pc, p0, pq]


def pb(f, g):
    s = 0
    for x, pm in zip(coords, moms):
        s += sp.diff(f, x) * sp.diff(g, pm) - sp.diff(f, pm) * sp.diff(g, x)
    return sp.expand(s)


out["Cs_Cc"] = str(pb(Cs, Cc))
C2e = pb(Cc, Hbar5)
C3e = pb(Cs, Hbar5)
out["C2_from_Cc_Hbar"] = str(C2e)
out["C3_from_Cs_Hbar"] = str(C3e)
pol2 = sp.Poly(sp.expand(C2e), P0, P1, Pc, p0, pq, domain='EX')
out["C2_monomials"] = [[str(mm), str(cf)] for mm, cf in zip(pol2.monoms(), pol2.coeffs())]
pol3 = sp.Poly(sp.expand(C3e), P0, P1, Pc, p0, pq, domain='EX')
out["C3_monomials"] = [[str(mm), str(cf)] for mm, cf in zip(pol3.monoms(), pol3.coeffs())]
out["C2_total_degree"] = str(sp.Poly(C2e, P0, P1, Pc, p0, pq, domain='EX').total_degree())
out["C2_on_Cc0_Cs0"] = str(sp.expand(C2e.subs({P1: sp.Rational(3, 2) * p0 * q, Pc: 0})))
dump("4")

# (5) Klammer {Cc,C2}, {Cs,C2}
Kcc = pb(Cc, C2e)
Kcs = pb(Cs, C2e)
out["Cc_C2"] = str(Kcc)
polK = sp.Poly(sp.expand(Kcc), P0, P1, Pc, p0, pq, domain='EX')
out["Cc_C2_monomials"] = [[str(mm), str(cf)] for mm, cf in zip(polK.monoms(), polK.coeffs())]
out["Cs_C2_nonzero"] = str(sp.expand(Kcs) != 0)
out["Cs_C2"] = str(Kcs)
dump("5")

# (6) homogener Zweig q=P1=Pc=pq=0: Reihenentwicklung in a1
C2hom = sp.expand(C2e.subs({q: 0, P1: 0, Pc: 0, pq: 0}))
out["C2_hom_not_identically_zero"] = str(C2hom != 0)
fh = sp.lambdify((a0, a1, c, P0, p0, be, m), C2hom, modules='mpmath')
ser = sp.series(C2hom, a1, 0, 9).removeO()
out["C2_hom_series_a1_upto8"] = str(sp.factor(sp.simplify(ser)))
for pw in range(1, 9):
    term = sp.simplify(ser.coeff(a1, pw))
    if term != 0:
        out["C2_hom_lead_power"] = pw
        out["C2_hom_lead_coeff"] = str(term)
        out["C2_hom_lead_coeff_factored"] = str(sp.factor(sp.together(term)))
        # numerische Gegenpruefung der Fuehrungskoeffizienten
        fc = sp.lambdify((a0, c, P0, p0, be, m), term, modules='mpmath')
        cval = fc(mpf('0.3'), mpf('0.6'), mpf('-1.1'), mpf('0.8'), mpf('1.7'), mpf('0.9'))
        out["C2_hom_lead_coeff_numeric"] = str(cval)
        out["C2_hom_ratio_check"] = str((fh(mpf('0.3'), mpf('1e-4'), mpf('0.6'), mpf('-1.1'), mpf('0.8'), mpf('1.7'), mpf('0.9')) / mpf('1e-4') - cval) / cval)
        break
fh = sp.lambdify((a0, a1, c, P0, p0, be, m), C2hom, modules='mpmath')
sc = []
for s in ['0.1', '0.01', '0.001']:
    v = fh(mpf('0.3'), mpf(s), mpf('0.6'), mpf('-1.1'), mpf('0.8'), mpf('1.7'), mpf('0.9'))
    sc.append([s, str(v)])
out["C2_hom_numeric_scaling"] = sc
out["C2_hom_lead_over_alpha1_numeric"] = str(fh(mpf('0.3'), mpf('0.01'), mpf('0.6'), mpf('-1.1'), mpf('0.8'), mpf('1.7'), mpf('0.9')) / mpf('0.01'))
dump("6")

# (7) Schliessbarkeitstest auf Cs=Cc=Hbar5=0 (ac=0)
P1v = mpf('1.5') * Wv[8] * Wv[4]
Hc = lambda P0v: H_closed0(Wv[0], Wv[1], Wv[3], Wv[4], P0v, P1v, mpf('0'), Wv[8], Wv[9], Wv[10], Wv[11])
P0root = findroot(Hc, mpf('-1.0'))
out["Hbar_root_P0"] = str(P0root)
out["Hbar_root_residual"] = str(Hc(P0root))
fC2 = sp.lambdify((a0, a1, c, q, P0, P1, Pc, p0, pq, be, m), C2e, modules='mpmath')
fC3 = sp.lambdify((a0, a1, c, q, P0, P1, Pc, p0, pq, be, m), C3e, modules='mpmath')
surf = (Wv[0], Wv[1], Wv[3], Wv[4], P0root, P1v, mpf('0'), Wv[8], Wv[9], Wv[10], Wv[11])
out["C2_on_constraint_surface"] = str(fC2(*surf))
out["C3_on_constraint_surface"] = str(fC3(*surf))
# zweiter Flaechenpunkt mit anderem Vorzeichen von P0, um Nichtverschwinden nicht punktuell zu behaupten
P0root2 = findroot(Hc, mpf('1.5'))
surf2 = (Wv[0], Wv[1], Wv[3], Wv[4], P0root2, P1v, mpf('0'), Wv[8], Wv[9], Wv[10], Wv[11])
out["Hbar_root2_P0"] = str(P0root2)
out["Hbar_root2_residual"] = str(Hc(P0root2))
out["C2_on_constraint_surface_2"] = str(fC2(*surf2))
out["Cc_on_surface_residual"] = str(P1v / 3 - Wv[8] * Wv[4] / 2)
dump("7")

# (8) ac != 0: {Cc, Hbar} aus Quadratur, zwei unabhaengige Routen
# Route A: analytisch nach Parameter abgeleitete Integranden, dann Quadratur
# Route B: zentrale Differenzen der Direktquadratur
dens = (be / 2) * sp.exp(-3 * alpha) * (-P ** 2 + p ** 2) \
    - sp.exp(alpha) * sp.diff(alpha, th) ** 2 / (6 * be) \
    + sp.exp(alpha) * sp.diff(phi, th) ** 2 / 2 \
    + (m ** 2 / 2) * sp.exp(3 * alpha) * phi ** 2
varnum = (a0, a1, ac, c, q, P0, P1, Pc, p0, pq, be, m, th)
fda1 = sp.lambdify(varnum, sp.expand(sp.diff(dens, a1)), modules='mpmath')
fdpq = sp.lambdify(varnum, sp.expand(sp.diff(dens, pq)), modules='mpmath')
quadmean = lambda f: quad(f, [0, 2 * mp.pi]) / (2 * mp.pi)
par = Wv[:12]
dHda1_A = quadmean(lambda t: fda1(*(par + [t])))
dHdpq_A = quadmean(lambda t: fdpq(*(par + [t])))
dCcdq = -Wv[8] / 2
dCcdP1 = mpf(1) / 3
routeA = dCcdq * dHdpq_A - dCcdP1 * dHda1_A
h = mpf('1e-5')
def Hnum(vs):
    return H_direct(*vs)
b1 = Wv[:]; b2 = Wv[:]
b1[1] += h; b2[1] -= h
dHda1_B = (Hnum(b1) - Hnum(b2)) / (2 * h)
b1 = Wv[:]; b2 = Wv[:]
b1[9] += h; b2[9] -= h
dHdpq_B = (Hnum(b1) - Hnum(b2)) / (2 * h)
routeB = dCcdq * dHdpq_B - dCcdP1 * dHda1_B
out["ac_nonzero_routeA_CcHbar"] = str(routeA)
out["ac_nonzero_routeB_CcHbar"] = str(routeB)
out["ac_nonzero_routes_rel_diff"] = str(abs(routeA - routeB) / max(abs(routeA), mpf('1e-30')))
out["ac_nonzero_Hbar_quadrature_value"] = str(Hnum(Wv))

# (8b) C3 und {Cs,C2} fuer ac != 0 

fda_c = sp.lambdify(varnum, sp.expand(sp.diff(dens, ac)), modules='mpmath')
fdc = sp.lambdify(varnum, sp.expand(sp.diff(dens, c)), modules='mpmath')
out["ac_nonzero_routeA_components"] = [str(dHda1_A), str(dHdpq_A)]


def C3_quad(vs):
    return -mpf(1) / 3 * quadmean(lambda t: fda_c(*(vs[:12] + [t])))


def C2_quad(vs):
    d_a1 = quadmean(lambda t: fda1(*(vs[:12] + [t])))
    d_pq = quadmean(lambda t: fdpq(*(vs[:12] + [t])))
    d_c = quadmean(lambda t: fdc(*(vs[:12] + [t])))
    return (-vs[8] / 2) * d_pq - (mpf(1) / 3) * d_a1 + (vs[4] / 2) * d_c


out["C3_quad_ac_nonzero"] = str(C3_quad(Wv))
Wz0 = Wv[:]
Wz0[2] = mpf('0')
out["C3_quad_ac_zero"] = str(C3_quad(Wz0))
kk = mpf('1e-4')
w1 = Wv[:]
w2 = Wv[:]
w1[2] += kk
w2[2] -= kk
out["Cs_C2_quad_ac_nonzero"] = str(-(mpf(1) / 3) * (C2_quad(w1) - C2_quad(w2)) / (2 * kk))
z1 = Wz0[:]
z2 = Wz0[:]
z1[2] += kk
z2[2] -= kk
out["Cs_C2_quad_ac_zero"] = str(-(mpf(1) / 3) * (C2_quad(z1) - C2_quad(z2)) / (2 * kk))
out["C2_quad_ac_zero_vs_symbolic"] = str(C2_quad(Wz0))
dump("8b")

# (8c) Anschluss der Quadratur-Route an die symbolische Closed-Form bei ac = 0
Wz = Wv[:]
Wz[2] = mpf('0')
routeA0 = C2_quad(Wz)
fC2num = sp.lambdify((a0, a1, ac, c, q, P0, P1, Pc, p0, pq, be, m), C2e, modules='mpmath')
sym0 = fC2num(*Wz)
out["ac0_routeA_vs_symbolic_C2"] = str(routeA0)
out["ac0_symbolic_C2"] = str(sym0)
out["ac0_routeA_vs_symbolic_reldiff"] = str(abs(routeA0 - sym0) / max(abs(sym0), mpf('1e-30')))
dump("8c")
print("ALL DONE", flush=True)
