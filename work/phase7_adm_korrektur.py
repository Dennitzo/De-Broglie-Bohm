# -*- coding: utf-8 -*-
"""
Exakte symbolische Gegentraeger der korrigierten konformen ADM-Projektion
(Phase 1) und der davon betroffenen Moden-, Dirac- und Rangrechnungen
(Phase 2) zu Kapitel 19.6 - 19.10.

Keine Simulation, keine Plots: ausschliesslich exakte Polynom-, Fourier- und
kanonische Algebra mit rationalen Koeffizienten.

Konventionen (Manuskript Kap. 5 und 11.1):
  Signatur (-,+,+,+), kappa = 8*pi*G_N, beta = kappa/6, c = 1,
  h_ij = exp(2 alpha) delta_ij, sqrt h = exp(3 alpha),
  P := 2 h_ij pi^ij  =>  pi_i^j = P delta_i^j/6,  pi^ij = (P/6) e^{-2a} delta^ij,
  phi = sqrt(beta) Phi,  p := Pi_phi/beta  =>  Pi_Phi = sqrt(beta) p.
pi^ij ist Dichte vom Gewicht +1; die kovariante Divergenz von pi_i^j enthaelt
deshalb den Gewichtsterm -Gamma^m_{jm} pi_i^j.

Ergenispapier: work/phase7_adm_korrektur.json
"""
import datetime
import itertools
import json
import os

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "phase7_adm_korrektur.json")
CHECKS = {}


def rec(name, ok, detail=""):
    CHECKS[name] = {"pass": bool(ok), "detail": str(detail)}
    print(("PASS " if ok else "FAIL "), "|", name, "|", detail)


def all_ok(items):
    return all(items)


# =====================================================================
# 0) Konforme Geometrie: Christoffelsymbole und Ricci-Skalar
# =====================================================================
x, y, z = sp.symbols("x y z", real=True)
XS = (x, y, z)
al = sp.Function("alpha")(*XS)
Pl = sp.Function("P")(*XS)
fll = sp.Function("phi")(*XS)
pml = sp.Function("p")(*XS)
Nsh = [sp.Function("N%d" % k)(*XS) for k in range(3)]   # getrennte Shift-Komponenten

dd = lambda f, i: sp.diff(f, XS[i])

h_lo = sp.diag(*[sp.exp(2 * al)] * 3)
h_up = h_lo.inv()

G = {}
for k in range(3):
    for i in range(3):
        for j in range(3):
            s = sum((dd(h_lo[j, m], i) + dd(h_lo[i, m], j) - dd(h_lo[i, j], m)) * h_up[k, m]
                    for m in range(3))
            G[(k, i, j)] = sp.simplify(s / 2)

# Vollstaendiger Abgleich ALLER 27 Komponenten gegen die geschlossene Form
# Gamma^k_ij = delta_{ki} a_,j + delta_{kj} a_,i - delta_{ij} a_,k   (euklidisch,
# Indizes 0,1,2 = x,y,z; hoch und tief fallen hier zusammen).
Gform_ok, Gform_bad = True, []
for k in range(3):
    for i in range(3):
        for j in range(3):
            expv = ((1 if k == i else 0) * dd(al, j) + (1 if k == j else 0) * dd(al, i)
                    - (1 if i == j else 0) * dd(al, k))
            if sp.simplify(G[(k, i, j)] - expv) != 0:
                Gform_ok = False
                Gform_bad.append((k, i, j))
rec("Christoffel konform: alle 27 Komponenten erfuellen "
    "Gamma^k_ij = delta_ki a_,j + delta_kj a_,i - delta_ij a_,k",
    Gform_ok and sp.simplify(G[(0, 0, 0)] - dd(al, 0)) == 0
    and sp.simplify(G[(0, 1, 1)] + dd(al, 0)) == 0 and sp.simplify(G[(0, 1, 2)]) == 0,
    "abweichende Komponenten: %s ; Stichproben Gamma^x_xx=a_,x ; Gamma^x_yy=-a_,x ; Gamma^x_yz=0" % Gform_bad)

Gc = [sp.simplify(sum(G[(j, j, i)] for j in range(3))) for i in range(3)]
rec("Kontraktion Gamma^j_ij = 3 a_,i (n = 3)",
    all_ok([sp.simplify(Gc[i] - 3 * dd(al, i)) == 0 for i in range(3)]),
    str(Gc))

lap = sum(sp.diff(al, XS[i], XS[i]) for i in range(3))
g2 = sum(dd(al, i) ** 2 for i in range(3))

# R^i_{jkl} = d_k Gamma^i_lj - d_l Gamma^i_kj + Gamma^i_km Gamma^m_lj - Gamma^i_lm Gamma^m_kj
Rie = {}
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                s = dd(G[(i, l, j)], k) - dd(G[(i, k, j)], l)
                for m in range(3):
                    s = s + G[(i, k, m)] * G[(m, l, j)] - G[(i, l, m)] * G[(m, k, j)]
                Rie[(i, j, k, l)] = sp.simplify(s)
Ric = sp.Matrix(3, 3, lambda i, j: sp.simplify(sum(Rie[(k, i, k, j)] for k in range(3))))
Rscal = sp.expand(sum(h_up[i, j] * Ric[i, j] for i in range(3) for j in range(3)))
Rform = sp.simplify(-sp.exp(-2 * al) * (4 * lap + 2 * g2))
rec("Ricci-Skalar konform: R = -e^{-2a}(4 Lap a + 2|grad a|^2)",
    sp.simplify(Rscal - Rform) == 0, "R = %s" % Rform)

# =====================================================================
# 1) Weg A: Impulsconstraint aus der kovarianten Divergenz (Gewicht +1)
# =====================================================================
pi_up = {(i, j): (sp.simplify(Pl / 6 * h_up[i, j]) if i == j else sp.Integer(0))
         for i in range(3) for j in range(3)}
pmix = sp.Matrix(3, 3, lambda i, j: sp.simplify(sum(h_lo[i, m] * pi_up[(m, j)] for m in range(3))))
rec("pi_i^j = (P/6) delta_i^j",
    all_ok([sp.simplify(pmix[i, j] - (Pl / 6 if i == j else 0)) == 0
            for i in range(3) for j in range(3)]),
    "P = 2 h_ij pi^ij und pi^ij rein spurhaft")

DivW, DivN = {}, {}
for i in range(3):
    sw = sn = sp.Integer(0)
    for j in range(3):
        sw += dd(pmix[i, j], j)
        sn += dd(pmix[i, j], j)
        for m in range(3):
            # gewichtete Kovariante (Dichte +1): + Gamma^j_jm pi_i^m - Gamma^m_ji pi_m^j
            sw += G[(j, j, m)] * pmix[i, m] - G[(m, j, i)] * pmix[m, j]
            # ungewichtete (uebliche) Kovariante: laesst den Gewichtsterm weg
            sn += G[(j, j, m)] * pmix[i, m] - G[(m, j, i)] * pmix[m, j]
        # Dichtegewicht +1: - Gamma^m_{jm} pi_i^j mit CONTRACTEDem Ableitungsindex j
        sw -= Gc[j] * pmix[i, j]
    DivW[i] = sp.simplify(sp.expand(sw))
    DivN[i] = sp.simplify(sp.expand(sn))

tgt = {i: sp.simplify((dd(Pl, i) - 3 * Pl * dd(al, i)) / 6) for i in range(3)}
rec("Weg A korrigiert: D_j pi_i^j = (d_i P - 3 P d_i a)/6  (komponentenweise i=1,2,3)",
    all_ok([sp.simplify(DivW[i] - tgt[i]) == 0 for i in range(3)]), str(DivW))
rec("Weg A ungewichtete Kovariante (bisherige Manuskriptfassung): D_j pi_i^j = d_i P/6",
    all_ok([sp.simplify(DivN[i] - dd(Pl, i) / 6) == 0 for i in range(3)]), str(DivN))

# =====================================================================
# 2) Weg B: Lie-Ableitungs-Funktional int (pi^ij L_N h_ij + p L_N phi)
# =====================================================================
def lie_h(i, j):
    s = sum(Nsh[k] * dd(h_lo[i, j], k) for k in range(3))
    for k in range(3):
        s += h_lo[k, j] * dd(Nsh[k], i) + h_lo[i, k] * dd(Nsh[k], j)
    return sp.expand(s)


integ = sp.Integer(0)
for i in range(3):
    for j in range(3):
        integ += pi_up[(i, j)] * lie_h(i, j)
integ += pml * sum(Nsh[k] * dd(fll, k) for k in range(3))
integ = sp.expand(integ)

HB = {}
for i in range(3):
    # Koeffizient von N^i selbst und von d_j N^i (Integrand ist linear im Shift)
    Lc = sp.simplify(sp.expand(integ).coeff(Nsh[i]))
    Msum = sp.Integer(0)
    for j in range(3):
        coef = sp.simplify(sp.expand(integ).coeff(dd(Nsh[i], j)))
        Msum += dd(coef, j)
    HB[i] = sp.simplify(sp.expand(Lc - Msum))

Hcorr = {i: sp.simplify(-2 * DivW[i] + pml * dd(fll, i)) for i in range(3)}
Hold = {i: sp.simplify(-2 * DivN[i] + pml * dd(fll, i)) for i in range(3)}
rec("Weg B == Weg A (korrigiert): H_i = -d_i P/3 + P d_i a + p d_i phi",
    all_ok([sp.simplify(HB[i] - Hcorr[i]) == 0 for i in range(3)]), str(HB))
rec("Weg B minus alte Fassung = + P d_i a (die alte Constraint fehlt genau dieses Glied)",
    all_ok([sp.simplify(HB[i] - Hold[i] - Pl * dd(al, i)) == 0 for i in range(3)]),
    str([sp.simplify(HB[i] - Hold[i]) for i in range(3)]))

# =====================================================================
# 3) Shift-Anteil der kanonischen Bewegungsgleichungen
# =====================================================================
def func_delta(F, var):
    """delta F / delta var fuer linear in var und d var abhaengiges F (Torus, randfrei)."""
    eta = sp.Function("eta")(*XS)
    eps = sp.Symbol("eps")
    dF = sp.expand(sp.diff(F.subs(var, var + eps * eta), eps).subs(eps, 0))
    out = sp.Integer(0)
    rest = dF
    for i in range(3):
        co = rest.coeff(dd(eta, i))
        if co != 0:
            out -= dd(co, i)
            rest = sp.expand(rest - co * dd(eta, i))
    out += rest.coeff(eta)
    return sp.simplify(sp.expand(out))


Fs = sum(Nsh[i] * Hcorr[i] for i in range(3))
want_a = sum(Nsh[i] * dd(al, i) for i in range(3)) + sum(dd(Nsh[i], i) for i in range(3)) / 3
want_p = sum(Nsh[i] * dd(fll, i) for i in range(3))
got_a = func_delta(Fs, Pl)
got_p = func_delta(Fs, pml)
rec("Shift-Anteil alpha: {alpha, int N^i H_i} = N^i d_i alpha + (d_i N^i)/3",
    sp.simplify(got_a - want_a) == 0, got_a)
rec("Shift-Anteil phi: {phi, int N^i H_i} = N^i d_i phi",
    sp.simplify(got_p - want_p) == 0, got_p)
got_a_old = func_delta(sum(Nsh[i] * Hold[i] for i in range(3)), Pl)
rec("Alte Constraint liefert fuer alpha nur (d_i N^i)/3; Advektion N^i d_i alpha fehlt",
    sp.simplify(got_a_old - sum(dd(Nsh[i], i) for i in range(3)) / 3) == 0, got_a_old)

# =====================================================================
# 4) Negativprobe: konstantes P != 0, konstantes phi, inhomogenes alpha
# =====================================================================
a1f = sp.Function("alpha_inhom")(*XS)
Pc0 = sp.Symbol("P_const")
probe_new = [sp.simplify(Hcorr[i].subs({al: a1f, Pl: Pc0, fll: sp.Symbol("phi_c"), pml: sp.Symbol("p_c")}))
             for i in range(3)]
probe_old = [sp.simplify(Hold[i].subs({al: a1f, Pl: Pc0, fll: sp.Symbol("phi_c"), pml: sp.Symbol("p_c")}))
             for i in range(3)]
rec("Negativprobe (P konst != 0, phi konst, alpha inhomogen): alt = 0, neu = P d_i alpha",
    all_ok([probe_old[i] == 0 for i in range(3)])
    and all_ok([sp.simplify(probe_new[i] - Pc0 * dd(a1f, i)) == 0 for i in range(3)]),
    "alt = %s ; neu = %s" % (probe_old, probe_new))

# =====================================================================
# 5) Kruemmungsterm bei allgemeinem N: partielle Integration auf T^3
# =====================================================================
kappa = sp.symbols("kappa", positive=True)
Nvar = sp.Function("Nfun")(*XS)
cur = sp.expand(Nvar * sp.exp(al) * (4 * lap + 2 * g2) / (2 * kappa))
tot = sp.expand(sum(dd(4 * Nvar * sp.exp(al) * dd(al, i), i) for i in range(3)) / (2 * kappa))
den1 = sp.expand(-Nvar * sp.exp(al) * g2 / kappa)
den2 = sp.expand(-2 * sp.exp(al) * sum(dd(Nvar, i) * dd(al, i) for i in range(3)) / kappa)
rec("Kruemmung bei allgemeinem N: int N e^a(4 Lap a + 2|grad a|^2)/(2 kappa) auf T^3 "
    "= int [-N e^a|grad a|^2/kappa - 2 e^a grad N.grad a/kappa]",
    sp.simplify(cur - tot - den1 - den2) == 0,
    "Totdiv = d_i(4 N e^a d^i a)/(2 kappa) ; Zielglieder %s und %s" % (den1, den2))
Nconst = sp.Symbol("N_t", real=True)
mixed = sp.simplify(sum(sp.diff(Nconst, XS[i]) * dd(al, i) for i in range(3)))
rec("Raeumlich konstante Lapse laesst den gemischten Beitrag exakt wegfallen",
    sp.simplify(mixed) == 0, "gemischter Beitrag bei N = konst = %s" % mixed)

# =====================================================================
# 6) Materienormierung phi = sqrt(beta) Phi bis zur Modengleichung
# =====================================================================
bet = sp.symbols("beta", positive=True)
Ph = fll / sp.sqrt(bet)
gradADM = sp.simplify(sp.expand(sp.exp(3 * al) / 2 * sum(h_up[i, i] * dd(Ph, i) ** 2 for i in range(3))))
gph = sum(dd(fll, i) ** 2 for i in range(3))
tgt_grad = sp.exp(al) * gph / (2 * bet)
rec("Materiegradient: sqrt(h)/2 h^ij d_i Phi d_j Phi = e^a |grad phi|^2/(2 beta)",
    sp.simplify(gradADM - tgt_grad) == 0, "%s gegen %s" % (gradADM, tgt_grad))
rec("Das ist NOT die gedruckte Dichte e^a|grad phi|^2/2: Faktor 1/beta fehlt dort",
    sp.simplify(tgt_grad - sp.exp(al) * gph / 2) != 0,
    "Differenz = %s" % sp.simplify(tgt_grad - sp.exp(al) * gph / 2))
rec("Impulsumrechnung aus kanonischer Einsform: Pi_Phi dPhi = p dphi mit Pi_Phi = sqrt(beta) p",
    sp.simplify(sp.sqrt(bet) * sp.Symbol("p_m") * sp.diff(fll, x) / sp.sqrt(bet)
                - sp.Symbol("p_m") * sp.diff(fll, x)) == 0,
    "Pi_Phi = Pi_phi/sqrt(beta) = sqrt(beta) p")
one_g = sp.simplify(sum(pi_up[(i, j)] * (2 * sp.Symbol("dalpha") * h_lo[i, j])
                        for i in range(3) for j in range(3)))
rec("Geometrie-Einsform: pi^ij (2 d_alpha h_ij) = P d_alpha",
    sp.simplify(one_g - Pl * sp.Symbol("dalpha")) == 0, one_g)
th = sp.symbols("theta", real=True)
mean = lambda e: sp.simplify(sp.integrate(e, (th, 0, 2 * sp.pi)) / (2 * sp.pi))
q1 = sp.symbols("q", real=True)
a0s = sp.symbols("a0", real=True)
gd = sp.exp(a0s) / (2 * bet) * (sp.diff(q1 * sp.sin(th), th)) ** 2
rec("Fourierprojektion des Gradientterms bei beta != 1: < e^a (phi')^2/(2 beta) > = e^{a0} q^2/(4 beta)",
    sp.simplify(mean(gd) - sp.exp(a0s) * q1 ** 2 / (4 * bet)) == 0, mean(gd))
gd_beta2 = sp.simplify(gd.subs(bet, 2))
rec("Explizite Gegenprobe beta = 2: Term = e^{a0} q^2 cos^2(theta)/4, Mittel = e^{a0} q^2/8",
    sp.simplify(mean(gd_beta2) - sp.exp(a0s) * q1 ** 2 / 8) == 0, mean(gd_beta2))

# =====================================================================
# 7) K = 1, eps_m = 0: korrigierte Projektionen aus r = P'/3 - P a' - p phi'
# =====================================================================
A, B, cc = sp.symbols("A B c", real=True)
P0, P1, Pc, p0, pq = sp.symbols("P0 P1 Pc p0 pq", real=True)
m2 = sp.symbols("m", real=True)
alphaK = a0s + A * sp.sin(th) + B * sp.cos(th)
phiK = cc + q1 * sp.sin(th)
PK = P0 + 2 * P1 * sp.sin(th) + 2 * Pc * sp.cos(th)
pK = p0 + 2 * pq * sp.sin(th)
D = lambda e: sp.diff(e, th)
r_new = sp.expand(sp.diff(PK, th) / 3 - PK * D(alphaK) - pK * D(phiK))
r_old = sp.expand(sp.diff(PK, th) / 3 - pK * D(phiK))
C0, Cs, Cc = (sp.simplify(mean(r_new)), sp.simplify(mean(r_new * sp.sin(th))),
              sp.simplify(mean(r_new * sp.cos(th))))
rec("K=1 korrigiert: C0 = B P1 - A Pc", sp.simplify(C0 - (B * P1 - A * Pc)) == 0, C0)
rec("K=1 korrigiert: Cs = -Pc/3 + P0 B/2", sp.simplify(Cs - (-Pc / 3 + P0 * B / 2)) == 0, Cs)
rec("K=1 korrigiert: Cc = P1/3 - (P0 A + p0 q)/2",
    sp.simplify(Cc - (P1 / 3 - (P0 * A + p0 * q1) / 2)) == 0, Cc)
C0o, Cso, Cco = (sp.simplify(mean(r_old)), sp.simplify(mean(r_old * sp.sin(th))),
                 sp.simplify(mean(r_old * sp.cos(th))))
rec("Vergleichsmodell mit r = P'/3 - p phi': C0 = 0, Cs = -Pc/3, Cc = P1/3 - p0 q/2",
    C0o == 0 and sp.simplify(Cso + Pc / 3) == 0
    and sp.simplify(Cco - (P1 / 3 - p0 * q1 / 2)) == 0,
    "C0=%s Cs=%s Cc=%s" % (C0o, Cso, Cco))

Qs, Ps_ = [a0s, A, B, cc, q1], [P0, P1, Pc, p0, pq]
pb = lambda f, g: sp.simplify(sp.expand(
    sum(sp.diff(f, qq) * sp.diff(g, pp) - sp.diff(f, pp) * sp.diff(g, qq)
        for qq, pp in zip(Qs, Ps_))))
br0s = pb(C0, Cs)
diff_slice = sp.simplify(br0s - (-Cc - p0 * q1 / 2))
rec("Klammer allgemein: {C0,Cs} = -P1/3 + P0 A/2",
    sp.simplify(br0s - (-P1 / 3 + P0 * A / 2)) == 0, br0s)
rec("{C0,Cs} = -Cc - p0 q/2 gilt IDENTISCH, also auch VOR jeder Einschraenkung auf eine Scheibe",
    sp.simplify(diff_slice) == 0, "Differenz = %s (freie A, B, c, q, P0..pq)" % diff_slice)
rec("Auf der Scheibe B = alpha_c = 0 ist {C0,Cs} = -Cc - p0 q/2 exakt",
    sp.simplify((br0s - (-Cc - p0 * q1 / 2)).subs(B, 0)) == 0, sp.simplify(br0s.subs(B, 0)))
rec("Weitere Primklaemmer: {Cs,Cc} und {C0,Cc}", True,
    "{Cs,Cc} = %s ; {C0,Cc} = %s" % (pb(Cs, Cc), pb(C0, Cc)))

ref = {q1: sp.Rational(1, 3), a0s: 0, A: 0, B: 0, cc: 0, Pc: 0, pq: 0,
       p0: 1, P1: sp.Rational(1, 2), P0: sp.sqrt(5) / 3}
# Der Referenzpunkt des Auftrags ist beta = 1; beta wird getrennt gehalten,
# damit dieselbe Belegung auch mit beta != 1 gerechnet werden kann.
ref1 = dict(ref)
ref1[bet] = 1
ref2 = dict(ref)
ref2[bet] = 2
pv = [sp.simplify(v.subs(ref)) for v in (C0, Cs, Cc)]
rec("Referenzpunkt (beta=1,m=0,U=0,a0=A=B=c=Pc=pq=0,p0=1,q=1/3,P1=1/2,P0=sqrt5/3): C0=Cs=Cc=0",
    all_ok([v == 0 for v in pv]), str(pv))
brref = sp.simplify(br0s.subs(ref))
rec("Referenzpunkt: {C0,Cs} = -1/6, also nicht null", sp.simplify(brref + sp.Rational(1, 6)) == 0, str(brref))

Hdens = (bet * sp.exp(-3 * alphaK) / 2 * (-PK ** 2 + pK ** 2)
         - sp.exp(alphaK) / (6 * bet) * D(alphaK) ** 2
         + sp.exp(alphaK) / (2 * bet) * D(phiK) ** 2
         + m2 ** 2 * sp.exp(3 * alphaK) * phiK ** 2 / (2 * bet))
# Mittelung NACH Einsetzen der Referenzdaten: mit freien A, B waeren es
# Bessel-Integrale (exp(3 A sin theta)) und nicht mehr trigonometrisch.
Hbar_ref = sp.simplify(mean(sp.expand(Hdens.subs(ref1).subs(m2, 0))))
kin = sp.simplify(mean((bet * sp.exp(-3 * alphaK) / 2 * (-PK ** 2 + pK ** 2)).subs(ref1)))
gA = sp.simplify(mean((-sp.exp(alphaK) / (6 * bet) * D(alphaK) ** 2).subs(ref1)))
gp = sp.simplify(mean((sp.exp(alphaK) / (2 * bet) * D(phiK) ** 2).subs(ref1)))
pot = sp.simplify(mean((m2 ** 2 * sp.exp(3 * alphaK) * phiK ** 2 / (2 * bet)).subs(ref1).subs(m2, 0)))
rec("Referenzpunkt: gemittelter Hamilton-Constraint Hbar = 0 (korrigierte Dichte)",
    Hbar_ref == 0, "kin = %s ; geom-Grau %s ; mat-Grau %s ; Pot %s ; Summe = %s"
    % (kin, gA, gp, pot, sp.simplify(kin + gA + gp + pot)))
rec("Postenkontrolle: kin = -1/36 und Materiegradient = +1/36 (Doppelnormierung <P^2> = P0^2 + 2 P1^2)",
    sp.simplify(kin + sp.Rational(1, 36)) == 0 and sp.simplify(gp - sp.Rational(1, 36)) == 0,
    "<P^2> = %s" % sp.simplify(mean((PK ** 2).subs(ref))))
Hbar_beta2 = sp.simplify(mean(sp.expand(Hdens.subs(ref2).subs(m2, 0))))
kin2 = sp.simplify(mean((bet * sp.exp(-3 * alphaK) / 2 * (-PK ** 2 + pK ** 2)).subs(ref2)))
gp2 = sp.simplify(mean((sp.exp(alphaK) / (2 * bet) * D(phiK) ** 2).subs(ref2)))
rec("Gegenprobe beta = 2: derselbe Punkt erfuellt die gemittelte Hamilton-Constraint NICHT mehr "
    "(kin = -1/18, Materiegradient = +1/72, Summe = -1/24)",
    sp.simplify(Hbar_beta2 + sp.Rational(1, 24)) == 0
    and sp.simplify(kin2 + sp.Rational(1, 18)) == 0
    and sp.simplify(gp2 - sp.Rational(1, 72)) == 0,
    "beta=2: kin = %s ; mat-Grau = %s ; Hbar = %s ; Form: -beta/36 + 1/(36 beta)" % (kin2, gp2, Hbar_beta2))

# =====================================================================
# 8) Dirac-Bilanz: starke Nullheit gehoert zur Dirac-Klammer
# =====================================================================
W, Dv, E = sp.symbols("W D E", real=True)
omega = sp.Matrix([[0, W], [-W, 0]])
ainv = omega.inv()
row_Cc = sp.Matrix([[0, W]])              # {C_c, chi_a}, chi = (C_c, D)
col_H = sp.Matrix([[Dv], [E]])             # {chi_a, Hbar} = ({C_c,Hbar}, {D,Hbar}) = (D, E)
corr = (row_Cc * ainv * col_H)[0]
rec("Dirac-Klammer {C_c,Hbar}_D = D - W*(1/W)*D = 0 exakt; {C_c,Hbar} selbst ist D (nicht 0)",
    sp.simplify(corr - Dv) == 0 and sp.simplify(Dv - corr) == 0,
    "Korrekturterm = %s ; {C_c,Hbar}_D = %s" % (sp.simplify(corr), sp.simplify(Dv - corr)))
ac = sp.symbols("alpha_c", real=True)
Hfc = sp.Function("Hfc")(ac)
Cs_old = -Pc / 3
br_cs = sp.simplify(sp.diff(Cs_old, Pc) * 0 - (-sp.Rational(1, 3)) * 0)
br_cs = sp.simplify(-sp.diff(Cs_old, Pc) * sp.diff(Hfc, ac))
rec("Altes Cs = -Pc/3: {Cs,Hbar} = +(1/3) dHbar/d alpha_c (Manuskript 19.7(iii) druckt Minus)",
    sp.simplify(br_cs - sp.diff(Hfc, ac) / 3) == 0, br_cs)
k1_, k2_ = sp.Rational(465194, 10 ** 6), sp.Rational(465242, 10 ** 6)
rel = sp.simplify((k2_ - k1_) / k1_)
rec("Gedruckte Zahlen 0,465194 und 0,465242: relative Abweichung", True,
    "rel = %s = %.6e (Manuskript sagt 1,0e-8)" % (rel, float(rel)))

# =====================================================================
# 9) Alter S4-Block K = (+)_{m=1..K} (m/3) J: Inverse, det, Kondition
# =====================================================================
for K in [1, 2, 3, 4]:
    KK = sp.Matrix.zeros(2 * K, 2 * K)
    for b in range(K):
        mm = sp.Rational(b + 1, 3)
        KK[2 * b, 2 * b + 1] = -mm
        KK[2 * b + 1, 2 * b] = mm
    INV = sp.Matrix.zeros(2 * K, 2 * K)
    for b in range(K):
        mm = sp.Rational(3, b + 1)
        INV[2 * b, 2 * b + 1] = mm          # -(3/m) J  =  [[0, 3/m],[-3/m,0]]
        INV[2 * b + 1, 2 * b] = -mm
    blk_inv_ok = sp.simplify(KK * INV - sp.eye(2 * K)) == sp.zeros(2 * K, 2 * K)
    detK = sp.simplify(KK.det())
    KG = sp.simplify(KK.T * KK)
    eigs = sorted([sp.simplify(e) for e in KG.eigenvals()])
    rec("S4-Block K=(+)_{m=1..K}(m/3)J, K=%d: inverse = (+)(-3/m)J, det = (K!)^2/9^K, cond_2 = K"
        % K,
        blk_inv_ok and detK == sp.simplify(sp.factorial(K) ** 2 / 9 ** K)
        and eigs[0] == sp.Rational(1, 9) and eigs[-1] == sp.Rational(K ** 2, 9),
        "det = %s ; K^T K Eigenwerte = %s (doppelt) ; cond = sqrt(max/min) = %d"
        % (detK, eigs, K))
    Z = sp.Matrix.zeros(2 * K, 2 * K)
    big = sp.Matrix.vstack(sp.Matrix.hstack(Z, KK), sp.Matrix.hstack(-KK, Z))
    sq_ok = sp.simplify(big * big - sp.Matrix.vstack(
        sp.Matrix.hstack(-(KK * KK), Z), sp.Matrix.hstack(Z, -(KK * KK)))) == sp.zeros(4 * K, 4 * K)
    rec("det([[0,K],[-K,0]]) = (det K)^2 durch Blockmultiplikation, K=%d" % K,
        sp.simplify(big.det() - detK ** 2) == 0 and sq_ok,
        "det(big) = %s ; (det K)^2 = %s ; big^2 Blockform bestaetigt: %s"
        % (sp.simplify(big.det()), detK ** 2, sq_ok))

# =====================================================================
# 10) K = 2-Kette (eps_m = 0): Pfaffsche vs. 2x2-Minor, Rang ueber det
# =====================================================================
u1, u2, pi1, pi2 = sp.symbols("u1 u2 pi1 pi2", real=True)
H1, H2, Pi1, Pi2 = sp.symbols("H1 H2 Pi1 Pi2", real=True)
a2_, b2_ = sp.symbols("a2 b2", real=True)
Qf = [a0s, A, B, a2_, b2_, cc, u1, u2]
Pf = [P0, Pi1, H1, Pi2, H2, p0, pi1, pi2]
U = {1: u1, 2: u2}
PI = {1: pi1, 2: pi2}
HP = {1: H1, 2: H2}
PP = {1: Pi1, 2: Pi2}


def Ccos(m):
    return sp.Rational(m, 3) * PP[m] - sp.Rational(m, 2) * p0 * U[m]


def Csin(m):
    s = -sp.Rational(m, 3) * HP[m]
    tot = sp.Integer(0)
    for j in (1, 2):
        lo = abs(j - m)
        ulo = U[lo] if lo in (1, 2) else sp.Integer(0)
        hi = j + m
        uhi = U[hi] if hi in (1, 2) else sp.Integer(0)
        tot += PI[j] * (sp.Integer(lo) * ulo - sp.Integer(hi) * uhi)
    return sp.simplify(s - tot / 2)


elems = [Csin(1), Ccos(1), Csin(2), Ccos(2)]
print("... K=2-Kette: Poisson-Klammern der vier Kettenelemente", flush=True)
omeg = sp.Matrix.zeros(4, 4)
for i in range(4):
    for j in range(4):
        omeg[i, j] = sp.simplify(sp.expand(
            sum(sp.diff(elems[i], qq) * sp.diff(elems[j], pp) - sp.diff(elems[i], pp) * sp.diff(elems[j], qq)
                for qq, pp in zip(Qf, Pf))))
antisym_ok = sp.simplify(omeg + omeg.T) == sp.zeros(4, 4)
pf = sp.simplify(omeg[0, 1] * omeg[2, 3] - omeg[0, 2] * omeg[1, 3] + omeg[0, 3] * omeg[1, 2])
det4 = sp.simplify(omeg.det())
rec("K=2-Kette eps_m=0, Reihenfolge (C_sin1, C_cos1, C_sin2, C_cos2): antisymmetrisch, "
    "Pfaffsche = -p0^2 u1^2/8, det = Pf^2",
    antisym_ok and sp.simplify(pf + p0 ** 2 * u1 ** 2 / 8) == 0 and sp.simplify(det4 - pf ** 2) == 0,
    "Pf = %s ; det4 = %s" % (pf, det4))
minors = {}
for rs in itertools.combinations(range(4), 2):
    for cs in itertools.combinations(range(4), 2):
        dv = sp.simplify(omeg.extract(list(rs), list(cs)).det())
        if dv != 0:
            minors["rows%s_cols%s" % (rs, cs)] = str(dv)
M12_34 = sp.simplify(omeg.extract([0, 1], [2, 3]).det())
M13_24 = sp.simplify(omeg.extract([0, 2], [1, 3]).det())
rec("K=2-Kette: der 2x2-Minor Zeilen (C_sin1,C_cos1) / Spalten (C_sin2,C_cos2) ist +p0^2 u1^2/8 und damit "
    "verschieden von der Pfaffschen -p0^2 u1^2/8 (Unterschied nur das Vorzeichen, nicht der Betrag)",
    sp.simplify(M12_34 - p0 ** 2 * u1 ** 2 / 8) == 0 and sp.simplify(pf + M12_34) == 0,
    "M(12,34) = %s ; M(13,24) = %s ; Pf = %s ; det4 = %s ; nichtverschwindende Minor = %d"
    % (M12_34, M13_24, pf, det4, len(minors)))
beleg = {p0: 4, u1: -2, u2: -2, pi1: 5, pi2: 2}
det4_b = sp.simplify(det4.subs(beleg))
rec("K=2-Kette: Rang 4 folgt aus der separat berechneten 4x4-Determinante (Belegung p0=4,"
    " u=(-2,-2), pi=(5,2))",
    det4_b != 0,
    "det4 = %s ; Bedingung: p0 != 0 und u1 != 0 (Pf = -p0^2 u1^2/8)" % det4_b)
pf_zero = sp.factor(sp.simplify(pf))
rec("Nichtverschwindensbedingung des 4x4-Minors: p0 != 0 und u1 != 0",
    pf_zero == -p0 ** 2 * u1 ** 2 / 8, "Pf faktorisiert = %s" % pf_zero)

# Paarmuster 19.10 A gegen eigene Rechnung
pat = {("(s1,c1)"): sp.simplify(omeg[0, 1]), ("(s1,c2)"): sp.simplify(omeg[0, 3]),
       ("(c1,s2)"): sp.simplify(omeg[1, 2]), ("(s1,s2)"): sp.simplify(omeg[0, 2])}
rec("Paarmuster 19.10 A bestaetigt: {S1,C1}=p0 u2/2, {S1,C2}=-p0 u1/2, {C1,S2}=p0 u1/4",
    sp.simplify(pat[("(s1,c1)")] - p0 * u2 / 2) == 0
    and sp.simplify(pat[("(s1,c2)")] + p0 * u1 / 2) == 0
    and sp.simplify(pat[("(c1,s2)")] - p0 * u1 / 4) == 0,
    str(pat))

report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "sympy_version": sp.__version__,
    "python": "CPython",
    "convention": "(-+++), kappa=8 pi G_N, beta=kappa/6, h_ij=e^{2a} delta_ij, "
                  "P=2 h_ij pi^ij, phi=sqrt(beta) Phi, p=Pi_phi/beta, "
                  "pi^ij Dichte Gewicht +1 (Gewichtsterm in D_j pi_i^j)",
    "checks": CHECKS,
    "all_pass": all(v["pass"] for v in CHECKS.values()),
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)
print("\n==== SUMMARY ====")
print("checks   =", len(CHECKS))
print("all_pass =", report["all_pass"])
print("failed   =", [k for k, v in CHECKS.items() if not v["pass"]])
print("json     =", OUT)
raise SystemExit(0 if report["all_pass"] else 1)
