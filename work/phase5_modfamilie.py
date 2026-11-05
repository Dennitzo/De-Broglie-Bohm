"""Phase 5: Modalfamilie des regulierten konformen Modells und Konditionierung der
Constraintreduktion beim Entfernen der Modetrimmung.

Modell (Dichte identisch zu work/phase4_dirac_v2.py, dort gegen Direktquadratur und
gegen die geschlossene Bessel-Form validiert):
    H    = (beta/2) e^{-3 alpha} (-P^2 + p^2) - e^{alpha} (alpha')^2/(6 beta)
           + e^{alpha} (phi')^2/2 + (m^2/2) e^{3 alpha} phi^2,   Hbar = <H>
    r    = P'/3 - p phi',        C_A = <r w_A>   (Galerkin-Testfunktionen w_A)
    <f>  = (2 pi)^-1 int_0^{2 pi} f d theta   (hier: Trapetzquadratur, M Punkte)

Modalfamilie der Ordnung K (Doppelnormierung der Impulse wie in §19.7/§19.8, damit
(a_j, P_j) und (c_j, p_j) kanonische Paare sind):
    alpha = a0 + sum_k (a_k sin k th [+ b_k cos k th])
    phi   = c  + sum_k (u_k sin k th [+ v_k cos k th])
Die eckigen Klammern sind die Sektorschalter geom_cos / matter_cos. K=1 mit
geom_cos=True, matter_cos=False ist die Trimung von §19.7/§19.8 (Phasendim 10);
K=1 mit geom_cos=False ist das zugehoerige alpha_c-freie Modell (Phasendim 8).

Zweck: Zu entscheiden, welche Groessen beim Entfernen der Trimmung konvergieren
muessten, wird die Constraintreduktion (Rang der Poisson-Matrix der Constraintmenge,
Zahl der Zweitklassenpaare, Konditionierung der Dirac-Inversen, Groesse des
Dirac-Korrekturtensors) als Funktion der Modenzahl gemessen.

Rechenweise: alle Feldgroessen sind linear in den Moduskoeffizienten, die Integranden
punktweise auf dem Gitter; erste und zweite Ableitungen der Integrale exakt aus Ketten-
und Produktregel (Jet-Rechnung, keine numerische Differentiation, keine simplify-Kaskade).
Punkte sind vektoerisiert. Fortschrittsdruck mit Laufzeit, JSON-Zwischenstand je Stufe.
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np

sys.path.insert(0, "work")
import phase4_model as M  # noqa: E402

OUT = "work/phase5_modfamilie.json"
T0 = time.time()
res: dict = {}
BETA = 1.7
MASS = 0.9
MPTS = 512


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.2f}s] {msg}", flush=True)


def dump(tag: str) -> None:
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, sort_keys=True)


# ---------------------------------------------------------------------------
# Modalitaets-Apparat
# ---------------------------------------------------------------------------
def mode_list(K: int, cos: bool):
    ms = []
    for k in range(1, K + 1):
        ms.append((k, False))
        if cos:
            ms.append((k, True))
    return ms


def basis(th, modes, doubled=False):
    cols, dcols = [np.ones(th.size)], [np.zeros(th.size)]
    f = 2.0 if doubled else 1.0
    for k, is_cos in modes:
        if is_cos:
            cols.append(f * np.cos(k * th)); dcols.append(-f * k * np.sin(k * th))
        else:
            cols.append(f * np.sin(k * th)); dcols.append(f * k * np.cos(k * th))
    return np.stack(cols, axis=1), np.stack(dcols, axis=1)


class Family:
    def __init__(self, K: int, geom_cos: bool = True, matter_cos: bool = True,
                 Mpts: int = MPTS, beta: float = BETA, mass: float = MASS):
        self.K, self.geom_cos, self.matter_cos = K, geom_cos, matter_cos
        self.M, self.beta, self.mass = Mpts, beta, mass
        th = 2.0 * np.pi * np.arange(Mpts) / Mpts
        self.gm, self.mm = mode_list(K, geom_cos), mode_list(K, matter_cos)
        self.Ag, self.dAg = basis(th, self.gm)
        self.Af, self.dAf = basis(th, self.mm)
        self.APg, self.dAPg = basis(th, self.gm, doubled=True)
        self.APf, _ = basis(th, self.mm, doubled=True)
        self.ng, self.nf = self.Ag.shape[1], self.Af.shape[1]
        self.W = self.Ag                        # Galerkin-Testfunktionen = Geometriebasis
        self.nW = self.W.shape[1]
        self.nq = self.ng + self.nf
        self.ndim = 2 * self.nq
        self.nchi = 1 + self.nW                 # Hbar plus Impulsprojektionen
        self.label = f"K={K},geom_cos={int(geom_cos)},matter_cos={int(matter_cos)}"

    def kernels(self) -> np.ndarray:
        KK = np.zeros((7, self.M, self.ndim))
        KK[0, :, : self.ng] = self.Ag
        KK[1, :, : self.ng] = self.dAg
        KK[2, :, self.ng:self.nq] = self.Af
        KK[3, :, self.ng:self.nq] = self.dAf
        KK[4, :, self.nq:self.nq + self.ng] = self.APg
        KK[5, :, self.nq:self.nq + self.ng] = self.dAPg
        KK[6, :, self.nq + self.ng:] = self.APf
        return KK

    def grid_values(self, Y: np.ndarray) -> np.ndarray:
        ng, nf, nq = self.ng, self.nf, self.nq
        za, zc, Pa, Pf = Y[:, :ng], Y[:, ng:nq], Y[:, nq:nq + ng], Y[:, nq + ng:]
        out = np.empty((Y.shape[0], self.M, 7))
        out[:, :, 0] = za @ self.Ag.T
        out[:, :, 1] = za @ self.dAg.T
        out[:, :, 2] = zc @ self.Af.T
        out[:, :, 3] = zc @ self.dAf.T
        out[:, :, 4] = Pa @ self.APg.T
        out[:, :, 5] = Pa @ self.dAPg.T
        out[:, :, 6] = Pf @ self.APf.T
        return out


def gh_hbar(fam: Family, S: np.ndarray):
    a, ap, ph, php, P, _, p = (S[:, :, i] for i in range(7))
    beta, mass = fam.beta, fam.mass
    Ea, E3a, E3 = np.exp(a), np.exp(-3.0 * a), np.exp(3.0 * a)
    val = (0.5 * beta * E3a * (-P**2 + p**2) - Ea * ap**2 / (6.0 * beta)
           + 0.5 * Ea * php**2 + 0.5 * mass**2 * E3 * ph**2)
    g1 = np.zeros(S.shape[:2] + (7,))
    g1[:, :, 0] = (-1.5 * beta * E3a * (-P**2 + p**2) - Ea * ap**2 / (6.0 * beta)
                   + 0.5 * Ea * php**2 + 1.5 * mass**2 * E3 * ph**2)
    g1[:, :, 1] = -Ea * ap / (3.0 * beta)
    g1[:, :, 2] = mass**2 * E3 * ph
    g1[:, :, 3] = Ea * php
    g1[:, :, 4] = -beta * E3a * P
    g1[:, :, 6] = beta * E3a * p
    g2 = np.zeros(S.shape[:2] + (7, 7))
    # dg/da = -1.5 beta E3a (-P^2+p^2) - Ea ap^2/(6 beta) + Ea php^2/2 + 1.5 m^2 E3 phi^2
    # d/da davon: ( -1.5 * -3 ) = 4.5 und ( 1.5 * 3 ) = 4.5
    g2[:, :, 0, 0] = (4.5 * beta * E3a * (-P**2 + p**2) - Ea * ap**2 / (6.0 * beta)
                      + 0.5 * Ea * php**2 + 4.5 * mass**2 * E3 * ph**2)
    g2[:, :, 0, 1] = g2[:, :, 1, 0] = -Ea * ap / (3.0 * beta)
    g2[:, :, 0, 2] = g2[:, :, 2, 0] = 3.0 * mass**2 * E3 * ph
    g2[:, :, 0, 3] = g2[:, :, 3, 0] = Ea * php
    g2[:, :, 0, 4] = g2[:, :, 4, 0] = 3.0 * beta * E3a * P
    g2[:, :, 0, 6] = g2[:, :, 6, 0] = -3.0 * beta * E3a * p
    g2[:, :, 1, 1] = -Ea / (3.0 * beta)
    g2[:, :, 2, 2] = mass**2 * E3
    g2[:, :, 3, 3] = Ea
    g2[:, :, 4, 4] = -beta * E3a
    g2[:, :, 6, 6] = beta * E3a
    return val, g1, g2


def gh_momentum(fam: Family, S: np.ndarray, A: int):
    w = fam.W[:, A]
    php, Pp, p = S[:, :, 3], S[:, :, 5], S[:, :, 6]
    val = w[None, :] * (Pp / 3.0 - p * php)
    g1 = np.zeros(S.shape[:2] + (7,))
    g1[:, :, 3] = -w[None, :] * p
    g1[:, :, 5] = w[None, :] / 3.0
    g1[:, :, 6] = -w[None, :] * php
    g2 = np.zeros(S.shape[:2] + (7, 7))
    g2[:, :, 3, 6] = g2[:, :, 6, 3] = -w[None, :]
    return val, g1, g2


def jets(fam: Family, Y: np.ndarray, kind: str, A: int = 0, want_hess: bool = False):
    S = fam.grid_values(Y)
    val_g, g1, g2 = (gh_hbar(fam, S) if kind == "H" else gh_momentum(fam, S, A))
    KK = fam.kernels()
    Minv = 1.0 / fam.M
    val = Minv * val_g.sum(axis=1)
    grad = Minv * np.einsum("nms,smu->nu", g1, KK)
    if not want_hess:
        return val, grad, None
    hess = Minv * np.einsum("nist,siu,tiv->nuv", g2, KK, KK, optimize=True)
    return val, grad, hess


def pbf(fam: Family, gF: np.ndarray, gG: np.ndarray) -> np.ndarray:
    """Poisson-Klammer aus Gradienten (..., ndim)."""
    nq = fam.nq
    return (np.einsum("...u,...u->...", gF[..., :nq], gG[..., nq:])
            - np.einsum("...u,...u->...", gF[..., nq:], gG[..., :nq]))


def grad_bracket(fam: Family, hessF, gF, hessG, gG) -> np.ndarray:
    """Gradient der Klammer {F,G} nach allen Variablen (Produktregel, beide Hessians)."""
    nq = fam.nq
    t1 = np.einsum("...uv,...u->...v", hessF[..., :nq, :], gG[..., nq:])
    t2 = np.einsum("...u,...uv->...v", gF[..., :nq], hessG[..., nq:, :])
    t3 = np.einsum("...uv,...u->...v", hessF[..., nq:, :], gG[..., :nq])
    t4 = np.einsum("...u,...uv->...v", gF[..., nq:], hessG[..., :nq, :])
    return t1 + t2 - t3 - t4


def chi_all(fam: Family, Y: np.ndarray):
    h, gH = jets(fam, Y, "H")[0:2]
    vals, rows = [h[:, None]], [gH[:, None, :]]
    for A in range(fam.nW):
        c, gC = jets(fam, Y, "C", A)[0:2]
        vals.append(c[:, None])
        rows.append(gC[:, None, :])
    return np.concatenate(vals, axis=1), np.concatenate(rows, axis=1)


def newton_surface(fam: Family, Y0: np.ndarray, iters: int = 80, tol: float = 1e-13):
    """Gedaempfter Gauss-Newton auf chi = 0; rueckgabe (Y, konvergenzmaske, rest_max, used)."""
    Y = Y0.copy()
    n = Y.shape[0]
    used = np.zeros(n, dtype=int)
    eye = np.eye(fam.nchi)[None]
    for _ in range(iters):
        chi, Jac = chi_all(fam, Y)
        nrm = np.linalg.norm(chi, axis=1)
        if np.all(nrm < tol):
            break
        used += (nrm >= tol)
        JJT = np.einsum("naj,nbj->nab", Jac, Jac)
        scale = np.max(np.abs(JJT.reshape(n, -1)), axis=1)[:, None, None]
        lam = np.linalg.solve(JJT + (1e-12 * scale + 1e-300) * eye, chi[..., None])[..., 0]
        step = np.einsum("naj,na->nj", Jac, lam)
        for theta in (1.0, 0.5, 0.25, 0.1, 0.02):
            Yn = Y - theta * step
            chi_n, _ = chi_all(fam, Yn)
            better = np.linalg.norm(chi_n, axis=1) < nrm
            Y[better] = Yn[better]
            if bool(better.all()):
                break
    chi, _ = chi_all(fam, Y)
    nrm = np.linalg.norm(chi, axis=1)
    return Y, nrm < 1e-9, float(nrm.max()), used


def analyse_omega(stack: np.ndarray, tol_rel: float = 1e-9) -> dict:
    """Rang- und Konditionsstatistik ueber einen Stapel (N, nchi, nchi) antisymmetrischer Matrices."""
    if stack.ndim == 2:
        stack = stack[None]
    ranks, smaxs, smins, kond = [], [], [], []
    for om in stack:
        sv = np.linalg.svd(om, compute_uv=False)
        smax = float(sv[0]) if sv.size else 0.0
        nz = sv[sv > tol_rel * max(1.0, smax)]
        ranks.append(nz.size)
        smaxs.append(smax)
        smins.append(float(nz[-1]) if nz.size else 0.0)
        kond.append(smax / nz[-1] if nz.size else float("inf"))
    return {
        "rang_med": float(np.median(ranks)),
        "rang_verteilung": {str(int(r)): ranks.count(r) for r in sorted(set(ranks))},
        "sigma_max_med": float(np.median(smaxs)),
        "sigma_min_nz_med": float(np.median(smins)),
        "sigma_min_nz_min": float(np.min(smins)),
        "kondition_med": float(np.median(kond)),
        "kondition_max": float(np.max(kond)),
        "omega_max": float(np.max(np.abs(stack))),
    }


def klammer_covector(fam: Family, Jac: np.ndarray) -> np.ndarray:
    """Klammerfeld cb[n,a,u] = {z_u, chi_a} aus der Constraint-Jacobi-Matrix.

    Jac[n,a,u] = d_{z_u} chi_a mit Variablenblock (q_1..q_nq, p_1..p_nq). Fuer
    u = q_j bleibt {q_j, chi_a} = d_{p_j} chi_a, fuer u = p_j wird es
    -d_{q_j} chi_a; beide Bloecke werden deshalb entlang der VARIABLENachse
    (Achse -1) umgespeist und vertauscht. Ein Schneiden entlang der Achse 1 (wie im
    ersten Entwurf) haette die Constraintachse verdoppelt und die nachfolgende
    Matrixmultiplikation in der Achsenzahl kollidieren lassen.
    """
    return np.concatenate([Jac[..., fam.nq:], -Jac[..., :fam.nq]], axis=-1)   # (N,nchi,ndim)


def dirac_korrektur(fam: Family, JC: np.ndarray, omega: np.ndarray) -> dict:
    """Frobeniusnorm des Dirac-Korrekturtensors C_jk = {z_j,chi_a}(omega^+)^{ab}{chi_b,z_k}.

    JC ist die Jacobi-Matrix GENAU der Menge, die reduziert wird (Form (N,m,ndim)),
    omega die zugehoerige Klammermatrix (N,m,m). Die Beschraenkung ist noetig, weil
    Hbar nie zur Dirac-Menge gehoeren darf: sonst ist {F,Hbar}_D identisch Null und
    jede gemessene "Korrektur" ist die Tautologie des Verfahrens, nicht die Dynamik.
    """
    cb = klammer_covector(fam, JC)
    pinv = omega_pinv(omega)
    C = np.einsum("naj,nab,nbk->njk", cb, pinv, cb, optimize=True)
    a = np.linalg.norm(C, axis=(1, 2))
    return {"korrektur_norm_max": float(a.max()), "korrektur_norm_med": float(np.median(a))}


def sample(fam: Family, n: int, rng, amp: float = 0.2):
    """Zufaellige Punkte mit geddaempftem Modespektrum (Amplitude amp/k, glatte Daten).

    Die Daempfung ist keine mathematische Voraussetzung des Modells, sondern haelt die
    Randdaten ueber alle Modenordnungen gleichmaessig glatt. Sonst waechse mit K der
    Hochfrequenzanteil der Stichprobe und die gemessene Konditionierung der
    Klammermatrix waere eine Eigenschaft der Stichprobe und nicht des Modells.
    """
    def fill(blk: np.ndarray, mlist) -> None:
        for j in range(1, blk.shape[1]):
            k = mlist[j - 1][0]
            blk[:, j] = rng.normal(0.0, amp / k, size=n)
    za = np.zeros((n, fam.ng)); za[:, 0] = 0.3 + 0.1 * rng.normal(size=n)
    zc = np.zeros((n, fam.nf)); zc[:, 0] = 0.6 + 0.1 * rng.normal(size=n)
    Pa = np.zeros((n, fam.ng)); Pa[:, 0] = -1.3 + 0.3 * rng.normal(size=n)
    Pf = np.zeros((n, fam.nf)); Pf[:, 0] = 0.8 + 0.2 * rng.normal(size=n)
    fill(za, fam.gm); fill(Pa, fam.gm)
    fill(zc, fam.mm); fill(Pf, fam.mm)
    return np.concatenate([za, zc, Pa, Pf], axis=1)


# ---------------------------------------------------------------------------
# Stufen
# ---------------------------------------------------------------------------
def referenzpunkte():
    base = {"a0": 0.3, "a1": 0.25, "c": 0.6, "q": 0.4, "P1": 1.5 * 0.8 * 0.4,
            "Pc": 0.0, "p0": 0.8, "pq": -0.35}
    pts = {}
    for tag, P0v in (("wurzel_1", -1.34014194788918786619041671333),
                     ("wurzel_2", 2.01384482125922654421633828245)):
        pts[tag] = [base["a0"], base["a1"], 0.0, base["c"], base["q"],
                    P0v, base["P1"], base["Pc"], base["p0"], base["pq"]]
    return pts


def stage_V() -> None:
    log("Stufe V: Quadraturapparat gegen geschlossene Bessel-Form (§19.7/§19.8)")
    Hbar = M.hbar_closed()
    Cs, Cc = M.constraints()
    Dsym = M.pb_sym(Cc, Hbar)
    Wsym = M.pb_sym(Cc, Dsym)
    layer = M.FastLayer({"Hbar": Hbar, "Cs": Cs, "Cc": Cc, "D": Dsym, "Dcd": Wsym},
                        {M.beta: BETA, M.m: MASS})
    y = np.array([referenzpunkte()["wurzel_1"]])
    closed = {n: float(layer.eval_np(n, y)[0][0]) for n in ("Hbar", "Cs", "Cc", "D", "Dcd")}
    out = {}
    for tag, fam in (("K1_trim_alpha_c_vollstaendig", Family(1, True, False)),):
        pg = (fam.APg.T @ fam.Ag) / fam.M
        pf = (fam.APf.T @ fam.Af) / fam.M
        pair = max(float(np.max(np.abs(pg - np.eye(fam.ng)))),
                   float(np.max(np.abs(pf - np.eye(fam.nf)))))
        yy = y.copy()
        hV, hG, hH = jets(fam, yy, "H", want_hess=True)
        jetsC = {A: jets(fam, yy, "C", A, want_hess=True) for A in range(fam.nW)}
        iCc = fam.nW - 1                      # Cc = <r cos> bzw. letzte Projektion
        cV, cG, cH = jetsC[iCc]
        Dq = pbf(fam, cG, hG)[0]
        gD = grad_bracket(fam, cH, cG, hH, hG)
        Wq = pbf(fam, cG, gD)[0]
        row = {
            "phasendim": fam.ndim, "freiheitsgrade": fam.nq,
            "paarung_max_abweichung": pair,
            "Hbar": float(hV[0]), "Cs": float(jetsC[1][0][0]), "Cc": float(cV[0]),
            "D": float(Dq), "W": float(Wq),
            "diff": {"Hbar": abs(closed["Hbar"] - float(hV[0])),
                     "Cs": abs(closed["Cs"] - float(jetsC[1][0][0])),
                     "Cc": abs(closed["Cc"] - float(cV[0])),
                     "D": abs(closed["D"] - float(Dq)),
                     "W": abs(closed["Dcd"] - float(Wq))},
        }
        out[tag] = row
        log(f"   {tag}: Phasendim {fam.ndim}, Paarung {pair:.2e}, "
            f"D Diff {row['diff']['D']:.2e}, W Diff {row['diff']['W']:.2e} "
            f"(W quad {row['W']:.12f}, W geschlossen {closed['Dcd']:.12f})")
    out["closed_form"] = closed
    res["V"] = out
    dump("V")


def stage_V2() -> None:
    """Unabhaengige Validierung: D und W aus reiner mpmath-Quadratur mit Differentenzahlen.

    Referenz ist allein das Integral <H> (keine Ableitungen, keine Bessel-Identitaeten),
    alle Klammerableitungen entstehen aus zentralen Differenzen desselben Integrals.
    """
    from mpmath import mp, mpf, quad, sin as ms, cos as mc, exp as me

    mp.dps = 32
    log("Stufe V2: D und W aus mpmath-Quadratur des Integrals (Differentenzahlen)")
    fam = Family(1, geom_cos=True, matter_cos=False)
    gm, mm, ng, nf, nq = fam.gm, fam.mm, fam.ng, fam.nf, fam.nq
    nd = fam.ndim
    beta_mp, mass_mp = mpf(str(BETA)), mpf(str(MASS))
    two = mpf(2)

    def Hbar_mp(z):
        za, zc = z[:ng], z[ng:nq]
        Pa, Pf = z[nq:nq + ng], z[nq + ng:]

        def dens(t):
            a = mpf(za[0]); ap = mpf(0)
            for j, (k, isc) in enumerate(gm):
                s = mc(k * t) if isc else ms(k * t)
                d = -mpf(k) * ms(k * t) if isc else mpf(k) * mc(k * t)
                a += mpf(za[j + 1]) * s
                ap += mpf(za[j + 1]) * d
            ph = mpf(zc[0]); php = mpf(0)
            for j, (k, isc) in enumerate(mm):
                s = mc(k * t) if isc else ms(k * t)
                d = -mpf(k) * ms(k * t) if isc else mpf(k) * mc(k * t)
                ph += mpf(zc[j + 1]) * s
                php += mpf(zc[j + 1]) * d
            P = mpf(Pa[0]); Pp = mpf(0)
            for j, (k, isc) in enumerate(gm):
                s = mc(k * t) if isc else ms(k * t)
                d = -mpf(k) * ms(k * t) if isc else mpf(k) * mc(k * t)
                P += two * mpf(Pa[j + 1]) * s
                Pp += two * mpf(Pa[j + 1]) * d
            p = mpf(Pf[0])
            for j, (k, isc) in enumerate(mm):
                p += two * mpf(Pf[j + 1]) * (mc(k * t) if isc else ms(k * t))
            return (beta_mp / 2) * me(-3 * a) * (-P**2 + p**2) - me(a) * ap**2 / (6 * beta_mp) \
                + me(a) * php**2 / 2 + (mass_mp**2 / 2) * me(3 * a) * ph**2

        return quad(dens, [0, 2 * mp.pi], error=False) / (2 * mp.pi)

    h = mpf("1e-5")

    def gradH(z):
        g = []
        for u in range(nd):
            zp = list(z); zm = list(z)
            zp[u] += h; zm[u] -= h
            g.append((Hbar_mp(zp) - Hbar_mp(zm)) / (2 * h))
        return g

    def bracket_mp(dF, gG):
        s = mpf(0)
        for u in range(nq):
            s += dF[u] * gG[nq + u] - dF[nq + u] * gG[u]
        return s

    def dCc_vec(z):
        # Cc = P1/3 - p0 q/2 in Indizes: q = ng+1, P1 = nq+1, p0 = nq+ng
        d = [mpf(0)] * nd
        d[ng + 1] = -mpf(z[nq + ng]) / 2
        d[nq + 1] = mpf(1) / 3
        d[nq + ng] = -mpf(z[ng + 1]) / 2
        return d

    def D_mp(z):
        return bracket_mp(dCc_vec(z), gradH(z))

    def W_mp(z):
        gD = []
        for u in range(nd):
            zp, zm = list(z), list(z)
            zp[u] += h; zm[u] -= h
            gD.append((D_mp(zp) - D_mp(zm)) / (2 * h))
        return bracket_mp(dCc_vec(z), gD)

    y = referenzpunkte()["wurzel_1"]
    ym = [mpf(repr(v)) for v in y]
    Dmp = D_mp(ym)
    Wmp = W_mp(ym)
    hV, hG, hH = jets(fam, np.array([y], dtype=float), "H", want_hess=True)
    jetsC = {A: jets(fam, np.array([y], dtype=float), "C", A, want_hess=True)
             for A in range(fam.nW)}
    cV, cG, cH = jetsC[fam.nW - 1]
    Dq = float(pbf(fam, cG, hG)[0])
    gD = grad_bracket(fam, cH, cG, hH, hG)
    Wq = float(pbf(fam, cG, gD)[0])
    # alpha_c-Richtungsableitungen von D (zeigen, dass die Beschraenkung nicht zufaellig ist)
    gD_np = np.asarray(gD).ravel()
    out = {
        "punkt": "wurzel_1",
        "D_mp": float(Dmp), "D_jet": Dq, "D_diff": abs(float(Dmp) - Dq),
        "W_mp": float(Wmp), "W_jet": Wq, "W_diff": abs(float(Wmp) - Wq),
        "dD_dalphac": float(gD_np[2]), "dD_dPc": float(gD_np[fam.nq + 2]),
    }
    res["V2"] = out
    log(f"   D: mp {out['D_mp']:.12f}, Jet {out['D_jet']:.12f}, Diff {out['D_diff']:.2e}")
    log(f"   W: mp {out['W_mp']:.12f}, Jet {out['W_jet']:.12f}, Diff {out['W_diff']:.2e}")
    log(f"   Richtungsableitungen von D: nach alpha_c {out['dD_dalphac']:.4e}, "
        f"nach P_c {out['dD_dPc']:.4e}")
    dump("V2")


def stage_V3() -> None:
    """Diagnose: Jet-Hessian gegen zentrale Differenzen des validierten Gradienten."""
    fam = Family(1, geom_cos=True, matter_cos=False)
    y = np.array([referenzpunkte()["wurzel_1"]], dtype=float)
    _, _, hH = jets(fam, y, "H", want_hess=True)
    _, _, hC = jets(fam, y, "C", fam.nW - 1, want_hess=True)
    h = 1e-6
    g_plus = jets(fam, y + h * np.eye(fam.ndim)[0:1].T if False else _bump(y, h, fam), "H")[1]
    FDH = np.zeros((fam.ndim, fam.ndim))
    for u in range(fam.ndim):
        gp = jets(fam, _bump(y, h, fam, u), "H")[1][0]
        gm = jets(fam, _bump(y, -h, fam, u), "H")[1][0]
        FDH[:, u] = (gp - gm) / (2 * h)
    FDC = np.zeros((fam.ndim, fam.ndim))
    for u in range(fam.ndim):
        gp = jets(fam, _bump(y, h, fam, u), "C", fam.nW - 1)[1][0]
        gm = jets(fam, _bump(y, -h, fam, u), "C", fam.nW - 1)[1][0]
        FDC[:, u] = (gp - gm) / (2 * h)
    DH = np.abs(hH[0] - FDH)
    DC = np.abs(hC[0] - FDC)
    idx = np.unravel_index(np.argmax(DH), DH.shape)
    res["V3"] = {
        "hessH_max_abweichung": float(DH.max()),
        "hessH_max_stelle": [int(idx[0]), int(idx[1])],
        "hessH_strahl": [float(hH[0].max()), float(FDH.max())],
        "hessC_max_abweichung": float(DC.max()),
    }
    log(f"   Hessian Hbar: max Abweichung {DH.max():.3e} an Stelle {idx}, "
        f"Betragsspanne Jet {hH[0].max():.3e} / FD {FDH.max():.3e}")
    log(f"   Hessian Cc  : max Abweichung {DC.max():.3e}")
    dump("V3")


def _bump(y: np.ndarray, h: float, fam: Family, u: int = -1) -> np.ndarray:
    out = y.copy()
    if u >= 0:
        out[0, u] += h
    return out


def omega_pinv(omega: np.ndarray, tol_rel: float = 1e-9) -> np.ndarray:
    """Pseudoinverse eines Stapels antisymmetrischer Matrices (SVD, vektoerisiert)."""
    U, s, Vt = np.linalg.svd(omega)
    keep = s > tol_rel * np.maximum(1.0, s[:, :1])
    sinv = np.where(keep, 1.0 / np.where(s > 0.0, s, 1.0), 0.0)
    return (Vt.transpose(0, 2, 1) * sinv[:, None, :]) @ U.transpose(0, 2, 1)


def pbf_vec(gA: np.ndarray, gB: np.ndarray, nq: int) -> np.ndarray:
    """Klammer jedes Eintrags von gA (N,m,ndim) mit der Funktion gB (N,ndim)."""
    return (np.einsum("nau,nu->na", gA[..., :nq], gB[..., nq:])
            - np.einsum("nau,nu->na", gA[..., nq:], gB[..., :nq]))


def flow_masse(fam: Family, Y: np.ndarray, JC: np.ndarray, omega: np.ndarray) -> dict:
    """Naives Feld {z,Hbar} gegen das Dirac-Feld auf der Constraintflaeche.

    Das Dirac-Feld ist das naive Feld abzgl. {z_j,chi_a}(omega^+)^{ab}{chi_b,Hbar}.
    Wächst der relative Korrekturanteil mit der Modenzahl K, wird die Reduktion im
    nicht getrimmten Limes NICHT gleichmaessig klein; die Constraintloesung haengt
    dann staerker an der Inversen der Klammermatrix als an der Dynamik selbst. Das ist
    die messbare Konditionsaussage zur Constraintreduktion (Nullwert heisst: Hbar
    commute schwach mit allen Constraints, Reduktion korrigiert die Fuehrung nicht).
    """
    nq = fam.nq
    _, gH, _ = jets(fam, Y, "H")
    naive = np.concatenate([gH[..., nq:], -gH[..., :nq]], axis=1)      # (N, ndim)
    cb = klammer_covector(fam, JC)                                      # (N, m, ndim)
    chiH = pbf_vec(JC, gH, nq)                                          # (N, m)
    pinv = omega_pinv(omega)
    corr = np.einsum("naj,nab,nb->nj", cb, pinv, chiH)
    nn = np.linalg.norm(naive, axis=1)
    cn = np.linalg.norm(corr, axis=1)
    rel = cn / np.maximum(nn, 1e-300)
    return {
        "feld_naive_med": float(np.median(nn)), "feld_naive_max": float(nn.max()),
        "feld_korrektur_rel_med": float(np.median(rel)),
        "feld_korrektur_rel_max": float(rel.max()),
        "klammer_chi_Hbar_max": float(np.max(np.abs(chiH))),
        "omega_inv_norm_max": float(np.max(np.linalg.norm(pinv, axis=(1, 2)))),
    }


def stage_F(variants, npts: int = 48, seed: int = 20261105) -> None:
    log(f"Stufe F: Modalfamilie, {len(variants)} Varianten, {npts} Randdaten je Variante")
    rng = np.random.default_rng(seed)
    rows: list = []
    for K, gcos, mcos in variants:
        tf = time.time()
        fam = Family(K, geom_cos=gcos, matter_cos=mcos)
        Y0 = sample(fam, npts, rng)
        Y, conv, restmax, used = newton_surface(fam, Y0)
        Ys = Y[conv]
        row = {
            "variante": fam.label, "K": K, "geom_cos": int(gcos), "matter_cos": int(mcos),
            "phasendim": fam.ndim, "constraints": fam.nchi, "randdaten": npts,
            "newton_konvergiert": int(conv.sum()), "newton_rest_max": restmax,
            "newton_iterationen_med": float(np.median(used)),
        }
        if Ys.shape[0] == 0:
            row.update({"rang": float("nan"), "bemerkung": "keine Constraintloesung"})
            rows.append(row)
            res["F"] = rows
            dump("F")
            log(f"   {fam.label}: Newton ohne Loesung (Rest {restmax:.1e})")
            continue
        chi, Jac = chi_all(fam, Ys)
        # Reduziert wird ausschliesslich die Constraintschar chi = C_A (Zeile 0 ist
        # Hbar und bleibt Generator). Die volle Menge (Hbar, C_A) wird nur als
        # Referenz ausgewiesen, weil ihre Klammermatrix wegen Hbar als Zeile/Spalte
        # nie die fuer die Reduktion massgebliche Inverse liefert.
        JC = Jac[:, 1:, :]
        omega = (np.einsum("nAu,nBu->nAB", JC[..., :fam.nq], JC[..., fam.nq:])
                 - np.einsum("nAu,nBu->nAB", JC[..., fam.nq:], JC[..., :fam.nq]))
        an = analyse_omega(omega)
        omega_full = (np.einsum("nAu,nBu->nAB", Jac[..., :fam.nq], Jac[..., fam.nq:])
                      - np.einsum("nAu,nBu->nAB", Jac[..., fam.nq:], Jac[..., :fam.nq]))
        an_full = analyse_omega(omega_full)
        corr = dirac_korrektur(fam, JC, omega)
        fm = flow_masse(fam, Ys, JC, omega)
        row.update({
            "rang": an["rang_med"], "rang_defizit": fam.nW - an["rang_med"],
            "rang_verteilung": an["rang_verteilung"],
            "rang_gesetztheit": fam.nW,
            "zweitklassenpaare": int(an["rang_med"] // 2),
            "erstklassige": fam.nW - int(an["rang_med"]),
            "sigma_max_med": an["sigma_max_med"], "sigma_min_nz_med": an["sigma_min_nz_med"],
            "sigma_min_nz_min": an["sigma_min_nz_min"],
            "kondition_med": an["kondition_med"], "kondition_max": an["kondition_max"],
            "omega_max": an["omega_max"],
            "rang_gesamt": an_full["rang_med"], "kondition_gesamt_med": an_full["kondition_med"],
            "dirac_korrektur_max": corr["korrektur_norm_max"],
            "dirac_korrektur_med": corr["korrektur_norm_med"],
            "constraint_rest_max": float(np.max(np.abs(chi))),
            "zeit_s": time.time() - tf,
        })
        row.update(fm)
        rows.append(row)
        log(f"   {fam.label}: Phasendim {fam.ndim}, Constraints {fam.nchi}, "
            f"Newton {row['newton_konvergiert']}/{npts} (Rest {restmax:.1e}), "
            f"Rang {row['rang']:.0f}, Zweitklassenpaare {row['zweitklassenpaare']}, "
            f"sigma_min {row['sigma_min_nz_med']:.3e}, Kondition {row['kondition_med']:.3e}, "
            f"Dirac-Korrektur {row['dirac_korrektur_med']:.3e}, "
            f"relative Feldkorrektur {row['feld_korrektur_rel_med']:.3e} "
            f"({row['zeit_s']:.1f}s)")
        res["F"] = rows
        dump("F")
    res["F_meta"] = {"npts": npts, "seed": seed}
    dump("F")


def load_existing() -> None:
    """Vorhandene Stufenergebnisse laden, damit einzelne Etappen getrennt lauffaehig sind."""
    global res
    try:
        with open(OUT, encoding="utf-8") as f:
            res = json.load(f)
        log(f"geladener Zwischenstand: {sorted(res)}")
    except (OSError, ValueError):
        res = {}


ALL_STAGES = ("V", "V2", "V3", "F")
VARIANTS = [(1, True, False), (1, True, True), (2, True, True), (3, True, True),
            (4, True, True), (2, True, False), (3, True, False)]


def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    picked = {a.upper() for a in argv if not a.startswith("-")}
    wanted = [s for s in ALL_STAGES if not picked or s.upper() in picked]
    load_existing()
    log(f"Ausgefuehrte Stufen: {wanted}")
    if "V" in wanted:
        stage_V()
    if "V2" in wanted:
        stage_V2()
    if "V3" in wanted:
        stage_V3()
    if "F" in wanted:
        stage_F(VARIANTS)
    log(f"fertig: {wanted}")


if __name__ == "__main__":
    main()
