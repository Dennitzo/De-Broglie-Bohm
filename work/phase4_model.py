"""Gemeinsames Modell und schneller numerischer Apparat fuer Phase 4 (S19.7/S19.8).

Modell (identisch zu work/phase4_dirac_v2.py, dort gegen Direktquadratur validated):
  Konfiguration  alpha = a0 + a1 sin(theta) + ac cos(theta),  phi = c + q sin(theta)
  Impulse        P = P0 + 2 P1 sin(theta) + 2 Pc cos(theta),  p = p0 + 2 pq sin(theta)
  kanonische Paare (a0,P0), (a1,P1), (ac,Pc), (c,p0), (q,pq)
  Hbar = <H> in geschlossener Bessel-Form auf der Scheibe ac = 0
  Constraints aus r = P'/3 - p phi':  Cs = <r sin> = -Pc/3,  Cc = <r cos> = P1/3 - p0 q/2

Zweck des Apparats: die fruehere Fassung hat Klammern und Felder vollsymbolisch mit
simplify gebaut und anschiessend mit mpmath seriell integriert (Stundenlauf ohne
Ausgabe). Hier wird deshalb getrennt:

  * symbolisch nur einmal ableiten und expandieren (keine simplify-Kaskade),
  * Bessel-Funktionen nach DER Differentiation zu Hilfs Symbolen umschreiben - die
    Kettenregel nach alpha_1 bleibt exakt erhalten, weil die Ableitung vorher erfolgt,
  * Werte und Gradienten aller Constraint-Funktionen einmal lambdifizieren, danach
    laeuft alles als reine Arithmetik auf numpy-Feldern (vektoerierte Trajektorien),
  * derselbe Ausdruck additionally als mpmath-Pfad ohne Hilfs Symbole, damit die
    schnelle Route unabhaengig gegen die genaue Route geprueft werden kann.

Konvention der Poisson-Klammer: {f,g} = sum_i (d_{q_i} f d_{p_i} g - d_{p_i} f d_{q_i} g).
"""

from __future__ import annotations

import mpmath
import numpy as np
import sympy as sp
from scipy.special import iv

# Symbolnamen wie in work/phase4_dirac_v2.py, damit gespeicherte Ausdruecke exakt
# verglichen werden koennen.
alpha0, alpha1, alpha_c, c, q = sp.symbols("alpha0 alpha1 alpha_c c q", real=True)
P0, P1, Pc, p0, pq = sp.symbols("P0 P1 Pc p0 pq", real=True)
beta, m = sp.symbols("beta m", positive=True)
QS = (alpha0, alpha1, alpha_c, c, q)
PS = (P0, P1, Pc, p0, pq)
Z = QS + PS
NCAN = len(Z)


def bessel_basis(z):
    """Besselfunktionen I_n(z) fuer n = 0..4 als Liste."""
    return [sp.besseli(n, z) for n in range(5)]


def hbar_closed() -> sp.Expr:
    """Geschlossene Form von Hbar auf der Scheibe alpha_c = 0 (Doppelnormierung)."""
    i = lambda n: sp.besseli(n, 3 * alpha1)
    j = lambda n: sp.besseli(n, alpha1)
    cz1 = (j(0) - j(2)) / 2
    geo = -(
        P0**2 * i(0)
        - 4 * P0 * P1 * i(1)
        + 2 * P1**2 * (i(0) + i(2))
        + 2 * Pc**2 * (i(0) - i(2))
    )
    matter = p0**2 * i(0) - 4 * p0 * pq * i(1) + 2 * pq**2 * (i(0) + i(2))
    pot = (
        -sp.exp(alpha0) * alpha1**2 * cz1 / (6 * beta)
        + sp.exp(alpha0) * q**2 * cz1 / 2
        + m**2 * sp.exp(3 * alpha0) * (c**2 * i(0) + 2 * c * q * i(1) + q**2 * (i(0) + i(2)) / 2) / 2
    )
    return sp.expand(beta * sp.exp(-3 * alpha0) * (geo + matter) / 2 + pot)


def constraints() -> tuple[sp.Expr, sp.Expr]:
    return sp.Rational(-1, 3) * Pc, P1 / 3 - p0 * q / 2


def pb_sym(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    """Poisson-Klammer symbolisch, nur expandieren (kein simplify: Bessel-Ausdruecke)."""
    acc = 0
    for fq, fp in zip(QS, PS):
        acc = acc + sp.diff(f, fq) * sp.diff(g, fp) - sp.diff(f, fp) * sp.diff(g, fq)
    return sp.expand(acc)


def monomials(expr: sp.Expr, keys=(P0, P1, Pc, p0, pq)):
    """Koeffiziententabelle ueber den Impulsen, domain='EX' (exakte Koeffizienten)."""
    poly = sp.Poly(sp.expand(expr), *keys, domain="EX")
    return [[str(mm), str(cf)] for mm, cf in zip(poly.monoms(), poly.coeffs())]


# ---------------------------------------------------------------------------
# schneller numerischer Apparat
# ---------------------------------------------------------------------------
def _collect(expr: sp.Expr, pairs: set) -> None:
    """Sammelt alle Bessel-Termen (Ordnung, Vorfaktor von alpha_1 im Argument)."""
    for e in sp.preorder_traversal(expr):
        if isinstance(e, sp.besseli):
            order = int(e.args[0])
            arg = sp.sympify(e.args[1])
            coeff = arg.coeff(alpha1)
            if not coeff or sp.simplify(arg - coeff * alpha1) != 0:
                raise ValueError(f"unerwartetes Bessel-Argument: {arg}")
            pairs.add((order, sp.nsimplify(coeff)))


class FastLayer:
    """Wert und Gradient jeder Funktion, zweimal: numpy (schnell) und mpmath (genau).

    Die numpy-Route enthaelt keine Bessel-Aufrufe mehr: alle Bessel-Termen sind zu
    Hilfs Symbolen geworden und werden aus scipy.special.iv einmal pro Auswertung
    eingesetzt. Die mpmath-Route uses die Originalausdruecke, ist also unabhaengig.
    """

    def __init__(self, funcs: dict[str, sp.Expr], const_values: dict | None = None):
        self.funcs = {k: sp.expand(v) for k, v in funcs.items()}
        # Modellparameter als feste Werte: beta und m sind im Modell Konstanten und
        #muessen in den lambdifizierten Argumentlisten stehen, sonst fehlen sie im
        #generierten Code (NameError). Beide Pfade erhalten dieselben Werte.
        items = sorted((const_values or {}).items(), key=lambda kv: str(kv[0]))
        self.const_syms = [k for k, _ in items]
        self.const_vals = [float(v) for _, v in items]
        self.const_mp = [mpmath.mpf(repr(v)) for _, v in items]
        pairs: set = set()
        self.grad_sym: dict[str, list] = {}
        for name, f in self.funcs.items():
            _collect(f, pairs)
            self.grad_sym[name] = [sp.diff(f, z) for z in Z]
            for g in self.grad_sym[name]:
                _collect(g, pairs)
        self.pairs = sorted(pairs, key=lambda t: (float(t[1]), t[0]))
        self.syms = [sp.Symbol(f"B{order}_{int(coeff * 100)}") for order, coeff in self.pairs]
        repl = {sp.besseli(order, coeff * alpha1): s for (order, coeff), s in zip(self.pairs, self.syms)}
        self.repl = repl
        self.args = list(Z) + self.const_syms + self.syms

        def fast(expr: sp.Expr) -> sp.Expr:
            out = expr.xreplace(repl)
            if any(isinstance(t, sp.besseli) for t in sp.preorder_traversal(out)):
                raise ValueError("Bessel-Rest nach Umschreibung der Hilfs Symbole")
            return out

        self._np = {}
        self._mp = {}
        for name, f in self.funcs.items():
            outs_fast = [fast(f)] + [fast(g) for g in self.grad_sym[name]]
            self._np[name] = sp.lambdify(self.args, outs_fast, "numpy")
            outs_mp = [f] + list(self.grad_sym[name])
            self._mp[name] = sp.lambdify(list(Z) + self.const_syms, outs_mp, "mpmath")

    # --- numpy-Route -------------------------------------------------------
    def aux_values(self, a1_values: np.ndarray) -> list[np.ndarray]:
        a1 = np.asarray(a1_values, dtype=float)
        return [iv(order, float(coeff) * a1) for order, coeff in self.pairs]

    def eval_np(self, name: str, Y: np.ndarray):
        """Wert (N,) und Gradient (N,10) einer Funktion an N Punkten."""
        Y = np.atleast_2d(np.asarray(Y, dtype=float))
        aux = self.aux_values(Y[:, 1])
        cols = [Y[:, i] for i in range(NCAN)]
        outs = self._np[name](*cols, *self.const_vals, *aux)
        n = Y.shape[0]
        # Ableitungen nach Variablen, die in einer Funktion nicht vorkommen, liefert
        # lambdify als skalare Null. Ohne Broadcast bekamen sie Form (1,) und passten
        # nicht zu den uebrigen Komponenten der Form (N,).
        val = np.broadcast_to(np.asarray(outs[0], dtype=float), (n,))
        grad = np.stack([np.broadcast_to(np.asarray(o, dtype=float), (n,)) for o in outs[1:]], axis=1)
        if grad.shape[0] != Y.shape[0]:
            raise ValueError(f"Forminkonsistenz bei {name}: {grad.shape} vs {Y.shape}")
        return val, grad

    # --- mpmath-Route (Referenz) -------------------------------------------
    def eval_mp(self, name: str, y):
        """Wert und Gradient (ein Punkt, hohe Praezision) aus dem Originalausdruck."""
        return self._mp[name](*list(y), *self.const_mp)


# ---------------------------------------------------------------------------
# Klammern und Vektorfelder aus Gradienten (reine Arithmetik)
# ---------------------------------------------------------------------------
def pb_np(gF: np.ndarray, gG: np.ndarray) -> np.ndarray:
    """Elementweise Poisson-Klammer aus Gradientenfeldern (N,10)."""
    return np.sum(gF[:, :5] * gG[:, 5:], axis=1) - np.sum(gF[:, 5:] * gG[:, :5], axis=1)


def coord_bracket_np(gF: np.ndarray) -> np.ndarray:
    """(N,10)-Feld mit Eintrag j = {z_j, F}; nutzt grad(z_j) = e_j."""
    out = np.empty_like(gF)
    out[:, :5] = gF[:, 5:]
    out[:, 5:] = -gF[:, :5]
    return out


def pb_dirac_np(gF, gG, g1, g2, W) -> np.ndarray:
    """Dirac-Klammer zum Paar (chi1,chi2) aus Gradienten.

    omega = [[0,W],[-W,0]] mit W = {chi1,chi2} hat die Inverse [[0,-1/W],[1/W,0]],
    also {F,G}_D = {F,G} - ({F,chi2}{chi1,G} - {F,chi1}{chi2,G})/W. Damit gilt
    {chi_a,G}_D = 0 identisch, was im Manuskript die starke Halierung der beiden
    zweiten-klassigen Constraints bedeutet.
    """
    base = pb_np(gF, gG)
    corr = (pb_np(gF, g2) * pb_np(g1, gG) - pb_np(gF, g1) * pb_np(g2, gG)) / W
    return base - corr


def dirac_pair_np(layer: FastLayer, Y: np.ndarray, chi1: str, chi2: str, ham: str) -> dict:
    """Algebra des zweiten-klassigen Paares (chi1, chi2) an N Punkten.

    Retourneiert W = {chi1,chi2}, B1 = {chi1,H}, B2 = {chi2,H}, das naive Feld
    z_dot_j = {z_j,H} und das Dirac-Feld z_dot_j = {z_j,H}_D.
    """
    _, gH = layer.eval_np(ham, Y)
    _, g1 = layer.eval_np(chi1, Y)
    _, g2 = layer.eval_np(chi2, Y)
    W = pb_np(g1, g2)
    B1 = pb_np(g1, gH)
    B2 = pb_np(g2, gH)
    naive = coord_bracket_np(gH)
    corr = (coord_bracket_np(g2) * B1[:, None] - coord_bracket_np(g1) * B2[:, None]) / W[:, None]
    return {"W": W, "B1": B1, "B2": B2, "naiv": naive, "dirac": naive - corr}


def rk4_np(field, Y0: np.ndarray, steps: int, tau: float) -> np.ndarray:
    """Vektorierte RK4: Startfeld (N,10), Feld-Funktion (N,10) -> Endpunkt (N,10)."""
    Y = np.array(Y0, dtype=float, copy=True)
    h = float(tau) / float(steps)
    f1 = field(Y)
    for _ in range(steps):
        f2 = field(Y + 0.5 * h * f1)
        f3 = field(Y + 0.5 * h * f2)
        f4 = field(Y + h * f3)
        Y = Y + h * (f1 + 2 * f2 + 2 * f3 + f4) / 6.0
        f1 = field(Y)
    return Y
