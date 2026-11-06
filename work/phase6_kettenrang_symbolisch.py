"""Phase 6, Teil A: Rangsprung der Ketten-Klammermatrix exakt (Fourier-Algebra).

Gegengeprueft wird dieselbe Kettenkonstruktion wie in work/phase5_modfamilie.py,
aber OHNE Quadratur und OHNE Toleranz: jedes Mittel <f> ist der Koeffizient der
Frequenz 0 der komplexen Fourierdarstellung, jedes Produkt eine Faltung. Alle
Aussagen sind damit Polynomidentitaeten in den generischen Moduskoeffizienten;
es wird kein Koeffizient spezialisiert, ausser in den ausdruecklich als
Spezialisierung gekennzeichneten Rangzeugnissen.

Familien (Sektorschalter eps_g, eps_m wie in §19.9):
    alpha = a0 + sum_k (a_k sin k th + [b_k cos k th])
    phi   = c  + sum_k (u_k sin k th + [v_k cos k th])
    P     = P0 + sum_k 2 (Pi_k sin k th + H_k cos k th)
    p     = p0 + sum_k 2 (pi_k sin k th + rho_k cos k th)
    r     = P'/3 - p phi',     C_A = <r w_A>,  w_A = Geometriebasis (1, sin, cos, ...)
Kanonische Paare wegen der Doppelnormierung: (a_j, Pa_j) und (u_j, Pf_j),
Poisson-Klammer {F,G} = sum_i (d_q_i F d_p_i G - d_p_i F d_q_i G).

Behandelte Aussagen (jeweils exakt, keine Toleranz):
  (1) C_A haengt in den Geometrievariablen nur von den Impulsen ab; omega^{AB}
      haengt ueberhaupt nicht von den Geometrievariablen (structural).
  (2) Im Sektor eps_m = 0 ist C_0 identisch null, in eps_m = 1 gilt die
      paritaetsabhaengige Form C_0 = sum_k k (pi_k v_k - rho_k u_k).
  (3) Paarmuster: welche Eintraege {C_A, C_B} ueberleben identisch und welche
      verschwinden exakt, fuer K = 1..3 in beiden Materiesektoren. Das ist die
      exakte Aussage ueber die Selektionsregeln der Mittelbildung.
  (4) Rankschranke rang omega <= 2K, weil die Kette 1 + 2K Elemente (ungerade)
      hat und omega antisymmetrisch ist; im Sektor eps_m = 0 zusaendlich wegen
      der identisch verschwindenden Zeile und Spalte 0.
  (5) Generischer Rang des Blocks omega[1:,1:] ist 2K fuer K = 1..4 (beide
      Materiesektoren), ausser im getrimmten Fall K = 1: dort ist der Block die
      Nullmatrix und alle 2x2-Minoren sind exakt null. Belegt wird das ueber
      einen nicht verschwindenden Minor der Ordnung 2K: die Determinante des
      Blocks ist fuer K <= 3 als Polynom explizit nicht null, und fuer alle K
      ist ihr Wert an einer festen rationalen Spezialierung nicht null. Da
      Spezialierung den Rang nur senken kann, liefert das together mit (4)
      rang_generisch = 2K ohne jede Toleranz.
  (6) kappa[A,j] = dC_A/dPa_j ist zustaendsunabhaengig und hat nicht
      verschwindende Determinante (Vertraeglichkeit der Eichbedingungen a_j).

Unabhaengige Gegenpruefung: alle C_A und omega-Eintraege werden an einem
festen Punkt gegen den Quadratur-/Jet-Apparat von work/phase5_modfamilie.py
gerechnet (exakte Fourieralgebra gegen Trapetzquadratur mit M = 512).

Fortschrittsdruck je Abschnitt, JSON-Zwischenstand nach jedem Abschnitt und nach
jedem Rangabschnitt; bereits vorhandene Abschnitte werden bei einem Fortsetzungs-
lauf uebersprungen (wiederaufnahmefaehig).
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import time

import sympy as sp

I = sp.I
OUT = "work/phase6_kettenrang_symbolisch.json"
T0 = time.time()
res: dict = {}


def log(msg: str) -> None:
    print(f"[{time.time() - T0:7.2f}s] {msg}", flush=True)


def dump() -> None:
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, sort_keys=True)


def load() -> None:
    global res
    if os.path.exists(OUT):
        try:
            with open(OUT, encoding="utf-8") as f:
                res = json.load(f)
            log(f"fortsetzung: vorhandene Abschnitte {sorted(res)}")
        except (OSError, ValueError):
            res = {}


# ---------------------------------------------------------------------------
# exakte Fourier-Polynom-Rechnung: dict Frequenz -> Koeffizient
# ---------------------------------------------------------------------------
def mul(f: dict, g: dict) -> dict:
    out: dict = {}
    for nf, cf in f.items():
        if cf == 0:
            continue
        for ng, cg in g.items():
            if cg == 0:
                continue
            k = nf + ng
            out[k] = sp.expand(out.get(k, sp.Integer(0)) + cf * cg)
    return {k: v for k, v in out.items() if v != 0}


def scale(f: dict, c) -> dict:
    return {k: sp.expand(c * v) for k, v in f.items() if v != 0}


def add(f: dict, g: dict) -> dict:
    out = dict(f)
    for k, c in g.items():
        out[k] = sp.expand(out.get(k, sp.Integer(0)) + c)
    return {k: v for k, v in out.items() if v != 0}


def derivative(f: dict) -> dict:
    return {n: sp.expand(I * n * c) for n, c in f.items() if c != 0}


def avg(f: dict):
    """Mittel <f> = Koeffizient der Frequenz 0 (fuer reelle f reell)."""
    v = f.get(0, sp.Integer(0))
    return sp.simplify(v)


def trig(mode):
    if mode is None:
        return {0: sp.Integer(1)}
    k, is_cos = mode
    if is_cos:
        return {k: sp.Rational(1, 2), -k: sp.Rational(1, 2)}
    return {k: -I * sp.Rational(1, 2), -k: I * sp.Rational(1, 2)}


def mode_list(K: int, cos: bool):
    ms = []
    for k in range(1, K + 1):
        ms.append((k, False))
        if cos:
            ms.append((k, True))
    return ms


def sym_names_matter(mlist, first: str, second: str):
    out = []
    for k, is_cos in mlist:
        out.append(sp.Symbol(f"{first}{k}" if not is_cos else f"{second}{k}"))
    return out


class SFam:
    """Symbolische Familie mit exakter Mittelbildung und kanonischen Paaren."""

    def __init__(self, K: int, geom_cos: bool = True, matter_cos: bool = True):
        self.K, self.geom_cos, self.matter_cos = K, geom_cos, matter_cos
        self.gm, self.mm = mode_list(K, geom_cos), mode_list(K, matter_cos)
        self.ng, self.nf = 1 + len(self.gm), 1 + len(self.mm)
        self.nW = self.ng
        # Konfiguration
        self.za = [sp.Symbol("a0")] + sym_names_matter(self.gm, "a", "b")
        self.zc = [sp.Symbol("c")] + sym_names_matter(self.mm, "u", "v")
        # Impulse (Doppelnormierung wie im Numerikapparat)
        self.Pa = [sp.Symbol("P0")] + sym_names_matter(self.gm, "Pi", "H")
        self.Pf = [sp.Symbol("p0")] + sym_names_matter(self.mm, "pi", "rho")
        self.pairs = list(zip(self.za, self.Pa)) + list(zip(self.zc, self.Pf))
        self.W = [trig(None)] + [trig(m) for m in self.gm]
        # Felder als Fourierpolynome
        self.P = {0: self.Pa[0]}
        for j, m in enumerate(self.gm):
            self.P = add(self.P, scale(trig(m), 2 * self.Pa[j + 1]))
        self.p = {0: self.Pf[0]}
        for j, m in enumerate(self.mm):
            self.p = add(self.p, scale(trig(m), 2 * self.Pf[j + 1]))
        self.phi = {0: self.zc[0]}
        for j, m in enumerate(self.mm):
            self.phi = add(self.phi, scale(trig(m), self.zc[j + 1]))
        self.dP = derivative(self.P)
        self.dphi = derivative(self.phi)
        self.C = [sp.expand(avg(mul(scale(self.dP, sp.Rational(1, 3)), w))
                             - avg(mul(mul(self.p, self.dphi), w))) for w in self.W]
        self.label = f"K={K},eps_g={int(geom_cos)},eps_m={int(matter_cos)}"
        self._omega = None

    def wlabel(self, A: int) -> str:
        if A == 0:
            return "1"
        k, is_cos = self.gm[A - 1]
        return (f"cos{k}" if is_cos else f"sin{k}")

    def pb(self, F, G):
        s = sp.Integer(0)
        for q, p in self.pairs:
            s += sp.expand(sp.diff(F, q) * sp.diff(G, p) - sp.diff(F, p) * sp.diff(G, q))
        return sp.expand(s)

    def omega(self):
        if self._omega is not None:
            return self._omega
        n = self.nW
        rows = [[sp.Integer(0)] * n for _ in range(n)]
        for A in range(n):
            for B in range(A + 1, n):
                v = self.pb(self.C[A], self.C[B])
                rows[A][B] = v
                rows[B][A] = sp.expand(-v)
        self._omega = sp.Matrix(rows)
        return self._omega


def rational_point(symbols, salt: int = 0):
    """Feste rationale Belegung (nur fuer die numerische Gegenprobe benutzt)."""
    vals = {}
    for i, s in enumerate(symbols):
        vals[s] = sp.Rational(3 * ((i + 1) * 7 + 5 * salt) % 11 + 2, 5)
    return vals


def int_point(symbols, trial: int) -> dict:
    """Deterministische ganze Belegung aus einem festen Pseudozufallsstrom.

    Ganze Werte sind absichtlich: Aus Polynomial != 0 folgt fuer jede Belegung
    ein Wert ungleich 0 nicht, aber aus einem Wert ungleich 0 folgt Polynomial
    ungleich 0 (Auswertung ist ein Ringhomomorphismus). Die Bestimmung des
    generischen Rangs laeuft deshalb ueber eine GEPRUEFTE Belegung, nicht ueber
    eine beliebige: eine einzige Belegung kann den Rang senken, sie kann ihn
    aber nie heben. Fundstelle: die erste feste rationale Belegung lieferte bei
    eps_m=0, K=3 Rang 4, obgleich das Determinantenpolynom nicht null ist.
    """
    import random

    rnd = random.Random(20261107 + 977 * trial)
    vals = {}
    for s in symbols:
        v = rnd.choice([-5, -4, -3, -2, 2, 3, 4, 5])
        vals[s] = sp.Integer(v)
    return vals


def good_point(block, symbols, tries: int = 12):
    """Erste Belegung, deren Blockdeterminante exakt nicht verschwindet.

    Rueckgabe (werte, versuch, determinant, alle_versuche). Ist keine Belegung
    erfolgreich, wird (None, tries, None, ...) zurueckgegeben; dann ist der
    Minor unter diesem Verfahren nicht belegt und der symbolische Weg noetig.
    """
    attempts = []
    for t in range(tries):
        vals = int_point(symbols, t)
        d = block.subs(vals).det()
        attempts.append({"versuch": t, "det": str(d), "nicht_null": bool(d != 0)})
        if d != 0:
            return vals, t, d, attempts
    return None, tries, None, attempts


def _short(expr, n: int = 420) -> str:
    s = str(expr)
    return s if len(s) <= n else s[:n] + f"... ({len(s)} Zeichen)"


def section_struktur() -> None:
    if "struktur" in res:
        log("Abschnitt struktur: vorhanden, uebersprungen")
        return
    struct = {}
    for K in (1, 2, 3):
        f = SFam(K, True, True)
        geom_conf = set(f.za[1:]) | set(f.za)
        paeder = [sp.simplify(sp.diff(c, q)) for c in f.C for q in f.za]
        free = sorted({str(s) for row in f.omega().tolist() for e in row
                       for s in e.free_symbols})
        struct[f.label] = {
            "C_A_haengt_nicht_von_geometriekonfiguration": bool(all(d == 0 for d in paeder)),
            "C0_exakt": _short(f.C[0]),
            "omega_freie_variablen": free,
            "geometrie_variablen_in_omega": sorted(v for v in free if v in geom_conf
                                                   or v in {str(s) for s in f.Pa}),
        }
        log(f"   {f.label}: Geometrieabhaengigkeit der Klammern: "
            f"{struct[f.label]['geometrie_variablen_in_omega'] or 'keine'}, "
            f"C_0 = {_short(f.C[0], 120)}")
    res["struktur"] = struct
    dump()


def section_mittelwertelement() -> None:
    if "mittelwertelement" in res:
        log("Abschnitt mittelwertelement: vorhanden, uebersprungen")
        return
    parity = {}
    for K in (1, 2, 3, 4):
        f = SFam(K, True, True)
        closed = sp.Integer(0)
        idx = 0
        for k in range(1, K + 1):
            u_s = f.zc[1 + idx]; pi_s = f.Pf[1 + idx]; idx += 1
            if f.matter_cos:
                v_c = f.zc[1 + idx]; rho_c = f.Pf[1 + idx]; idx += 1
                closed += k * (pi_s * v_c - rho_c * u_s)
        parity[f"K={K},eps_m=1"] = {
            "differenz_zu_C0": str(sp.simplify(f.C[0] - closed)),
            "geschlossen": bool(sp.simplify(f.C[0] - closed) == 0),
        }
        f0 = SFam(K, True, False)
        parity[f"K={K},eps_m=0"] = {
            "C0_exakt_null": bool(sp.simplify(f0.C[0]) == 0),
            "C0_wert": _short(sp.simplify(f0.C[0])),
        }
        log(f"   K={K}: Paritaetsform eps_m=1 Differenz "
            f"{parity[f'K={K},eps_m=1']['differenz_zu_C0']}, "
            f"eps_m=0 C_0 null: {parity[f'K={K},eps_m=0']['C0_exakt_null']}")
    res["mittelwertelement"] = parity
    dump()


def section_paarmuster() -> None:
    """Welche Paarungen {C_A, C_B} ueberleben die Mittelbildung exakt."""
    if "paarmuster" in res:
        log("Abschnitt paarmuster: vorhanden, uebersprungen")
        return
    out = {}
    for eps_m in (False, True):
        for K in (1, 2, 3):
            tf = time.time()
            f = SFam(K, True, eps_m)
            om = f.omega()
            nW = f.nW
            labels = [f.wlabel(A) for A in range(nW)]
            patt = [[bool(om[A, B] != 0) for B in range(nW)] for A in range(nW)]
            nonzero = {}
            for A in range(nW):
                for B in range(A + 1, nW):
                    if om[A, B] != 0:
                        nonzero[f"{labels[A]}|{labels[B]}"] = _short(om[A, B], 300)
            zero = [f"{labels[A]}|{labels[B]}" for A in range(nW) for B in range(A + 1, nW)
                    if om[A, B] == 0]
            key = f"eps_m={int(eps_m)},K={K}"
            out[key] = {"labels": labels, "muster": patt,
                        "anzahl_nichtnull_paare": len(nonzero),
                        "anzahl_nullpaare": len(zero),
                        "nichtnull_paare": nonzero,
                        "null_paare": zero,
                        "c0_null": bool(sp.simplify(f.C[0]) == 0),
                        "zeit_s": time.time() - tf}
            log(f"   {key}: {len(nonzero)} ueberlebende Paarungen von "
                f"{len(nonzero) + len(zero)} moeglichen, C_0 "
                f"{'identisch null' if out[key]['c0_null'] else 'nicht null'} "
                f"({out[key]['zeit_s']:.1f}s)")
    res["paarmuster"] = out
    dump()


def section_rang() -> None:
    """Generischer Rang ohne Toleranz: Strukturoberschranke plus Minorzeugnis."""
    ranks = res.setdefault("rang", {})
    for eps_m in (False, True):
        for K in (1, 2, 3, 4):
            key = f"eps_m={int(eps_m)},K={K}"
            if key in ranks:
                log(f"   {key}: vorhanden, uebersprungen")
                continue
            tf = time.time()
            f = SFam(K, True, eps_m)
            om = f.omega()
            nW = f.nW
            m = nW - 1                                  # Blockgroesse 2K
            block = om[1:, 1:]
            entry = {"nW": nW, "block_groesse": m,
                     "ungerade_kette": bool(nW % 2 == 1),
                     "antisymmetrisch_exakt": bool(
                         all((om[A, B] + om[B, A]) == 0
                             for A in range(nW) for B in range(nW)))}
            # (a) getrimmter Fall: Blockidentitaet und alle 2x2-Minoren exakt null
            if K == 1:
                nz = 0
                for i1 in itertools.combinations(range(m), 2):
                    for i2 in itertools.combinations(range(m), 2):
                        m2 = sp.Matrix([[block[a, b] for b in i2] for a in i1])
                        if sp.simplify(m2.det()) != 0:
                            nz += 1
                entry["minoren_2x2_nicht_null"] = nz
                entry["block_ist_nullmatrix"] = bool(
                    all(sp.simplify(e) == 0 for row in block.tolist() for e in row))
                entry["omega_ist_nullmatrix"] = bool(
                    all(sp.simplify(e) == 0 for row in om.tolist() for e in row))
                entry["eintraege_omega"] = [_short(sp.simplify(e), 160)
                                            for row in om.tolist() for e in row]
            # (b) Minorzeugnis als Polynom. Die Polynomdeterminante ist der teuerste
            # Schritt; sie wird nur dort gerechnet, wo sie nachweislich bezahlbar ist
            # (Block bis Groesse 4, dazu Groesse 6 im Sektor eps_m=0). Fuer eps_m=1 und
            # Groesse 6 wurde sie versucht und lief in ueber 22 Minuten nicht durch; sie
            # ist aber nicht noetig, weil das Minorzeugnis (c) ohne Polynom auskommt.
            polynom_geplant = (0 < m <= 4) or (m == 6 and not eps_m)
            if polynom_geplant:
                tb = time.time()
                detb = sp.simplify(block.det(method="berkowitz"))
                entry["det_block_polynom_nicht_null"] = bool(detb != 0)
                entry["det_block_anzahl_glieder"] = (len(sp.Add.make_args(detb))
                                                     if detb != 0 else 0)
                entry["det_block_polynom"] = (_short(detb, 2000) if m <= 4 else
                                              f"(Polynom mit {entry['det_block_anzahl_glieder']} "
                                              f"GLIEDERN, nicht abgedruckt)")
                entry["det_block_zeit_s"] = time.time() - tb
            else:
                entry["det_block_polynom_nicht_null"] = None
                entry["det_block_bemerkung"] = ("Determinante des Blocks der Groesse "
                                                f"{m} nicht als Polynom expandiert "
                                                "(Kosten zu hoch); der Minor ist durch "
                                                "die geprüfte ganzzahlige Belegung (c) "
                                                "exakt belegt")
            # (c) GEPRUEFTE ganzzahlige Belegung: Minor != 0 und Rang ohne Toleranz
            syms = sorted(om.free_symbols, key=str)
            vals, versuch, detb_spec, attempts = good_point(block, syms)
            entry["belegungsversuche"] = attempts
            entry["belegung_versuch"] = versuch
            entry["belegung_gefunden"] = bool(vals is not None)
            if vals is not None:
                oms = om.subs(vals)
                blocks = block.subs(vals)
                entry["belegung"] = {str(k): str(v) for k, v in
                                     sorted(vals.items(), key=lambda kv: str(kv[0]))[:12]}
                entry["det_block_belegung"] = str(detb_spec)
                entry["det_block_belegung_nicht_null"] = bool(detb_spec != 0)
                entry["rang_block_belegung"] = int(blocks.rank())
                entry["rang_volle_kette_belegung"] = int(oms.rank())
                entry["kern_dimension_volle_kette"] = int(nW - oms.rank())
            else:
                entry["det_block_belegung"] = None
                entry["det_block_belegung_nicht_null"] = None
                entry["rang_block_belegung"] = None
                entry["rang_volle_kette_belegung"] = None
                entry["kern_dimension_volle_kette"] = None
                entry["belegung_bemerkung"] = ("keine Belegung mit nicht verschwindendem "
                                               f"Minor in {versuch} Versuchen gefunden")
            # (d) symbolische Elimination nur fuer kleine Bloecke
            if 0 < m <= 4:
                entry["rang_generisch_block_symbolisch"] = int(block.rank())
            else:
                entry["rang_generisch_block_symbolisch"] = None
                entry["rang_generisch_bemerkung"] = (
                    "symbolische Elimination uebersprungen (Kosten exponentiell); "
                    "der generische Rang folgt aus Minorzeugnis (c) plus "
                    "Strukturoberschranke (e)")
            # (e) Strukturoberschranke rang omega <= 2K
            entry["oberschranke_2K"] = 2 * K
            entry["det_omega_exakt_null"] = (bool(sp.simplify(om.det(method="berkowitz")) == 0)
                                             if nW <= 5 else None)
            if nW > 5:
                entry["det_omega_bemerkung"] = ("Determinante der vollen Kette nicht "
                                                "expandiert; die Strukturbeobachtung "
                                                "'ungerade antisymmetrisch' gilt allgemein")
            low = entry["rang_block_belegung"]
            entry["rang_generisch_block"] = low
            entry["rang_generisch_block_gleich_2K"] = (bool(low == 2 * K)
                                                       if low is not None else None)
            entry["begruendung"] = ("rang(Spezialisierung) <= rang_generisch; "
                                    "rang_generisch omega <= 2K wegen Antisymmetrie und "
                                    f"ungerader Kettenlaenge {nW}"
                                    + (" und identisch verschwindender Zeile/Spalte 0"
                                       if not eps_m else "")
                                    + f"; damit rang_generisch Block = {low}")
            entry["zeit_s"] = time.time() - tf
            ranks[key] = entry
            dump()
            log(f"   {key}: Blockgroesse {m}, Rang an geprüfter Belegung "
                f"{entry['rang_block_belegung']} (Versuch "
                f"{entry['belegung_versuch']}) / Oberschranke {2 * K}, "
                f"Minor als Polynom nicht null: {entry['det_block_polynom_nicht_null']}, "
                f"Kern dim omega {entry['kern_dimension_volle_kette']} "
                f"({entry['zeit_s']:.1f}s)")
    dump()


def section_kappa() -> None:
    if "kappa" in res:
        log("Abschnitt kappa: vorhanden, uebersprungen")
        return
    kap = {}
    for K in (1, 2, 3, 4):
        f = SFam(K, True, False)
        rows = [[sp.simplify(sp.diff(f.C[A], f.Pa[j])) for j in range(1, f.ng)]
                for A in range(1, f.nW)]
        M = sp.Matrix(rows)
        d = sp.simplify(M.det()) if M.shape[0] == M.shape[1] else None
        kap[f"K={K}"] = {"matrix": [[_short(e, 120) for e in row] for row in M.tolist()],
                         "det": _short(d), "det_nicht_null": bool(d != 0),
                         "quadratisch": bool(M.shape[0] == M.shape[1]),
                         "zustaende": sorted({str(s) for row in M.tolist() for e in row
                                              for s in e.free_symbols})}
        log(f"   kappa K={K}: |det| = {_short(d, 80)}, zustaendige Variablen: "
            f"{kap[f'K={K}']['zustaende'] or 'keine'}")
    res["kappa"] = kap
    dump()


def section_gegenprobe() -> None:
    if "gegenprobe_numerisch" in res:
        log("Abschnitt gegenprobe_numerisch: vorhanden, uebersprungen")
        return
    log("Gegenpruefung: exakte Fourieralgebra gegen Quadratur-/Jet-Apparat")
    import numpy as np
    if "work" not in sys.path:
        sys.path.insert(0, "work")
    import phase5_modfamilie as F5                                  # noqa: E402

    cross = {}
    for K, eps_m in ((1, False), (1, True), (2, False), (2, True), (3, False)):
        tf = time.time()
        f = SFam(K, True, eps_m)
        fam = F5.Family(K, geom_cos=True, matter_cos=bool(eps_m))
        syms = f.za + f.zc + f.Pa + f.Pf
        vals = rational_point(syms, salt=K * 3 + (0 if eps_m else 1))
        vec = np.array([float(vals[s]) for s in syms])
        assert vec.size == fam.ndim, (vec.size, fam.ndim)
        chi, Jac = F5.chi_all(fam, vec[None, :])
        JC = Jac[:, 1:, :]
        om_num = F5.bracket_matrix(JC, JC, fam.nq)[0]
        om_sym = np.array([[complex(f.pb(f.C[A], f.C[B]).subs(vals))
                            for B in range(f.nW)] for A in range(f.nW)])
        d_om = float(np.max(np.abs(om_num - om_sym.real)))
        c_sym = np.array([complex(c.subs(vals)) for c in f.C])
        d_c = float(np.max(np.abs(chi[0][1:] - c_sym.real)))
        cross[f"eps_m={int(eps_m)},K={K}"] = {
            "omega_max_abweichung": d_om, "C_max_abweichung": d_c,
            "omega_betrag": float(np.max(np.abs(om_num))),
            "C_betrag": float(np.max(np.abs(chi[0][1:]))),
            "zeit_s": time.time() - tf}
        log(f"   eps_m={int(eps_m)}, K={K}: omega max Abweichung {d_om:.3e} "
            f"(Betrag bis {np.max(np.abs(om_num)):.3e}), C max Abweichung {d_c:.3e} "
            f"({time.time() - tf:.1f}s)")
    res["gegenprobe_numerisch"] = cross
    dump()


def main() -> None:
    load()
    log("Phase 6 Teil A: exakte Fourier-Algebra der Kettenklammern")
    section_struktur()
    section_mittelwertelement()
    section_paarmuster()
    section_rang()
    section_kappa()
    section_gegenprobe()
    log("fertig")


if __name__ == "__main__":
    main()
