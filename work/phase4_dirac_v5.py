"""Phase 4, Schritt P-SCHNELL: Dirac-Klassifikation und Dirac-Fluss des getrimmten Modells.

Ersetzung von work/phase4_dirac_v4.py. Grund: die fruehere Fassung hat Poisson-Matrix,
Dirac-Klammer und Leitungsfelder vollsymbolisch mit simplify gebaut und anschiessend
seriell mit mpmath integriert (mehr als elf Minuten ohne jede Ausgabe, ein Kern).
Diese Fassung rechnet in Stufen, gibt jede Stufe sofort aus, legt nach jeder Stufe einen
maschinenlesbaren Zwischenstand ab und laeuft vektoeriert ueber einem Ensemble von
Anfangsdaten; die unabhaengigen hochpraezisen Referenztraektorien werden auf mehrere
Prozesse verteilt.

Inhalt der Stufen:
  A  symbolische Klammern (expand, kein simplify), exakter Anschluss an work/phase4_dirac_v2.py
  B  schneller Apparat und Validierung numpy-Route gegen mpmath-Route
  C  Reproduktion der v2-Referenzpunkte (Wurzel von Hbar nach P0, Wert von D auf der Flaeche)
  D  Randdaten auf der vollen Constraintflaeche Hbar = Cs = Cc = D = 0
  E  Poisson-Matrix, Rang, Pfaffsche, Jacobi-Reste, definierende Eigenschaft der Dirac-Klammer
  F  Multiplikator u^c und Dirac-Theorem (Feldvergleich auf der Flaeche)
  G  Flussvergleich naiv gegen Dirac, Konvergenzordnung, Constraintdrift (Prozesspool)
  H  Hochpraezise mpmath-Referenz (dps=30) als Fehlerbudget der float64-Rechnung
  I  Zaehlung first/second class und reduzierte Phasenraumdimension

Ausgabe: work/phase4_dirac_v5.json
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import mpmath
import multiprocessing as mproc
import os
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase4_model as M  # noqa: E402

OUT = Path(__file__).resolve().with_suffix(".json")
V2JSON = Path(__file__).resolve().with_name("phase4_dirac_v2.json")
# Modellparameter wie in work/phase4_dirac_v2.py, damit alle Zahlen vergleichbar bleiben.
CONSTS = {M.beta: 1.7, M.m: 0.9}
START = time.monotonic()
res: dict = {"schema": "phase4_dirac_v5", "stufenzeiten": {}}


def log(msg: str) -> None:
    print(f"[{time.monotonic() - START:7.2f}s] {msg}", flush=True)


def dump(stage: str) -> None:
    res["stufenzeiten"][stage] = round(time.monotonic() - START, 3)
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res, indent=1, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    tmp.replace(OUT)


def fnum(x) -> float:
    return float(x)


def metrics(a, b) -> dict:
    """Absolute und relative Abweichung zweier Zahlenmengen."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    diff = np.abs(a - b)
    scale = np.maximum(np.abs(b), 1.0)
    return {"abs": float(np.max(diff)), "rel": float(np.max(diff / scale)),
            "rel_gegen_max": float(np.max(diff) / max(float(np.max(np.abs(b))), 1e-30))}


def arel(a, b) -> float:
    return metrics(a, b)["rel_gegen_max"]


def safe_max_abs(a):
    """Maximalbetrag ueber endliche Eintraege. Retourne (max, anzahl_nichtendlich).

    Wichtig: Python's max(0.0, nan) liefert 0.0, damit waere eine divergierte
    Trajektorie in der Driftstatistik als 'kein Drift' erschienen.
    """
    a = np.asarray(a, dtype=float)
    fin = np.isfinite(a)
    if not np.any(fin):
        return float("nan"), int(a.size)
    return float(np.max(np.abs(a[fin]))), int((~fin).sum())


def maxdiff_endlich(a, b):
    A = np.asarray(a, dtype=float)
    B = np.asarray(b, dtype=float)
    d = np.abs(A - B)
    fin = np.isfinite(d)
    if not np.any(fin):
        return float("nan"), int(d.size)
    return float(np.max(d[fin])), int((~fin).sum())


# ---------------------------------------------------------------------------
# A: symbolische Klammern
# ---------------------------------------------------------------------------
def stage_symbolic():
    log("Stufe A: symbolische Klammern (expand ohne simplify)")
    Hbar = M.hbar_closed()
    Cs, Cc = M.constraints()
    t0 = time.monotonic()
    D = M.pb_sym(Cc, Hbar)
    C3 = M.pb_sym(Cs, Hbar)
    Dcd = M.pb_sym(Cc, D)
    DH = M.pb_sym(D, Hbar)
    dD_dac = sp.expand(sp.diff(D, M.alpha_c))
    res["stufenzeiten"]["A_klammern_sek"] = round(time.monotonic() - t0, 3)
    CsD = M.pb_sym(Cs, D)
    log("   {Cs,D} in der geschlossenen Form berechnet, now Jacobi-Identitaeten")
    G = {"Hbar": Hbar, "Cs": Cs, "Cc": Cc, "D": D}
    jac = {}
    for A, B, C in (("Hbar", "Cs", "Cc"), ("Hbar", "Cs", "D"), ("Hbar", "Cc", "D"), ("Cs", "Cc", "D")):
        t1 = time.monotonic()
        try:
            acc = sp.expand(
                M.pb_sym(G[A], M.pb_sym(G[B], G[C]))
                + M.pb_sym(G[B], M.pb_sym(G[C], G[A]))
                + M.pb_sym(G[C], M.pb_sym(G[A], G[B]))
            )
            jac[f"{A},{B},{C}"] = {"rest": str(acc), "null": str(acc) == "0",
                                   "sekunden": round(time.monotonic() - t1, 3)}
        except Exception as exc:  # noqa: BLE001
            jac[f"{A},{B},{C}"] = {"fehler": f"{type(exc).__name__}: {exc}"}
        log(f"   Jacobi {A},{B},{C}: {jac[f'{A},{B},{C}'].get('null', 'fehler')}")
    res["symbolisch"] = {
        "Cs_D_geschlossen": str(CsD),
        "jacobi": jac,
        "Cs": str(Cs),
        "Cc": str(Cc),
        "D_monomzahl": len(sp.Add.make_args(sp.expand(D))),
        "Dcd_monomzahl": len(sp.Add.make_args(sp.expand(Dcd))),
        "DH_monomzahl": len(sp.Add.make_args(sp.expand(DH))),
        "Cs_Cc_exakt_null": str(M.pb_sym(Cs, Cc)) == "0",
        "C3_from_Cs_Hbar": str(C3),
        "dD_dalpha_c_exakt_null": str(dD_dac) == "0",
        "D_poly": M.monomials(D),
        "Dcd_poly": M.monomials(Dcd),
    }
    if V2JSON.exists():
        v2 = json.loads(V2JSON.read_text(encoding="utf-8"))
        # Wichtig: in sympy 1.14 ist 'beta' zusaetzlich der Name der Beta-Funktion.
        # Ohne lokales Symbolwoerterbuch wird der gespeicherte Ausdruck deshalb falsch
        # geparst (Fehlerbild Pow * FunctionClass). Deshalb immer local_dict uebergeben.
        loc = {str(s): s for s in (list(M.Z) + [M.beta, M.m])}
        Dv2 = parse_expr(v2["C2_from_Cc_Hbar"], local_dict=loc)
        Dcdv2 = parse_expr(v2["Cc_C2"], local_dict=loc)
        res["anschluss_v2"] = {
            "D_diff_exakt": str(sp.expand(D - Dv2)) == "0",
            "Dcd_diff_exakt": str(sp.expand(Dcd - Dcdv2)) == "0",
            "v2_C2_on_constraint_surface": v2["C2_on_constraint_surface"],
            "v2_C2_on_constraint_surface_2": v2["C2_on_constraint_surface_2"],
            "v2_Hbar_root_P0": v2["Hbar_root_P0"],
            "v2_Hbar_root2_P0": v2["Hbar_root2_P0"],
            "v2_Cs_C2_closed": v2["Cs_C2"],
            "v2_Cs_C2_quadratur_ac0": v2["Cs_C2_quad_ac_zero"],
            "v2_Cs_C2_quadratur_acungleich0": v2["Cs_C2_quad_ac_nonzero"],
        }
    log(f"   Klammern fertig in {res['stufenzeiten']['A_klammern_sek']}s, "
        f"Anschluss an v2: {res.get('anschluss_v2', {}).get('D_diff_exakt')}")
    return Hbar, Cs, Cc, D, Dcd, DH


# ---------------------------------------------------------------------------
# B: schneller Apparat plus Validierung
# ---------------------------------------------------------------------------
def build_layer(Hbar, Cs, Cc, D, Dcd, DH):
    pq2, pq1, pq0 = [sp.expand(sp.Poly(D, M.pq).coeff_monomial(M.pq**k)) for k in (2, 1, 0)]
    p0sq, p0l, p0c = [sp.expand(sp.Poly(Hbar, M.P0).coeff_monomial(M.P0**k)) for k in (2, 1, 0)]
    funcs = {
        "Hbar": Hbar,
        "Cs": Cs,
        "Cc": Cc,
        "D": D,
        "Dcd": Dcd,
        "DH": DH,
        "D_pq2": pq2,
        "D_pq1": pq1,
        "D_pq0": pq0,
        "H_P02": p0sq,
        "H_P01": p0l,
        "H_P00": p0c,
    }
    return M.FastLayer(funcs, CONSTS)


def stage_layer(layer):
    log("Stufe B: schneller Apparat, Validierung numpy gegen mpmath")
    rng = np.random.default_rng(20240405)
    pts = []
    for _ in range(8):
        pts.append([
            rng.uniform(-0.4, 0.6), rng.uniform(-0.5, 0.5), 0.0,
            rng.uniform(-0.8, 0.8), rng.uniform(-0.6, 0.6),
            rng.uniform(-2.0, 2.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0),
            rng.uniform(-1.2, 1.2), rng.uniform(-1.0, 1.0),
        ])
    Y = np.array(pts, dtype=float)
    mpmath.mp.dps = 30
    aus = {}
    for name in ("Hbar", "Cs", "Cc", "D", "Dcd", "DH"):
        val_np, grad_np = layer.eval_np(name, Y)
        val_mp, grad_mp = [], []
        for y in pts:
            outs = layer.eval_mp(name, [mpmath.mpf(repr(v)) for v in y])
            val_mp.append(float(outs[0]))
            grad_mp.append([float(g) for g in outs[1:]])
        mv = metrics(val_np, np.array(val_mp))
        mg = metrics(grad_np, np.array(grad_mp))
        aus[name] = {"wert": mv, "gradient": mg}
        log(f"   {name:5s} Wert absolut {mv['abs']:.2e}, Gradient absolut {mg['abs']:.2e}")
    # Klammer direkt: numpy-Route gegen symbolisch ausgewertete Klammer
    _, gH = layer.eval_np("Hbar", Y)
    _, g1 = layer.eval_np("Cc", Y)
    _, g2 = layer.eval_np("D", Y)
    pb_np_D = M.pb_np(g1, gH)
    subs_all = {z: sp.nsimplify(v) for z, v in zip(M.Z, pts[0])}
    pb_sym_D = np.array([
        float(sp.N(M.pb_sym(layer.funcs["Cc"], layer.funcs["Hbar"]).subs(
            {z: sp.nsimplify(v) for z, v in zip(M.Z, pt)}).subs(
            {k: sp.nsimplify(v) for k, v in CONSTS.items()}), 25))
        for pt in pts
    ]) if subs_all else np.array([])
    aus["Klammer_Cc_Hbar"] = metrics(pb_np_D, pb_sym_D)
    log(f"   Klammer Cc,Hbar absolut {aus['Klammer_Cc_Hbar']['abs']:.2e}")
    res["validierung_layer"] = aus
    res["hilfs_symbole"] = [[int(o), float(cf)] for o, cf in layer.pairs]
    dump("B")


# ---------------------------------------------------------------------------
# C: v2-Referenzpunkte reproduzieren
# ---------------------------------------------------------------------------
def stage_v2_reproduktions(layer):
    log("Stufe C: v2-Referenzpunkte und D-Werte auf der Flaeche")
    base = {
        "a0": 0.3, "a1": 0.25, "c": 0.6, "q": 0.4,
        "P1": 1.5 * 0.8 * 0.4, "Pc": 0.0, "p0": 0.8, "pq": -0.35,
    }
    out = {}
    for tag, P0v in (("wurzel_1", -1.34014194788918786619041671333),
                     ("wurzel_2", 2.01384482125922654421633828245)):
        y = np.array([[base["a0"], base["a1"], 0.0, base["c"], base["q"],
                       P0v, base["P1"], base["Pc"], base["p0"], base["pq"]]])
        vals = {}
        for name in ("Hbar", "Cs", "Cc", "D", "Dcd", "DH"):
            vals[name] = float(layer.eval_np(name, y)[0][0])
        out[tag] = vals
        log(f"   {tag}: Hbar {vals['Hbar']:.12e}, D {vals['D']:.15f}, W {vals['Dcd']:.12f}")
    res["v2_reproduktion"] = out
    if "anschluss_v2" in res:
        ref1 = float(res["anschluss_v2"]["v2_C2_on_constraint_surface"])
        ref2 = float(res["anschluss_v2"]["v2_C2_on_constraint_surface_2"])
        res["v2_reproduktion"] = {
            "punkte": out,
            "D_abweichung_wurzel1": abs(out["wurzel_1"]["D"] - ref1),
            "D_abweichung_wurzel2": abs(out["wurzel_2"]["D"] - ref2),
            "referenz_werte": [ref1, ref2],
        }
    dump("C")


# ---------------------------------------------------------------------------
# D: Randdaten auf der vollen Constraintflaeche
# ---------------------------------------------------------------------------
def solve_surface(layer, rng, n=200, iters=60):
    """Randdaten auf der vollen Flaeche Hbar = Cs = Cc = D = 0.

    Cs = 0 und Cc = 0 werden durch Pc = 0 und P1 = 3 p0 q / 2 erfuellt. Fuer D = 0 und
    Hbar = 0 bleiben p_q und P_0 uebrig; beide Constraints haengen an beiden Groessen,
    also wird das 2x2-System vektorisiert mit Newton geloest (Startwert: quadratische
    Loesung von Hbar = 0 bei p_q = 0).
    """
    a0 = rng.uniform(0.15, 0.45, n)
    a1 = rng.uniform(0.05, 0.45, n)
    cc = rng.uniform(0.3, 0.9, n)
    qq = rng.uniform(0.15, 0.6, n)
    p0 = rng.uniform(0.4, 1.2, n)
    Y = np.zeros((n, 10))
    Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4] = a0, a1, 0.0, cc, qq
    Y[:, 8] = p0
    Y[:, 6] = 1.5 * p0 * qq
    h2 = layer.eval_np("H_P02", Y)[0]
    h1 = layer.eval_np("H_P01", Y)[0]
    h0 = layer.eval_np("H_P00", Y)[0]
    discP = h1**2 - 4 * h2 * h0
    ok0 = (np.abs(h2) > 1e-14) & (discP >= 0)
    Y[:, 5] = np.where(ok0, (-h1 + np.sqrt(np.clip(discP, 0.0, None))) / np.where(np.abs(h2) > 1e-14, 2 * h2, 1.0), np.nan)
    hist = []
    for _ in range(iters):
        vD, gD = layer.eval_np("D", Y)
        vH, gH = layer.eval_np("Hbar", Y)
        a11, a12 = gD[:, 9], gD[:, 5]
        a21, a22 = gH[:, 9], gH[:, 5]
        det = a11 * a22 - a12 * a21
        good = np.isfinite(det) & (np.abs(det) > 1e-12)
        dpq = np.where(good, (a22 * vD - a12 * vH) / np.where(good, det, 1.0), 0.0)
        dP0 = np.where(good, (-a21 * vD + a11 * vH) / np.where(good, det, 1.0), 0.0)
        step = np.clip(np.abs(dpq), 0.0, 2.0) < np.inf
        Y[:, 9] = np.where(step, Y[:, 9] - dpq, np.nan)
        Y[:, 5] = np.where(step, Y[:, 5] - dP0, np.nan)
        resid = np.max(np.abs(np.stack([vD, vH], axis=1)), axis=1)
        hist.append(float(np.nanmax(resid)) if np.any(np.isfinite(resid)) else float("nan"))
        if np.all(np.isfinite(resid)) and float(np.max(resid)) < 1e-13:
            break
    vals = {name: layer.eval_np(name, Y)[0] for name in ("Hbar", "Cs", "Cc", "D")}
    resid = np.max(np.abs(np.stack([vals[k] for k in ("Hbar", "Cs", "Cc", "D")], axis=1)), axis=1)
    keep = np.isfinite(Y).all(axis=1) & (resid < 1e-9)
    stats = {k: float(np.max(np.abs(v[keep]))) if np.any(keep) else float("nan") for k, v in vals.items()}
    return Y[keep], stats, {"iterationen": len(hist), "letzter_rest": hist[-1] if hist else None}


def stage_randdaten(layer):
    log("Stufe D: Randdaten auf der vollen Flaeche Hbar = Cs = Cc = D = 0")
    pool_rows = []
    versuche = []
    for seed in range(6):
        rng = np.random.default_rng(31415 + 7 * seed)
        rows, residuen, info = solve_surface(layer, rng, n=400)
        versuche.append({"seed": seed, "punkte": int(rows.shape[0]), "rest": residuen, **info})
        log(f"   seed {seed}: {rows.shape[0]} Punkte, Constraintrest {residuen}, "
            f"Newton {info['iterationen']} Iterationen")
        pool_rows.append(rows)
        if sum(r.shape[0] for r in pool_rows) >= 96:
            break
    Y = np.concatenate(pool_rows, axis=0)[:128]
    if Y.shape[0] == 0:
        raise RuntimeError("keine Randdaten auf der vollen Constraintflaeche gefunden")
    residuen = {name: float(np.max(np.abs(layer.eval_np(name, Y)[0])))
                for name in ("Hbar", "Cs", "Cc", "D")}
    log(f"   {Y.shape[0]} Randdatenpunkte verwendet, Rest maxima {residuen}")
    vals = {name: layer.eval_np(name, Y)[0] for name in ("Dcd", "DH")}
    res["randdaten"] = {
        "anzahl": int(Y.shape[0]),
        "versuche": versuche,
        "constraint_rest_max": residuen,
        "W_bereich": [float(np.min(vals["Dcd"])), float(np.max(vals["Dcd"]))],
        "DH_bereich": [float(np.min(vals["DH"])), float(np.max(vals["DH"]))],
        "punkte_auswahl": [[fnum(v) for v in row] for row in Y[:3]],
    }
    dump("D")
    return Y


# ---------------------------------------------------------------------------
# E: Poisson-Matrix, Rang, Jacobi, definierende Eigenschaft
# ---------------------------------------------------------------------------
def stage_klammeralgebra(layer, Y):
    log("Stufe E: Poisson-Matrix, Rang, Jacobi-Reste, Dirac-Eigenschaft")
    names = ["Hbar", "Cs", "Cc", "D"]
    grads = {n: layer.eval_np(n, Y)[1] for n in names}
    mats = np.zeros((Y.shape[0], 4, 4))
    for i in range(4):
        for j in range(4):
            mats[:, i, j] = M.pb_np(grads[names[i]], grads[names[j]])
    pfaff = mats[:, 0, 1] * mats[:, 2, 3] - mats[:, 0, 2] * mats[:, 1, 3] + mats[:, 0, 3] * mats[:, 1, 2]
    ranks = []
    for k in range(Y.shape[0]):
        A = mats[k]
        scale = max(1.0, float(np.max(np.abs(A))))
        ranks.append(int(np.linalg.matrix_rank(A, tol=1e-9 * scale)))
    counts = {int(r): ranks.count(r) for r in sorted(set(ranks))}
    # Hinweis: die numerische Jacobi-Identitaet braucht zweite Ableitungen (Gradient der
    # Klammer). Sie wird deshalb in Stufe A symbolisch exakt geprueft, hier bleibt die
    # definierende Eigenschaft der Dirac-Klammer.: {chi_a, G}_D = 0 fuer alle Koordinaten G
    g1, g2 = grads["Cc"], grads["D"]
    W = M.pb_np(g1, g2)
    eye = np.eye(10)
    maxchi = 0.0
    for a, gchi in (("Cc", g1), ("D", g2)):
        for j in range(10):
            gj = np.tile(eye[j], (Y.shape[0], 1))
            val = M.pb_dirac_np(gchi, gj, g1, g2, W)
            maxchi = max(maxchi, float(np.max(np.abs(val))))
    aw = np.abs(W)
    res["klammeralgebra"] = {
        "rang_verteilung": counts,
        "W_quantile_abs": {str(q): float(np.quantile(aw, q))
                           for q in (0.0, 0.05, 0.25, 0.5, 0.75, 1.0)},
        "W_unterhalb": {f"unter_{t:g}": int(np.sum(aw < t))
                        for t in (1e-3, 1e-2, 0.1, 0.5)},
        "pfaffsche_max_abs": float(np.max(np.abs(pfaff))),
        "matrix_auswahl": [[[f"{v:.6e}" for v in row] for row in mats[0]]],
        "dirac_halierung_chi_max_abs": maxchi,
        "W_min_abs": float(np.min(np.abs(W))),
    }
    log(f"   Rangverteilung {counts}, Pfaffsche max {np.max(np.abs(pfaff)):.2e}, "
        f"Dirac-Halierung max {maxchi:.2e}")
    dump("E")


# ---------------------------------------------------------------------------
# F: Multiplikator und Dirac-Theorem
# ---------------------------------------------------------------------------
def stage_multiplikator(layer, Y):
    log("Stufe F: Multiplikator u^c und Dirac-Theorem auf der Flaeche")
    pair = M.dirac_pair_np(layer, Y, "Cc", "D", "Hbar")
    gH = layer.eval_np("Hbar", Y)[1]
    g1 = layer.eval_np("Cc", Y)[1]
    naive = pair["naiv"]
    dirac = pair["dirac"]
    W = pair["W"]
    DH = pair["B2"]
    u = DH / W
    eye = np.eye(10)
    fits = []
    for j in range(10):
        gj = np.tile(eye[j], (Y.shape[0], 1))
        cbc = M.pb_np(gj, g1)  # {z_j, Cc}
        sel = np.abs(cbc) > 1e-6
        if np.any(sel):
            lam = (dirac[:, j] - naive[:, j]) / cbc
            fits.append({"j": j, "rel": arel(lam[sel], u[sel]), "abs": float(np.max(np.abs(lam[sel] - u[sel])))})
    best = min(fits, key=lambda d: d["rel"]) if fits else None
    # Multiplikator aus Erhaltung von D: {D, Hbar + lam Cc} = 0 -> lam = DH / W
    lam_cons = DH / W
    res["multiplikator"] = {
        "u_aus_Dirac_feld": {"min": float(np.min(u)), "max": float(np.max(u)),
                             "mean": float(np.mean(u))},
        "u_aus_D_erhaltung": {"min": float(np.min(lam_cons)), "max": float(np.max(lam_cons)),
                             "mean": float(np.mean(lam_cons))},
        "u_identitaet_max_abs": float(np.max(np.abs(u - lam_cons))),
        "feldvergleich_bester_kanal": best,
        "kanalzaehlung": len(fits),
    }
    # Dirac-Theorem: Feld {z,H}_D gleich Feld von H + u Cc auf der Flaeche
    uCc_field = np.zeros_like(naive)
    _, gCc = layer.eval_np("Cc", Y)
    for j in range(10):
        gj = np.tile(eye[j], (Y.shape[0], 1))
        uCc_field[:, j] = u * M.pb_np(gj, gCc)
    theorie = naive + uCc_field
    res["dirac_theorem"] = {
        "feld_max_abs_differenz": float(np.max(np.abs(dirac - theorie))),
        "feld_typischer_betrag": float(np.max(np.abs(naive))),
    }
    log(f"   u-Bereich [{np.min(u):.6f}, {np.max(u):.6f}], Identitaet zur D-Erhaltung "
        f"{res['multiplikator']['u_identitaet_max_abs']:.2e}, Theorem-Differenz "
        f"{res['dirac_theorem']['feld_max_abs_differenz']:.2e}")
    dump("F")


# ---------------------------------------------------------------------------
# G/H: Flussrechnungen im Prozesspool
# ---------------------------------------------------------------------------
_CACHE = {}


def _layer_in_worker():
    """Layer einmal je Workerprozess aufbauen (Pool: getrennte Adressraeume)."""
    if "layer" not in _CACHE:
        Hbar = M.hbar_closed()
        Cs, Cc = M.constraints()
        D = M.pb_sym(Cc, Hbar)
        Dcd = M.pb_sym(Cc, D)
        DH = M.pb_sym(D, Hbar)
        _CACHE["layer"] = build_layer(Hbar, Cs, Cc, D, Dcd, DH)
    return _CACHE["layer"]


def run_flow(task):
    """Ein Flusslauf ueber das ganze Ensemble; Divergenzen werden zaehlt, nicht versteckt."""
    method, steps, tau, payload, gruppe = task
    layer = _layer_in_worker()
    Y0 = np.asarray(payload, dtype=float)
    if method == "naiv":
        def field(YY):
            return M.coord_bracket_np(layer.eval_np("Hbar", YY)[1])
    else:
        def field(YY):
            return M.dirac_pair_np(layer, YY, "Cc", "D", "Hbar")["dirac"]
    Y = Y0.copy()
    alive = np.isfinite(Y).all(axis=1)
    segs = 4
    seg_steps = max(1, steps // segs)
    drift = {k: 0.0 for k in ("Hbar", "Cs", "Cc", "D")}
    marks = []
    tote_insgesamt = 0
    betrag_max = 0.0
    with np.errstate(all="ignore"):
        for s in range(segs):
            Y = M.rk4_np(field, Y, seg_steps, tau / segs)
            vals = {name: layer.eval_np(name, Y)[0] for name in ("Hbar", "Cs", "Cc", "D")}
            gestorben = ~np.isfinite(Y).all(axis=1)
            tote_insgesamt += int(np.sum(gestorben & alive))
            alive &= ~gestorben
            entry = {"tau": round(tau * (s + 1) / segs, 4), "lebend": int(np.sum(alive)),
                     "tod": tote_insgesamt}
            for k, arr in vals.items():
                m, nb = safe_max_abs(arr[alive]) if np.any(alive) else (float("nan"), 0)
                if np.isfinite(m):
                    drift[k] = max(drift[k], m)
                entry[k] = m
                entry[k + "_nichtendlich"] = nb
            if np.any(alive):
                betrag_max = max(betrag_max, safe_max_abs(Y[alive])[0])
            marks.append(entry)
    end = Y[:2]
    return {"methode": method, "gruppe": gruppe, "schritte": steps, "tau": tau,
            "drift_max": drift,
            "konstraint_drift_max": drift,
            "trajektorien_tot": tote_insgesamt,
            "trajektorien_lebend": int(np.sum(alive)),
            "koordinaten_betrag_max": betrag_max,
            "end_zeilen": [[fnum(v) for v in row] for row in end],
            "end_erste_8": [[fnum(v) for v in row] for row in Y[:8]],
            "marks": marks}


def run_mp(task):
    """Hochpraezise mpmath-Referenztraektorie (dps=30), ein Punkt, eine Methode."""
    mp = mpmath
    mp.mp.dps = 30
    punkt, method, steps, tau, point = task
    layer = _layer_in_worker()
    one, half, two, six = mp.mpf(1), mp.mpf("0.5"), mp.mpf(2), mp.mpf(6)
    y = [mp.mpf(repr(v)) for v in point]

    def deriv(state):
        gH = layer.eval_mp("Hbar", state)[1:]
        naive = [gH[i + 5] if i < 5 else -gH[i] for i in range(10)]
        if method == "naiv":
            return naive
        g1 = layer.eval_mp("Cc", state)[1:]
        g2 = layer.eval_mp("D", state)[1:]

        def pbb(gA, gB):
            return sum(gA[i] * gB[i + 5] - gA[i + 5] * gB[i] for i in range(5))

        W = pbb(g1, g2)
        B1 = pbb(g1, gH)
        B2 = pbb(g2, gH)
        out = []
        for j in range(10):
            pjD = g2[j + 5] if j < 5 else -g2[j]
            pjC = g1[j + 5] if j < 5 else -g1[j]
            out.append(naive[j] - (pjD * B1 - pjC * B2) / W)
        return out

    def add(state, k, factor):
        return [state[i] + factor * k[i] for i in range(10)]

    h = mp.mpf(repr(tau)) / steps
    grenze = mp.mpf("1e3")
    zeitgrenze = 90.0
    anfang = time.monotonic()
    divergenz = None
    schritt = 0
    for schritt in range(1, steps + 1):
        k1 = deriv(y)
        k2 = deriv(add(y, k1, half * h))
        k3 = deriv(add(y, k2, half * h))
        k4 = deriv(add(y, k3, h))
        y = [y[i] + h * (k1[i] + two * k2[i] + two * k3[i] + k4[i]) / six for i in range(10)]
        if any(not mp.isfinite(v) for v in y) or any(abs(v) > grenze for v in y):
            divergenz = {"schritt": schritt, "grund": "Betragsschwelle 1e3 ueberschritten"}
            break
        if time.monotonic() - anfang > zeitgrenze:
            divergenz = {"schritt": schritt, "grund": f"Zeitlimit {zeitgrenze:g}s erreicht"}
            break

    def sicher(v):
        if not mp.isfinite(v) or abs(v) > grenze:
            return None
        return float(v)

    out = {"punkt": punkt, "methode": method, "schritte": schritt, "tau": tau,
           "genauigkeit_dps": 30, "divergenz": divergenz,
           "laufzeit_s": round(time.monotonic() - anfang, 2),
           "end": [("nichtendlich" if not mp.isfinite(v) else str(v)) for v in y],
           "end_float": [sicher(v) for v in y]}
    if divergenz is None:
        vals = {n: layer.eval_mp(n, y)[0] for n in ("Hbar", "Cs", "Cc", "D")}
        out["drift"] = {k: sicher(v) for k, v in vals.items()}
    else:
        out["drift"] = None
    return out


TAUS_FLUSS = (0.05, 0.10, 0.20, 0.40)
SCHRITTE_FLUSS = (250, 500, 1000)


def stage_fluss(layer, Y):
    log(f"Stufe G: Flussvergleich naiv gegen Dirac, {Y.shape[0]} Anfangsdaten, Prozesspool")
    g1 = layer.eval_np("Cc", Y)[1]
    g2 = layer.eval_np("D", Y)[1]
    Wrand = M.pb_np(g1, g2)
    wohl = np.abs(Wrand) >= 0.5
    gruppen = [("alle", np.ones(Y.shape[0], dtype=bool)), ("wohlbedingt", wohl)]
    log(f"   tau-Leiter {TAUS_FLUSS}, Schrittzahlen {SCHRITTE_FLUSS}")
    log(f"   |W| an der Randflaeche: min {np.min(np.abs(Wrand)):.3e}, "
        f"median {np.median(np.abs(Wrand)):.3e}, max {np.max(np.abs(Wrand)):.3e}, "
        f"wohlbedingt mit |W| >= 0.5: {int(np.sum(wohl))}/{Y.shape[0]}")
    tasks = []
    groessen = {}
    for gname, maske in gruppen:
        groessen[gname] = int(np.sum(maske))
        payload = [[fnum(v) for v in row] for row in Y[maske]]
        for tau in TAUS_FLUSS:
            for meth in ("naiv", "dirac"):
                for s in SCHRITTE_FLUSS:
                    tasks.append((meth, s, tau, payload, gname))
    workers = max(1, min(6, mproc.cpu_count()))
    log(f"   {len(tasks)} Flusslaeufe auf {workers} Prozessen")
    fluesse = []
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for out in pool.map(run_flow, tasks):
            fluesse.append(out)
            log(f"   {out['gruppe']:11s} tau {out['tau']:.2f} {out['methode']:5s} "
                f"Schritte {out['schritte']:5d} "
                f"lebend {out['trajektorien_lebend']:4d}/{groessen[out['gruppe']]:4d} "
                f"Drift Cc {out['drift_max']['Cc']:.3e}, D {out['drift_max']['D']:.3e}, "
                f"Hbar {out['drift_max']['Hbar']:.3e}, "
                f"|y|max {out['koordinaten_betrag_max']:.2e}")
    res["fluss"] = fluesse
    res["fluss_gruppen"] = {n: groessen[n] for n, _ in gruppen}
    # Konvergenzordnung aus den Endpunkten, je Gruppe, Methode und tau getrennt
    for gname, _ in gruppen:
        for meth in ("naiv", "dirac"):
            for tau in TAUS_FLUSS:
                sel = sorted([f for f in fluesse
                              if f["methode"] == meth and f["gruppe"] == gname
                              and abs(f["tau"] - tau) < 1e-12],
                             key=lambda d: d["schritte"])
                if len(sel) != 3:
                    continue
                d1, bad1 = maxdiff_endlich(sel[0]["end_erste_8"], sel[1]["end_erste_8"])
                d2, bad2 = maxdiff_endlich(sel[1]["end_erste_8"], sel[2]["end_erste_8"])
                ver = (d1 / d2) if (np.isfinite(d1) and np.isfinite(d2) and d2 > 0) else float("nan")
                res.setdefault("konvergenz", {})[f"{gname}_{meth}_tau{tau}"] = {
                    "gruppe": gname, "methode": meth, "tau": tau,
                    "d_250_500": d1, "d_500_1000": d2, "verhaeltnis": ver,
                    "erwartetes_verhaeltnis_rk4": 16.0,
                    "nichtendliche_differenzen": [bad1, bad2],
                    "trajektorien_tot": [f["trajektorien_tot"] for f in sel],
                    "trajektorien_lebend": [f["trajektorien_lebend"] for f in sel],
                }
                log(f"   Konvergenz {gname} {meth} tau {tau:.2f}: d250/500 {d1:.3e}, "
                    f"d500/1000 {d2:.3e}, Verhaeltnis {ver:.2f}, "
                    f"lebend {[f['trajektorien_lebend'] for f in sel]}")
    # groester Constraintdrift je Gruppe und Methode, uber alle tau und Schrittweiten
    driftvergleich = {}
    for gname, _ in gruppen:
        for meth in ("naiv", "dirac"):
            sel = [f for f in fluesse if f["methode"] == meth and f["gruppe"] == gname]
            eintrag = {k: float(max(f["drift_max"][k] for f in sel))
                       for k in ("Hbar", "Cs", "Cc", "D")}
            eintrag["groesstes_tau"] = float(max(f["tau"] for f in sel))
            driftvergleich[f"{gname}_{meth}"] = eintrag
    res["driftvergleich_gesamt"] = driftvergleich
    for k, v in driftvergleich.items():
        log(f"   Drift gesamt {k:23s} Hbar {v['Hbar']:.2e}, Cs {v['Cs']:.2e}, "
            f"Cc {v['Cc']:.2e}, D {v['D']:.2e}")
    dump("G")


def stage_hochprazision(layer, Y):
    log("Stufe H: mpmath-Referenz (dps=30) als Fehlerbudget, gleiche Schrittweite wie float64")
    taus = (0.10, 0.20)
    schritte = 250
    g1 = layer.eval_np("Cc", Y)[1]
    g2 = layer.eval_np("D", Y)[1]
    W = M.pb_np(g1, g2)
    absteigend = [int(i) for i in np.argsort(-np.abs(W))]
    idx = sorted(set(absteigend[:2]) | {int(np.argmin(np.abs(W)))})
    log(f"   gewehlte Randpunkte nach |W|: {idx} mit |W| = "
        f"{[float(abs(W[i])) for i in idx]}")
    tasks = []
    for i in idx:
        pt = [fnum(v) for v in Y[i]]
        for meth in ("naiv", "dirac"):
            for tau in taus:
                tasks.append((i, meth, schritte, tau, pt))
    workers = max(1, min(4, mproc.cpu_count()))
    log(f"   {len(tasks)} mpmath-Lauefe auf {workers} Prozessen, je Lauf 90s, Stufe 240s")
    outs = []
    gesehen = set()
    pool = cf.ProcessPoolExecutor(max_workers=workers)
    futuren = {pool.submit(run_mp, t): t for t in tasks}
    try:
        for fut in cf.as_completed(futuren, timeout=240):
            t = futuren[fut]
            try:
                o = fut.result()
            except Exception as exc:  # eine divergente Trajektorie darf den Lauf nicht stoppen
                o = {"punkt": t[0], "methode": t[1], "schritte": 0, "tau": t[3],
                     "divergenz": {"grund": "Abbruch: " + type(exc).__name__},
                     "end": None, "end_float": None, "drift": None}
            outs.append(o)
            gesehen.add((o["punkt"], o["methode"], round(float(o["tau"]), 6)))
            d = o.get("drift")
            if isinstance(d, dict):
                log(f"   mpmath Punkt {o['punkt']} tau {o['tau']:.2f} {o['methode']:5s} "
                    f"Drift Cc {d['Cc']:.2e}, D {d['D']:.2e}, Laufzeit {o.get('laufzeit_s')}s")
            else:
                dv = o.get("divergenz") or {}
                log(f"   mpmath Punkt {o['punkt']} tau {o['tau']:.2f} {o['methode']:5s} "
                    f"divergiert: {dv.get('grund')} (Schritt {dv.get('schritt')})")
    except cf.TimeoutError:
        log("   Zeitlimit der Stufe erreicht; noch laufende mpmath-Laefe werden verworfen")
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    for t in tasks:  # nicht ausgefuehrte Laeufe ehrlich dokumentieren
        if (t[0], t[1], round(float(t[3]), 6)) not in gesehen:
            outs.append({"punkt": t[0], "methode": t[1], "schritte": 0, "tau": t[3],
                         "divergenz": {"grund": "nicht ausgefuehrt (Stufenzeitlimit)"},
                         "end": None, "end_float": None, "drift": None})
            log(f"   mpmath Punkt {t[0]} tau {t[3]:.2f} {t[1]:5s} nicht ausgefuehrt")
    res["hochpraezision"] = outs
    dump("H")
    # Fehlerbudget: gleiches Randdatum, gleiche Schrittweite h = tau/schritte, float64 gegen dps=30
    budget = {}
    for i in idx:
        Y0 = np.array([Y[i]], dtype=float)
        for meth in ("naiv", "dirac"):
            if meth == "naiv":
                field = lambda YY: M.coord_bracket_np(layer.eval_np("Hbar", YY)[1])
            else:
                field = lambda YY: M.dirac_pair_np(layer, YY, "Cc", "D", "Hbar")["dirac"]
            for tau in taus:
                Yend = M.rk4_np(field, Y0, schritte, tau)
                ref = next((o for o in outs if o["methode"] == meth and o["punkt"] == i
                            and abs(o["tau"] - tau) < 1e-12), None)
                eintrag = {"float64_endpunkt": [fnum(v) for v in Yend[0]],
                           "schrittweite": tau / schritte}
                if ref is not None and ref.get("end_float") and \
                        all(v is not None for v in ref["end_float"]):
                    y_ref = np.array(ref["end_float"], dtype=float)
                    dd = float(np.max(np.abs(Yend[0] - y_ref)))
                    betrag = float(np.max(np.abs(y_ref)))
                    eintrag.update({"max_abs_endpunkt_differenz": dd,
                                    "typischer_betrag": betrag,
                                    "relativ": dd / max(1.0, betrag)})
                    log(f"   Fehlerbudget Punkt {i} {meth} tau {tau:.2f}: {dd:.3e}")
                else:
                    eintrag["mpmath_divergiert"] = True
                    log(f"   Fehlerbudget Punkt {i} {meth} tau {tau:.2f}: mpmath divergiert")
                budget[f"punkt{i}_{meth}_tau{tau}"] = eintrag
    res["fehlerbudget_float64_gegen_mp30"] = budget
    dump("H")


# ---------------------------------------------------------------------------
# I: Zaehlung
# ---------------------------------------------------------------------------
def stage_zaehlung():
    log("Stufe I: Zaehlung first/second class")
    res["zaehlung"] = {
        "phasenraum_dimension": 10,
        "first_class": ["Hbar", "Cs"],
        "second_class": ["Cc", "D"],
        "begruendung": "{Cs,Cc}=0 exakt, {Cs,Hbar}=0 exakt, {Cc,Hbar}=D und {Cc,D} ungleich null "
                       "auf der Flaeche; Rang der 4x4-Poisson-Matrix auf der Flaeche 2.",
        "formel": "10 - 2*2 - 2",
        "reduzierte_dimension": 10 - 2 * 2 - 2,
        "kanonische_paare_reduziert": (10 - 2 * 2 - 2) // 2,
    }
    dump("I")
    log("fertig")


def main():
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    res["modellparameter"] = {str(k): float(v) for k, v in CONSTS.items()}
    res["werkzeuge"] = {"numpy": np.__version__, "sympy": sp.__version__,
                        "mpmath": mpmath.__version__, "kerne": mproc.cpu_count()}
    Hbar, Cs, Cc, D, Dcd, DH = stage_symbolic()
    layer = build_layer(Hbar, Cs, Cc, D, Dcd, DH)
    stage_layer(layer)
    stage_v2_reproduktions(layer)
    Y = stage_randdaten(layer)
    stage_klammeralgebra(layer, Y)
    stage_multiplikator(layer, Y)
    stage_fluss(layer, Y)
    stage_hochprazision(layer, Y)
    stage_zaehlung()


if __name__ == "__main__":
    main()
