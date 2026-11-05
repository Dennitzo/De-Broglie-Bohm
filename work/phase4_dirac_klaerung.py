"""Klaerung der Klassifikation: welche Constraintkombination ist wirklich erstklassig?

Hintergrund. work/phase4_dirac_v5.py zaehlt Hbar und Cs als erstklassig. Die gespeicherte
Poisson-Matrix mats[0] (Reihenfolge Hbar, Cs, Cc, D) zeigt aber {Hbar,D} = -49,79 am ersten
Randpunkt, also verschwindet diese Klammer auf der Constraintflaeche NICHT. Nach Diracs
Definition ist erstklassig, wer mit samtlichen Constraints schwach verschwindende Klammern
hat. Diese Rechnung prueft numerisch:
  1. Rang der 4x4-Matrix omega und Rang des Blocks (Cc,D),
  2. die analytische Richtung des Nullraums von omega (Kombination Hbar - ({Hbar,D}/W) Cc),
  3. ob gerade die totale Hamilton-Funktion H_T = Hbar + u Cc mit u = {D,Hbar}/W diese
     erstklassige Kombination ist, indem alle vier Klammern {H_T, chi_a} auf der Flaeche
     gemessen werden,
  4. Bereich von {Hbar,D} und {D,Hbar} auf den Randdaten.
Das Ergebnis ist der Beleg fuer den Manuskripttext, nicht eine Behauptung.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase4_dirac_v5 as V5  # noqa: E402
import phase4_model as M  # noqa: E402

V5.dump = lambda *a, **k: None
V5.log = lambda *a, **k: None

Hbar, Cs, Cc, D, Dcd, DH = V5.stage_symbolic()
layer = V5.build_layer(Hbar, Cs, Cc, D, Dcd, DH)
Y = V5.stage_randdaten(layer)

N = Y.shape[0]
namen = ["Hbar", "Cs", "Cc", "D"]
g = {n: layer.eval_np(n, Y)[1] for n in namen}
W = M.pb_np(g["Cc"], g["D"])
bHD = M.pb_np(g["Hbar"], g["D"])          # {Hbar, D}
DHb = M.pb_np(g["D"], g["Hbar"])          # {D, Hbar}
u = DHb / W

print("Randdaten:", N, " Punkte, |W|-Bereich [%.6g, %.6g]" % (np.abs(W).min(), np.abs(W).max()))
print("{Hbar,D} Bereich [%.6g, %.6g], |{Hbar,D}| min %.6g" % (bHD.min(), bHD.max(), np.abs(bHD).min()))
print("{D,Hbar} Bereich [%.6g, %.6g]" % (DHb.min(), DHb.max()))
print("Multiplikator u = {D,Hbar}/W Bereich [%.6g, %.6g]" % (u.min(), u.max()))

# 1) Rang der omega und des Zweitklassenblocks
mats = np.zeros((N, 4, 4))
for i in range(4):
    for j in range(4):
        mats[:, i, j] = M.pb_np(g[namen[i]], g[namen[j]])
rangs, rangs_block = {}, {}
for k in range(N):
    scale = max(1.0, float(np.max(np.abs(mats[k]))))
    r = int(np.linalg.matrix_rank(mats[k], tol=1e-9 * scale))
    rangs[r] = rangs.get(r, 0) + 1
    rb = int(np.linalg.matrix_rank(mats[k][np.ix_([2, 3], [2, 3])], tol=1e-9 * scale))
    rangs_block[rb] = rangs_block.get(rb, 0) + 1
print("Rang omega:", rangs, " Rang Block (Cc,D):", rangs_block)

# 2) Nullraumrichtung: v = (1, 0, u, 0)?  Analyse sagt v_Cc = -(b/W) = u (weil b = {Hbar,D}).
#    Geprüft wird omega^T v = 0 bzw. omega v = 0 fuer v = (1, 0, u, 0) spaltenweise.
v = np.stack([np.ones(N), np.zeros(N), u, np.zeros(N)], axis=1)
res_uv = np.stack([np.einsum("ab,b->a", mats[k], v[k]) for k in range(N)])
res_vu = np.stack([np.einsum("a,ab->b", v[k], mats[k]) for k in range(N)])
print("Nullraumtest omega v: max |Rest| %.3e, v^T omega: %.3e"
      % (np.max(np.abs(res_uv)), np.max(np.abs(res_vu))))
skalar = np.einsum("ni,nij,nj->n", v, mats, v)
print("v^T omega v (antisymmetrisch, muss exakt 0): %.3e" % np.max(np.abs(skalar)))

# 3) Klammern von H_T = Hbar + u Cc mit allen Constraints, direkt ausgewertet.
HT = Hbar + u * Cc  # u ist eine Zahl pro Punkt, hier symbolisch nicht noetig
print()
print("Klammer {H_T, chi_a} auf der Flaeche (direkt, mit Cc-Termen):")
# Direkt in numpy: {f, chi} = grad f J grad chi. Dafuer brauchen wir grad(Hbar) + u grad(Cc)
# plus den Term Cc * grad(u) (entfaellt auf der Flaeche exakt, weil Cc = 0).
gCc = g["Cc"]
for a, name in enumerate(namen):
    feld = M.pb_np(g["Hbar"], g[name]) + u * M.pb_np(gCc, g[name])
    zusatz = layer.eval_np("Cc", Y)[0]  # Cc-Wert pro Punkt (constraintrest)
    print(f"   {name:5s}: maxabs = {np.max(np.abs(feld)):.3e}   "
          f"(Cc-Rest max {np.max(np.abs(zusatz)):.3e})")
print()
print("Vergleich: Klammer {Hbar, chi_a} ohne Zusatzterm (nur Hbar ist nicht erstklassig):")
for a, name in enumerate(namen):
    feld = M.pb_np(g["Hbar"], g[name])
    print(f"   {name:5s}: maxabs = {np.max(np.abs(feld)):.3e}")
