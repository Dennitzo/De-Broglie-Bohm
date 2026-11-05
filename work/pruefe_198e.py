"""Unabhaengige Nachrechnung der Randdaten und von |W| je Punkt.

Das Skript ruft die Stadiumsfunktionen des Hauptskripts auf, ohne dass dabei das
Ergebnispapier uebgeschrieben wird (dump ist auf eine Nichtstun-Funktion gesetzt).
Es prueft damit zweierlei:
  1. Reproduzierbarkeit: Constraintreste und |W|-Verteilung stimmen mit dem gespeicherten
     Papier ueberein,
  2. welche Randpunkte die Hochpraezisionslaeufe wirklich sind (kleinstes und die zwei
     groessten |W|), weil das Manuskript die Zuordnung benutzt.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase4_dirac_v5 as V5  # noqa: E402
import phase4_model as M  # noqa: E402

dat = json.loads(Path("work/phase4_dirac_v5.json").read_text(encoding="utf-8"))

# Schreibfunktion des Hauptskripts stilllegen, damit nichts uebgeschrieben wird.
V5.dump = lambda *a, **k: None
V5.log = lambda *a, **k: None

Hbar, Cs, Cc, D, Dcd, DH = V5.stage_symbolic()
layer = V5.build_layer(Hbar, Cs, Cc, D, Dcd, DH)
Y = V5.stage_randdaten(layer)

g1 = layer.eval_np("Cc", Y)[1]
g2 = layer.eval_np("D", Y)[1]
W = M.pb_np(g1, g2)
a = np.abs(W)
o = np.sort(a)
n = a.size
print("Anzahl Randdaten:", n)
print("|W| min %.6g  Quantil 0.05 %.6g  Median %.6g  max %.6g"
      % (o[0], o[max(0, int(round(0.05 * (n - 1))))], o[n // 2], o[-1]))
print("unter 0.1:", int((a < 0.1).sum()), " unter 0.5:", int((a < 0.5).sum()),
      " ab 0.5:", int((a >= 0.5).sum()))
argmin = int(np.argmin(a))
absteigend = [int(i) for i in np.argsort(-a)]
print("argmin |W| =", argmin, "mit |W| = %.6g" % a[argmin])
print("zwei groesstes |W|:", [(i, float(a[i])) for i in absteigend[:2]])
for i in (37, 45, 96):
    rang = int((a > a[i]).sum()) + 1
    print(f"Punkt {i}: |W| = {float(a[i]):.6g}, Rang {rang} von {n} "
          f"(1 = groesster Betrag)")
print("Vergleich mit Papier: W_min_abs %.6g vs nachgerechnet %.6g"
      % (dat["klammeralgebra"]["W_min_abs"], o[0]))
rest = dat["randdaten"]["constraint_rest_max"]
for name in ("Hbar", "Cs", "Cc", "D"):
    wert = float(np.max(np.abs(layer.eval_np(name, Y)[0])))
    print(f"Constraint {name}: nachgerechnet {wert:.3e}, Papier {rest[name]:.3e}")
