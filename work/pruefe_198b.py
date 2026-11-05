"""Kompakter Ausdruck der fuer Abschnitt 19.8 benoetigten Zahlen aus phase4_dirac_v5.json.

Listen mit vielen Eintraegen (Randdaten, Einzellaufen) werden nur zusammengefasst:
Anzahl, Minimum/Maximum, und die Kennung des Punktes mit dem kleinsten |W|.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

dat = json.loads(Path("work/phase4_dirac_v5.json").read_text(encoding="utf-8"))


def tief(obj, vor="", tiefe=0, maxliste=4):
    if isinstance(obj, dict):
        for k in sorted(obj.keys()):
            v = obj[k]
            if isinstance(v, dict):
                print(f"{vor}{k}: dict({len(v)})")
                if tiefe < 4:
                    tief(v, vor + "  ", tiefe + 1, maxliste)
            elif isinstance(v, list):
                zahlen = [x for x in v if isinstance(x, (int, float)) and not isinstance(x, bool)]
                if len(v) <= maxliste:
                    print(f"{vor}{k}: {v}")
                elif len(zahlen) == len(v):
                    print(f"{vor}{k}: liste({len(v)}) min={min(zahlen):.6g} max={max(zahlen):.6g}")
                else:
                    print(f"{vor}{k}: liste({len(v)})")
            else:
                print(f"{vor}{k} = {v!r}")
    else:
        print(f"{vor}{type(obj).__name__}: {obj!r}"[:300])


for schluessel in ["stufenzeiten", "klammeralgebra", "multiplikator", "konvergenz",
                   "fluss", "fluss_gruppen", "hochpraezision", "fehlerbudget_float64_gegen_mp30",
                   "randdaten", "schema", "symbolisch"]:
    if schluessel not in dat:
        print(f"### {schluessel}: fehlt")
        continue
    print(f"### {schluessel}")
    tief(dat[schluessel])
    print()

# Punktespezifisch: Minimum |W| und zugehoerige indexbasierte Befunde.
rp = dat.get("randdaten")
if isinstance(rp, dict):
    print("### Randdaten-Schluessel:", sorted(rp.keys()))
elif isinstance(rp, list):
    print("### Randdaten Liste der Laenge", len(rp), "Schluessel:", sorted(rp[0].keys()))
    ws = []
    for i, p in enumerate(rp):
        for k, v in p.items():
            if isinstance(v, (int, float)) and ("W" in str(k) or "|W|" in str(k)):
                ws.append((abs(v), i, k, v))
    if ws:
        ws.sort()
        print("kleinste |W|:", [(f"{a:.6g}", i, k) for a, i, k, _ in ws[:4]])
        print("groesste |W|:", [(f"{a:.6g}", i, k) for a, i, k, _ in ws[-2:]])
        quant = sorted(a for a, _, _, _ in ws)
        n = len(quant)
        for q in (0.0, 0.05, 0.5, 1.0):
            idx = min(n - 1, int(round(q * (n - 1))))
            print(f"Quantil {q}: {quant[idx]:.6g}")
        print("unter 0.1:", sum(1 for a in quant if a < 0.1),
              "unter 0.5:", sum(1 for a in quant if a < 0.5),
              "ab 0.5:", sum(1 for a in quant if a >= 0.5))
        print("Median:", f"{quant[n // 2]:.6g}")
