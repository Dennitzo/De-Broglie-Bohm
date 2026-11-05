"""Zusammenfassung der Fluss- und Hochpraezisionslaeufe aus phase4_dirac_v5.json.

Die Manuskripttabellen fuer lebende Trajektorien, Constraintdrifts und die
Hochpraezisionsgegenpruefung sollenexakt aus diesen gespeicherten Werten entstehen.
"""

from __future__ import annotations

import json
from pathlib import Path

dat = json.loads(Path("work/phase4_dirac_v5.json").read_text(encoding="utf-8"))

print("### fluss: Eintraege", len(dat["fluss"]))
print("Schluessel eines Eintrags:", sorted(dat["fluss"][0].keys()))
for ein in dat["fluss"]:
    kurz = {}
    for k, v in ein.items():
        if isinstance(v, (int, float, str)) or v is None:
            kurz[k] = v
        elif isinstance(v, dict):
            kurz[k] = {kk: (f"{vv:.4g}" if isinstance(vv, float) else vv) for kk, vv in v.items()}
    print(kurz)

print()
print("### hochpraezision: Eintraege", len(dat["hochpraezision"]))
print("Schluessel eines Eintrags:", sorted(dat["hochpraezision"][0].keys()))
for ein in dat["hochpraezision"]:
    kurz = {}
    for k, v in ein.items():
        if isinstance(v, list) and len(v) > 4:
            kurz[k] = f"<liste {len(v)}>"
        elif isinstance(v, dict):
            kurz[k] = {kk: (f"{vv:.4g}" if isinstance(vv, float) else vv) for kk, vv in v.items()}
        else:
            kurz[k] = v
    print(kurz)

print()
print("### Validierung (Wert/Gradient numpy gegen mpmath)")
for name, block in sorted(dat["validierung_layer"].items()):
    zeilen = {}
    for art, d in sorted(block.items()):
        if isinstance(d, dict):
            zeilen[art] = {kk: (f"{vv:.3g}" if isinstance(vv, float) else vv) for kk, vv in d.items()}
        else:
            zeilen[art] = d
    print(name, zeilen)
