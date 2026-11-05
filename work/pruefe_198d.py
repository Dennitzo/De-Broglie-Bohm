"""Erzeugt die Tabellenwerte fuer Abschnitt 19.8 direkt aus phase4_dirac_v5.json.

Ausgegeben werden je Methode und Ensemble:
  * lebende Trajektorien bei 250, 500 und 1000 Schritten,
  * die Constraintdrifts bei derselben Schrittanzahl,
  * das Konvergenzverhaeltnis.
Zusaetzlich die Hochpraezisionslaeufe mit Punktindex, tau, Methode, Drift und Divergenz.
"""

from __future__ import annotations

import json
from pathlib import Path

dat = json.loads(Path("work/phase4_dirac_v5.json").read_text(encoding="utf-8"))

tab = {}
for e in dat["fluss"]:
    schluessel = (e["methode"], e["gruppe"], e["tau"])
    tab.setdefault(schluessel, {})[e["schritte"]] = e

for methode in ("naiv", "dirac"):
    for gruppe in ("alle", "wohlbedingt"):
        for tau in (0.05, 0.1, 0.2, 0.4):
            schl = (methode, gruppe, tau)
            if schl not in tab:
                continue
            print(f"--- {methode} / {gruppe} / tau={tau}")
            for schritte in (250, 500, 1000):
                e = tab[schl].get(schritte)
                if e is None:
                    print(f"   {schritte}: fehlt")
                    continue
                d = e["drift_max"]
                print(f"   {schritte}: lebend={e['trajektorien_lebend']} tot={e['trajektorien_tot']} "
                      f"betrag_max={e['koordinaten_betrag_max']:.4g} "
                      f"drift Hbar={d['Hbar']:.4g} Cs={d['Cs']:.4g} Cc={d['Cc']:.4g} D={d['D']:.4g}")

print()
print("--- Konvergenzverhaeltnisse (Erwartung 16 fuer RK4)")
for schl, bl in sorted(dat["konvergenz"].items()):
    print(f"{schl}: verhaeltnis={bl['verhaeltnis']:.4f} "
          f"lebend={bl['trajektorien_lebend']} tot={bl['trajektorien_tot']} "
          f"nichtendlich={bl['nichtendliche_differenzen']}")

print()
print("--- Hochpraezisionslaeufe")
for e in dat["hochpraezision"]:
    d = e["drift"]
    if d is None:
        print(f"punkt {e['punkt']} tau {e['tau']} {e['methode']}: divergenz={e['divergenz']}")
        continue
    print(f"punkt {e['punkt']} tau {e['tau']} {e['methode']} schritte={e['schritte']} "
          f"laufzeit={e['laufzeit_s']}s drift Hbar={d['Hbar']:.4g} Cs={d['Cs']:.4g} "
          f"Cc={d['Cc']:.4g} D={d['D']:.4g} divergenz={e['divergenz']}")

print()
print("--- Fehlerbudget float64 gegen mpmath (30 Stellen), identische Schrittweite")
for schl, bl in sorted(dat["fehlerbudget_float64_gegen_mp30"].items()):
    if bl.get("mpmath_divergiert"):
        print(f"{schl}: mpmath divergiert, Schrittweite={bl['schrittweite']}")
    else:
        print(f"{schl}: abs={bl['max_abs_endpunkt_differenz']:.4g} rel={bl['relativ']:.4g} "
              f"h={bl['schrittweite']} typisch={bl['typischer_betrag']:.4g}")
