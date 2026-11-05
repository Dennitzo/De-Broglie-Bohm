"""Gegenpruefung der in Abschnitt 19.8 gedruckten Zahlen gegen work/phase4_dirac_v5.json.

Das Skript aendert nichts, sondern druckt die tatsaechlich gespeicherten Werte in
derselben Reihenfolge, wie sie im Manuskriptabschnitt verwendet werden. Damit laesst
sicht pruefen, ob Manuskript und Rechnung uebereinstimmen (Akzeptanzkriterium ac6).
"""

from __future__ import annotations

import json
from pathlib import Path

PFAD = Path("work/phase4_dirac_v5.json")


def zeige_struktur(obj, vor="", tiefe=0):
    if tiefe > 1:
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                print(f"{vor}{k}: <{type(v).__name__}{' len=' + str(len(v)) if isinstance(v, list) else ''}>")
                zeige_struktur(v, vor + "  ", tiefe + 1)
            else:
                print(f"{vor}{k} = {v!r}")
    elif isinstance(obj, list):
        if obj and isinstance(obj[0], dict):
            print(f"{vor}[0] Schluessel: {sorted(obj[0].keys())}")
        else:
            print(f"{vor}Beispiele: {obj[:6]!r}")


def main() -> None:
    dat = json.loads(PFAD.read_text(encoding="utf-8"))
    print("=== Top-Level-Schluessel ===")
    print(sorted(dat.keys()))
    print()
    for schl in sorted(dat.keys()):
        print(f"=== {schl} ===")
        zeige_struktur(dat[schl], "  ")
        print()
    print("=== Publikationsstatus ===")
    text = Path("Publikation.md").read_text(encoding="utf-8", newline="")
    marke = "### 19.8 Dirac-Klammer"
    print("enthalt 19.8:", marke in text)
    print("CRLF:", "\r\n" in text, "Zeilen:", text.count("\n"), "Zeichen:", len(text))


if __name__ == "__main__":
    main()
