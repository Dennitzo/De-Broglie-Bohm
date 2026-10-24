# De-Broglie-Bohm

Dieses Repository enthält die deutschsprachige Fachpublikation zur de-Broglie-Bohm-Theorie sowie ausgewählte Python-Modelle, Lean-Beweise und Abbildungen.

## Publikation

- [Publikation als Markdown](Publikation.md)
- [Publikation als PDF](Publikation.pdf)

Die Arbeit behandelt Bohmsche Teilchen- und Felddynamik, relativistische Modelle sowie ausgewählte Anwendungen in Quantenfeldtheorie, Gravitation und Kosmologie. Die Aussagen gelten jeweils innerhalb der in der Publikation beschriebenen Modelle und Voraussetzungen. Numerische Simulationen und formale Beweise ersetzen keinen allgemeinen physikalischen oder empirischen Nachweis.

## Begleitmaterial

- `work/`: freigegebene Python-Modelle und Lean-Beweismodule unter `work/lean/`.
- `figures/`: Abbildungen und Plots als PNG- und SVG-Dateien.
- `Bohmian_Field_Model.py`: ein endliches Bohmsches Zweimoden-Feldmodell.

Die freigegebene Dateiauswahl enthält keine lokalen Quellenkopien, Rohdaten, Forschungsnotizen, Prüfberichte oder Build-Artefakte. Literaturangaben und fachliche Einordnung stehen in der Publikation.

## Python-Skripte

Die Python-Dateien sind eigenständige Modelle und Berechnungen; ein gemeinsamer Projekt-Runner ist nicht enthalten. Für das Feldmodell installiere Python und NumPy, SciPy sowie SymPy in einer lokalen Umgebung und starte es vom Repository-Hauptverzeichnis:

```powershell
python -m pip install numpy==2.5.2 scipy==1.17.1 sympy==1.14.0
python Bohmian_Field_Model.py
```

Andere Skripte können zusätzliche lokale Eingabedaten voraussetzen. Solche Daten sind nicht Bestandteil dieses Repositories.

## Lean-Beweise

Die Lean-Module befinden sich in `work/lean/` und verwenden Mathlib. Die lokale Lake-Konfiguration und Build-Artefakte werden nicht mitveröffentlicht; zum erneuten Kompilieren ist eine passende Lean-/Mathlib-Umgebung erforderlich.
