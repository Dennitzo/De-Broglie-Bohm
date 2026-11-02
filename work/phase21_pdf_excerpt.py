"""Unabhaengige Textextraktion des gerenderten Publikations-PDFs (Kapitel 21).

Zweck: Belegen, dass das neue Kapitel 21 (Globale zeitliche Vollstaendigkeit
und Kollapsanalyse des gefuehrten Skalenfaktors) tatsaechlich im erzeugten PDF
enthalten ist, also nicht nur in der Markdown-Quelle. Dies ersetzt die
unzuverlaessige Fensteraufnahme des Betrachters. Es wird kein inhaltlicher
Beweis geliefert, sondern allein die Praesenz des Kapitels, seiner
Unterabschnitte, seiner Tabelle und der Skript- und Lean-Referenzen im
gerenderten Dokument nachgewiesen. Schluesselwoerter sind umlautarm gewaehlt,
damit die Textextraktion zuverlaessig trifft.
"""
import sys
from PyPDF2 import PdfReader

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

pdf = sys.argv[1] if len(sys.argv) > 1 else "Publikation.pdf"
reader = PdfReader(pdf)
pages = reader.pages
total = len(pages)

texts = []
for pg in pages:
    texts.append(pg.extract_text() or "")

# Umlautarme Schluesselwoerter aus Kapitel 21.
keys = [
    "Globale zeitliche",
    "Fortsetzungssatz",
    "linearem Wachstum",
    "Kollapszeit",
    "Gronwall",
    "Superlinear",
    "Endlichkeit",
    "phase21_global_existence",
    "GlobalExistenceBarrier",
]

hit_pages = {}
for k in keys:
    hit_pages[k] = [i + 1 for i, t in enumerate(texts) if k in t]

print("pdf", pdf)
print("total_pages", total)
for k in keys:
    print("key", k, "->", hit_pages.get(k, []))

# Kapitel-21-Unterabschnitte (auto-nummeriert, Punktformat im Haupttext).
sub = {
    "21.1": [i + 1 for i, t in enumerate(texts) if "21.1 Fortsetzungssatz" in t],
    "21.2": [i + 1 for i, t in enumerate(texts) if "21.2 Wesentlichkeit" in t],
    "21.3": [i + 1 for i, t in enumerate(texts) if "21.3 Br" in t],
    "21.5": [i + 1 for i, t in enumerate(texts) if "21.5 Ergebnis und Reichweite" in t],
}
for kk, vv in sub.items():
    print("section", kk, "->", vv)

lean_ref = [i + 1 for i, t in enumerate(texts) if "work/lean/GlobalExistenceBarrier.lean" in t]
skript_ref = [i + 1 for i, t in enumerate(texts) if "analysis/phase21_global_existence.py" in t]
print("lean_referenz_seiten", lean_ref)
print("skript_referenz_seiten", skript_ref)

kapitel_haupttext = bool(sub["21.1"]) and bool(hit_pages["Fortsetzungssatz"])
kapitel_ueberschrift = bool(hit_pages["Globale zeitliche"])
belege = bool(lean_ref) and bool(skript_ref)
tabelle = bool(sub["21.5"])
gesamt = kapitel_haupttext and kapitel_ueberschrift and belege and tabelle
print("KAPITEL_21_IM_PDF", "JA" if gesamt else "NEIN")
print("  kapitel_haupttext", kapitel_haupttext)
print("  kapitel_ueberschrift", kapitel_ueberschrift)
print("  skript_und_lean_referenz", belege)
print("  ergebnis_tabelle", tabelle)
sys.exit(0 if gesamt else 4)
