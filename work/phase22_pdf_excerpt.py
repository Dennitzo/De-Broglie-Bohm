"""Unabhaengige Textextraktion des gerenderten Publikations-PDFs (Kapitel 22)."""
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

keys = [
    "Periodenbedingung",
    "knotenfreien",
    "Wirbelmodell",
    "Topologie",
    "Einheitskreis",
    "Phanomen",
    "phase22_periods",
    "work/lean/Periods.lean",
]

hit_pages = {}
for k in keys:
    hit_pages[k] = [i + 1 for i, t in enumerate(texts) if k in t]

print("pdf", pdf)
print("total_pages", total)
for k in keys:
    print("key", k, "->", hit_pages.get(k, []))

sub = {
    "22.1": [i + 1 for i, t in enumerate(texts) if "22.1 Topologie" in t],
    "22.2": [i + 1 for i, t in enumerate(texts) if "22.2 Herleitung" in t],
    "22.3": [i + 1 for i, t in enumerate(texts) if "22.3 Formaler" in t],
    "22.4": [i + 1 for i, t in enumerate(texts) if "22.4 Ergebnis" in t],
}
for kk, vv in sub.items():
    print("section", kk, "->", vv)

kapitel_haupttext = bool(sub["22.1"]) and bool(hit_pages["Periodenbedingung"])
kapitel_ueberschrift = bool(hit_pages["Periodenbedingung"])
lean_ref = [i + 1 for i, t in enumerate(texts) if "work/lean/Periods.lean" in t]
skript_ref = [i + 1 for i, t in enumerate(texts) if "phase22_periods" in t]
belege = bool(lean_ref) and bool(skript_ref)
tabelle = bool(sub["22.4"])
gesamt = kapitel_haupttext and kapitel_ueberschrift and belege and tabelle
print("KAPITEL_22_IM_PDF", "JA" if gesamt else "NEIN")
print("  kapitel_haupttext", kapitel_haupttext)
print("  kapitel_ueberschrift", kapitel_ueberschrift)
print("  skript_und_lean_referenz", belege)
print("  ergebnis_tabelle", tabelle)
sys.exit(0 if gesamt else 4)
