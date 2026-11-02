"""Unabhaengige Textextraktion des gerenderten Publikations-PDFs.

Zweck: Belegen, dass das neue Kapitel 20 (Austauschbilanz, Quellenkonsistenz
und Rueckwirkung fuer n Felder) tatsaechlich im erzeugten PDF enthalten ist,
also nicht nur in der Markdown-Quelle. Dies ersetzt die unzuverlaessige
Fensteraufnahme des Betrachters. Es wird kein inhaltlicher Beweis geliefert,
sondern allein die presence des Kapitels und seiner Unterabschnitte sowie der
Lean-Referenz im gerenderten Dokument nachgewiesen.
"""
import sys
from PyPDF2 import PdfReader

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

pdf = sys.argv[1] if len(sys.argv) > 1 else "Publikation.pdf"
reader = PdfReader(pdf)
pages = reader.pages
total = len(pages)

keys = [
    "Austauschbilanz",
    "Quellenkonsistenz",
    "Einzelkonfigurationsquelle",
    "Divergenzfreiheit",
    "R\u00fcckwirkung",
    "ExchangeBalanceCore",
]

hit_pages = {}
texts = []
for i, pg in enumerate(pages):
    t = pg.extract_text() or ""
    texts.append(t)
    for k in keys:
        if k in t:
            hit_pages.setdefault(k, []).append(i + 1)

print("pdf", pdf)
print("total_pages", total)
for k in keys:
    print("key", k, "->", hit_pages.get(k, []))

# Kapitel 20 Haupttext- und Inhaltsverzeichnisseite(n).
toc_pages = [i + 1 for i, t in enumerate(texts) if "20 Austauschbilanz" in t]
sec203 = [i + 1 for i, t in enumerate(texts) if "20.3 Die verallgemeinerte Austauschbilanz" in t]
sec204 = [i + 1 for i, t in enumerate(texts) if "20.4 Quellenkonsistenz" in t]
lean_ref = [i + 1 for i, t in enumerate(texts) if "work/lean/ExchangeBalanceCore.lean" in t]

print("kapitel20_kopf_oder_gliederung_seiten", toc_pages)
print("abschnitt_20_3_seiten", sec203)
print("abschnitt_20_4_seiten", sec204)
print("lean_referenz_seiten", lean_ref)

ok = bool(toc_pages and sec203 and sec204)
print("KAPITEL_20_IM_PDF", "JA" if ok else "NEIN")
sys.exit(0 if ok else 3)
