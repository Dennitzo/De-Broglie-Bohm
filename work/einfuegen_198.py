"""Fuegt die Sektion 19.8 in Publikation.md vor Kapitel 20 ein.

Die Datei enthaelt Zeilenenden CRLF; der neue Abschnitt wird an die vorhandene
Zeilenendkonvention angepasst, damit die deterministische PDF-Erzeugung ein
einheitliches Bild erhaelt. Ein zweiter Lauf ist ohne Rueckwirkung, wenn die
Sektion bereits vorhanden ist (Marke wird geprueft). Alle Zahlen im Abschnitt
stammen aus work/phase4_dirac_v5.json, work/pruefe_198.py bis work/pruefe_198d.py
und work/phase4_dirac_klaerung.py.
"""

from pathlib import Path

MARKE = "### 19.8 Dirac-Klammer, Fortpflanzung und geführter Fluss des getrimmten Systems"
KAPITEL = "## 20. Regulierter klassischer Grenzübergang, Bilanz und Konditionierung"

ABSCHNITT = r"""### 19.8 Dirac-Klammer, Fortpflanzung und geführter Fluss des getrimmten Systems

**Zweck.** §19.7 hatte als Punkt (vi) offene Rechnungsschritte benannt: die Determinante der Poisson-Matrix der Constraintmenge auf der Constraintfläche, die explizite Dirac-Klammer, die daraus folgende reduzierte Phasenraumdimension, die Fortpflanzung von $D$ und eine $\alpha_c$-vollständige geschlossene Form von $\bar{\mathcal H}$. Diese Sektion rechnet die ersten vier Punkte auf demselben getrimmten System wie §19.7 — gleiche Koordinaten, gleiche Doppelnormierung, gleiche Modellparameter $\beta=1{,}7$, $m=0{,}9$ — und führt die zugehörige Flussrechnung durch. Der fünfte Punkt bleibt offen. Tragende Aussagen sind symbolisch (geschlossene Form, keine numerische Differentiation), gegen eine unabhängige Hochpräzisionsroute geprüft und für die algebraischen Kernaussagen zusätzlich formalisiert; die Grenzen der Reichweite sind am Ende ausgewiesen.

**Constraintmenge und Klassifikation.** Gezählt wird in der zehndimensionalen Orthonormalphase $(\alpha_0,\alpha_1,\alpha_c,c,q\,;\,P_0,P_1,P_c,p_0,p_q)$ zu der Constraintmenge
$$\chi=(\bar{\mathcal H},\;C_s,\;C_c,\;D),\qquad C_s=-\tfrac13P_c,\qquad C_c=\tfrac13P_1-\tfrac12p_0q,\qquad D=\{C_c,\bar{\mathcal H}\}.$$
Auf $128$ Randdaten der vollen Fläche $\bar{\mathcal H}=C_s=C_c=D=0$ (Newton-Iteration, Constraintrest maximal $1{,}69\cdot10^{-15}$ bei $\bar{\mathcal H}$, exakt $0$ bei $C_s$, $2{,}78\cdot10^{-17}$ bei $C_c$, $1{,}57\cdot10^{-15}$ bei $D$) hat die $4\times4$-Poisson-Matrix $\omega_{ab}=\{\chi_a,\chi_b\}$ durchweg den Rang $2$, und ihre Pfaffsche Determinante ist exakt $0$; der Block $(C_c,D)$ hat an allen Punkten Rang $2$. Die Einträge sind in dieser Sektion numerisch, nicht symbolisch, also gilt die Rangangabe für die geprüfte Datenmenge und für die Scheibe $\alpha_c=0$ der geschlossenen Form.

Die Zerlegung in erst- und zweitklassig ist dabei feiner, als die Zählung allein nahelegt. Gemessen auf den Randdaten:
$$\{C_s,C_c\}=\{C_s,\bar{\mathcal H}\}=\{C_s,D\}=0,\qquad \{C_c,D\}=W\ne0,$$
$$\{\bar{\mathcal H},C_c\}=-D\approx0,\qquad \{\bar{\mathcal H},D\}\ \text{nimmt Werte in}\ [-5{,}03\cdot10^{1},\,1{,}13\cdot10^{-2}]\ \text{an},$$
mit kleinstem Betrag $1{,}13\cdot10^{-2}$, also verschwindet $\{\bar{\mathcal H},D\}$ auf der Fläche **nicht**. $\bar{\mathcal H}$ ist nach Diracs Definition nicht erstklassig, denn erstklassig heißt: mit allen Constraints schwach verschwindende Klammern. Erstklassig ist die Kombination
$$G_\star=\bar{\mathcal H}-\frac{\{\bar{\mathcal H},D\}}{W}\,C_c=\bar{\mathcal H}+\frac{\{D,\bar{\mathcal H}\}}{W}\,C_c,$$
und das ist gerade die total Hamilton-Funktion $\mathcal H_{\rm T}$ dieses Systems. Ihre vier Klammern mit den Constraints sind auf den Randdaten (schwacher Teil, der Term $C_c\{u,\chi_a\}$ ist wegen $C_c\le2{,}78\cdot10^{-17}$ nicht auflösbar) gemessen mit maximalem Betrag $3{,}33\cdot10^{-14}$ ($\bar{\mathcal H}$), exakt $0$ ($C_s$), $1{,}55\cdot10^{-15}$ ($C_c$) und $1{,}78\cdot10^{-15}$ ($D$); dieselbe Rechnung ohne den Zusatzterm $u\,C_c$ liefert für $\{\bar{\mathcal H},D\}$ den Maximalbetrag $5{,}03\cdot10^{1}$. Analytisch ist die Nullraumrichtung durch $\omega\,(1,0,u,0)^{\mathsf T}$ mit dem Rest $3{,}33\cdot10^{-14}$ gegenüber Matrixeinträgen bis $5{,}03\cdot10^{1}$ bestätigt.

Damit sind $(C_c,D)$ das zweite-klassige Paar und $\{C_s,G_\star\}$ die erstklassige Menge, und die Zählung
$$\dim_{\rm red}=10-2\cdot\underbrace{2}_{\text{erstklassig}}-\underbrace{2}_{\text{zweitklassig}}=4$$
liefert zwei kanonische Paare im reduzierten Phasenraum. Die vier Jacobi-Identitäten der Menge, $(\bar{\mathcal H},C_s,C_c)$, $(\bar{\mathcal H},C_s,D)$, $(\bar{\mathcal H},C_c,D)$ und $(C_s,C_c,D)$, sind in der geschlossenen Form symbolisch exakt null (Laufzeiten $0{,}019$ s bis $1{,}19$ s je Identität); das ist die Bedingung dafür, dass diese Klassifikation überhaupt algebraisch zulässig ist.

**Explizite Dirac-Klammer.** Mit $\chi_1=C_c$, $\chi_2=D$,
$$\omega=\begin{pmatrix}0&W\\-W&0\end{pmatrix},\qquad \omega^{-1}=\begin{pmatrix}0&-\frac1W\\[2pt]\frac1W&0\end{pmatrix},$$
ist die Dirac-Klammer auf der Scheibe $W\ne0$
$$\{F,G\}_D=\{F,G\}-\{F,\chi_a\}\,(\omega^{-1})^{ab}\,\{\chi_b,G\}
=\{F,G\}+\frac{\{F,D\}\{G,C_c\}-\{F,C_c\}\{G,D\}}{W}.$$
Die definierende Eigenschaft wurde koordinatenweise gemessen: $\{\chi_a,y^j\}_D=0$ für alle zehn Koordinaten $y^j$ und beide $a$ mit maximalem Rest $4{,}44\cdot10^{-16}$, also auf Maschinengenauigkeit. Die algebraischen Schritte dieser und der folgenden Absätze sind formalisiert in `work/phase4_dirac_klassifikation_lean.lean` (Poisson-Klammer als Struktur, Dirac-Klammer, ihre Halierung der beiden Zweitklassigen, die vier Klammern von $\mathcal H_{\rm T}$, Nichtersterstklassigkeit von $\bar{\mathcal H}$ bei $\{D,\bar{\mathcal H}\}\ne0$, und die unten angeschriebene Identität zwischen Dirac-Entwicklung und totaler Hamilton-Funktion).

**Totale Hamilton-Funktion und Fortpflanzung.** Die totale Hamilton-Funktion des getrimmten Systems ist
$$\mathcal H_{\rm T}=\bar{\mathcal H}+u\,C_c,\qquad u=\frac{\{D,\bar{\mathcal H}\}}{W},$$
und ihre Constraintfortpflanzung lautet, weil $\{C_c,C_c\}=0$ und $\{D,C_c\}=-W$,
$$\frac{dC_c}{d\tau}=\{C_c,\mathcal H_{\rm T}\}=0\ \text{stark},\qquad
\frac{dD}{d\tau}=\{D,\bar{\mathcal H}\}+u\{D,C_c\}+\{D,u\}C_c=\{D,u\}C_c,$$
$$\frac{d\bar{\mathcal H}}{d\tau}=-u\,D+C_c\{\bar{\mathcal H},u\},\qquad
\frac{dC_s}{d\tau}=C_c\{C_s,u\}.$$
Alle vier rechten Seiten liegen im von den Constraints erzeugten Ideal, $\mathcal H_{\rm T}$ ist also der Eichgenerator. Der Multiplikator ist durch die Erhaltung von $D$ festgelegt; die Identität $\{D,\bar{\mathcal H}\}-u\,W=0$ misst exakt $0$, und der aus dem Dirac-Feld zurückgerechnete Multiplikator stimmt mit dem aus der $D$-Erhaltung berechneten auf exakt denselben Wertebereich $[-31{,}908;\,29{,}653]$ überein. Wegen $u\sim W^{-1}$ enthält der Term $\{D,u\}C_c$ Potenzen von $W^{-1}$ und ist damit in den Randzeiten $|W|\to0$ verstärkt; das ist der algebraische Grund für die numerische Konditionierung unten — und, wie die Flussrechnung zeigt, nicht der einzige Mechanismus.

**Rechenapparat und seine Validierung.** Die geschlossenen Ausdrücke werden in `work/phase4_model.py` einmal symbolisch expandiert (ohne `simplify`) und dann als vektorisierte numpy-Auswertung mit den Bessel-Bestandteilen als Hilfstermen eingesetzt; dieselben Ausdrücke liegen unabhängig als mpmath-Route mit $30$ Stellen vor. Der Vergleich beider Routen über Wert und Gradient liefert absolute Reste zwischen $1{,}39\cdot10^{-17}$ und $2{,}27\cdot10^{-13}$ bei den Werten und bis $6{,}82\cdot10^{-13}$ bei den Gradienten; die Klammer $\{C_c,\bar{\mathcal H}\}$ selbst stimmt auf $1{,}87\cdot10^{-14}$. Die Reproduktion der Referenzpunkte aus §19.7 ist exakt: $D=0{,}789790412246887$ bzw. $D=-0{,}253058687154208$ an den beiden reellen Nullstellen der projizierten Hamilton-Constraint, zugehörig $W=-0{,}841400540750$ bzw. $W=-0{,}641445935174$. Die Randdaten und ihre Constraintreste sind unabhängig nachgerechnet worden (`work/phase4_dirac_klaerung.py`) und stimmen bitgenau mit dem Ergebnispapier überein.

**Konditionierung der zweiten Klassenklammer.** Über die $128$ Randdaten ist
$$|W|_{\min}=2{,}41\cdot10^{-2},\qquad |W|_{0{,}05}=0{,}161,\qquad |W|_{\rm med}=0{,}615,\qquad |W|_{\max}=2{,}425,$$
(Quantile linear interpoliert) mit $4$ Punkten unterhalb $0{,}1$ und $50$ Punkten unterhalb $0{,}5$. Die Flussrechnung wird deshalb über zwei Ensembles gefahren: alle $128$ Punkte und das wohlbedingte Teilensemble $|W|\ge0{,}5$ mit $78$ Punkten. Als Lebendzählung gilt eine Trajektorie, deren Koordinatenbeträge bis zum Ende $\tau$ unterhalb der Schwelle $10^3$ bleiben. Bei der Standardauflösung von $250$ Schritten:

| Ende $\tau$ | naiv, alle | naiv, wohlbedingt | Dirac, alle | Dirac, wohlbedingt |
|---|---|---|---|---|
| $0{,}05$ | $127/128$ | $77/78$ | $128/128$ | $78/78$ |
| $0{,}10$ | $121/128$ | $71/78$ | $127/128$ | $77/78$ |
| $0{,}20$ | $113/128$ | $65/78$ | $119/128$ | $70/78$ |
| $0{,}40$ | $83/128$ | $62/78$ | $108/128$ | $67/78$ |

**Constraintdrift der geführten und der naiven Strömung.** Maximaldrift über das jeweilige Ensemble bei $250$ Schritten:

| Methode | Ende $\tau$ | $\bar{\mathcal H}$ | $C_s$ | $C_c$ | $D$ |
|---|---|---|---|---|---|
| naiv, alle | $0{,}05$ | $1{,}01\cdot10^{-9}$ | $0$ exakt | $6{,}80\cdot10^{-2}$ | $5{,}11$ |
| naiv, alle | $0{,}10$ | $2{,}31\cdot10^{-3}$ | $0$ exakt | $1{,}87\cdot10^{-1}$ | $8{,}25\cdot10^{1}$ |
| naiv, alle | $0{,}20$ | $3{,}24\cdot10^{-3}$ | $0$ exakt | $2{,}75\cdot10^{-1}$ | $4{,}56\cdot10^{1}$ |
| naiv, alle | $0{,}40$ | $6{,}44\cdot10^{-4}$ | $0$ exakt | $4{,}58\cdot10^{-1}$ | $2{,}83\cdot10^{1}$ |
| Dirac, wohlbedingt | $0{,}05$ | $1{,}30\cdot10^{-1}$ | $0$ exakt | $4{,}95\cdot10^{-8}$ | $7{,}28\cdot10^{-1}$ |
| Dirac, wohlbedingt | $0{,}10$ | $1{,}05\cdot10^{1}$ | $0$ exakt | $1{,}62\cdot10^{-5}$ | $1{,}61\cdot10^{1}$ |
| Dirac, wohlbedingt | $0{,}20$ | $1{,}19\cdot10^{3}$ | $0$ exakt | $7{,}40\cdot10^{-3}$ | $1{,}44\cdot10^{3}$ |
| Dirac, wohlbedingt | $0{,}40$ | $4{,}19\cdot10^{24}$ | $0$ exakt | $1{,}76\cdot10^{17}$ | $3{,}64\cdot10^{24}$ |

Die Dirac-Zeilen sind für das Ensemble aller Punkte identisch, die Extremdrifte dieser Strömung liegen also im wohlbedingten Teilensemble. Bei der naiven Strömung ist es umgekehrt: ihr $D$-Drift fällt im wohlbedingten Ensemble für $\tau=0{,}20$ auf $5{,}74$ statt $4{,}56\cdot10^{1}$ (bei $\tau=0{,}40$ sind es $5{,}74$ gegenüber $2{,}83\cdot10^{1}$, dort jeweils mit $1000$ Schritten gemessen); der Maximaldrift gehört dort zu einem schlecht konditionierten Punkt. Drei Befunde sind sauber und wiederspruchsfrei. Erstens bleibt $C_s$ in allen Läufen beider Methoden exakt erhalten. Zweitens hält die Dirac-Strömung $C_c$ für kleines $\tau$ auf der Fläche, bei $\tau=0{,}05$ bis $4{,}95\cdot10^{-8}$ ($250$ Schritte) beziehungsweise $1{,}02\cdot10^{-6}$ ($1000$ Schritte), die naive Strömung nicht ($6{,}80\cdot10^{-2}$). Drittens verlässt die geführte Strömung die Fläche in $\bar{\mathcal H}$ und $D$ bereits für $\tau=0{,}10$ in der Größenordnung $10^1$ und für $\tau\ge0{,}2$ stark anwachsend. Das ist kein Rundungsartefakt: die mpmath-Referenz mit $30$ Stellen zeigt an denselben Randdaten dieselben Drifts, und einzelne Punkte laufen in endlicher Zeit über jede Betragsschranke — bei Punkt $45$ überschreitet die naive Trajektorie bei $\tau=0{,}10$ nach $154$ von $250$ Schritten die Schwelle $10^3$, und zwei weitere Läufe dort enden mit Überlauf.

**Konvergenzordnung des Schemas.** Für jede Kombination aus Methode, Ensemble und $\tau$ wurde mit $250$, $500$ und $1000$ Schritten gerechnet; das Verhältnis aufeinanderfolgender Endpunktdifferenzen muss für ein Verfahren vierter Ordnung nahe $2^4=16$ liegen. Gemessen:

| Ende $\tau$ | naiv, alle | naiv, wohlbedingt | Dirac, alle | Dirac, wohlbedingt |
|---|---|---|---|---|
| $0{,}05$ | $18{,}14$ | $18{,}14$ | $18{,}02$ | $16{,}30$ |
| $0{,}10$ | $16{,}32$ | $16{,}32$ | $18{,}75$ | $17{,}22$ |
| $0{,}20$ | $16{,}35$ | $16{,}55$ | $22{,}14$ | $1{,}01$ |
| $0{,}40$ | $16{,}52$ | $16{,}52$ | $2928{,}10$ | $19{,}81$ |

Die naive Strömung hat in allen vier Fällen das erwartete Verhältnis, die Diskretisierung ist also in Ordnung. Bei der Dirac-Strömung liegt es für $\tau\le0{,}1$ im Bereich $16$ bis $19$, darüber verlieren die Zahlen ihre Bedeutung, weil bei $10$ bis $30$ Trajektorien Endpunktdifferenzen nichtendlich werden. Die Abweichungen großer $\tau$-Werte stammen nicht vom Verfahren, sondern von der Divergenz einzelner Bahnen.

**Hochpräzisionsreferenz und Fehlerbudget.** Die drei Randpunkte mit dem kleinsten und den beiden größten $|W|$ sind Punkt $96$ ($|W|=2{,}41\cdot10^{-2}$, Rang $128$ von $128$), Punkt $37$ ($|W|=2{,}4248$) und Punkt $45$ ($|W|=2{,}3639$); sie sind in einer unabhängigen Nachrechnung genauso bestimmt worden. Für sie wurde die Flussrechnung mit $30$ Stellen und derselben Schrittweite $h=\tau/250$ wiederholt:

| Punkt | $|W|$ | Ende $\tau$ | naiv: $C_c$ / $D$ | Dirac: $C_c$ / $D$ |
|---|---|---|---|---|
| $96$ | $2{,}41\cdot10^{-2}$ | $0{,}10$ | $-7{,}48\cdot10^{-2}$ / $3{,}36$ | $6{,}02\cdot10^{-3}$ / $7{,}93\cdot10^{-2}$ |
| $96$ | | $0{,}20$ | nicht ausgeführt (Stufenzeitlimit) | $1{,}12\cdot10^{-2}$ / $2{,}56\cdot10^{1}$ |
| $37$ | $2{,}4248$ | $0{,}10$ | $-5{,}58\cdot10^{-2}$ / $8{,}04\cdot10^{-1}$ | $-1{,}14$ / $3{,}41$ |
| $37$ | | $0{,}20$ | $-8{,}37\cdot10^{-2}$ / $7{,}52\cdot10^{-1}$ | $-1{,}40$ / $9{,}61\cdot10^{-1}$ |
| $45$ | $2{,}3639$ | $0{,}10$ | divergiert (Schwelle $10^3$ bei Schritt $154$) | $2{,}97$ / $7{,}49$ |
| $45$ | | $0{,}20$ | divergiert (Überlauf) | divergiert (Überlauf) |

Der Vergleich von float64 und $30$ Stellen bei identischer Schrittweite liefert am Ende absolute Endpunktdifferenzen in der Größenordnung $1$: Punkt $37$ Dirac $1{,}87$ ($\tau=0{,}10$) und $4{,}39$ ($\tau=0{,}20$), naiv $1{,}22$ und $3{,}05$; Punkt $45$ Dirac $1{,}95$ ($\tau=0{,}10$), die übrigen Läufe dort divergieren in der Hochpräzisionsrechnung; Punkt $96$ Dirac $0{,}459$ und $1{,}43$, naiv $0{,}547$ ($\tau=0{,}10$) und divergierend ($\tau=0{,}20$). Relative Differenzen $0{,}30$ bis $1{,}22$ bedeuten: auf der Fläche ist durch Verfeinerung keine Bahngenauigkeit zu gewinnen. Bemerkenswert ist die Reihenfolge der Drifts. Für $C_c$ und $D$ zeigt der Punkt mit dem kleinsten $|W|=2{,}41\cdot10^{-2}$ die kleinsten Drifts ($6{,}02\cdot10^{-3}$ und $7{,}93\cdot10^{-2}$ bei $\tau=0{,}10$), während die Punkte mit $|W|\approx2{,}4$ weit stärker abweichen und einer von ihnen überläuft; bei $\bar{\mathcal H}$ ist die Reihenfolge nicht so geordnet (bei $\tau=0{,}10$ Dirac $0{,}318$ am Punkt $37$ gegenüber $1{,}19$ am Punkt $96$). Die Verstärkung über $u\sim W^{-1}$ aus der Fortpflanzungsgleichung ist damit als Kandidat für den zugrunde liegenden Mechanismus benannt, durch diese Läufe aber nicht belegt; die Driftverteilung ist in $|W|$ nicht monoton.

**Ergebnis.** Die in §19.7 offen gelassenen Punkte sind für das getrimmte gemeinsame System bis auf die $\alpha_c$-Vollständigkeit erledigt: Die Poisson-Matrix der Constraintmenge hat auf der vollen Fläche Rang $2$ mit nicht verschwindendem Block $\{C_c,D\}=W$; ihre Pfaffsche Determinante ist exakt $0$; die Dirac-Klammer ist in obiger geschlossener Form hingeschrieben und erfüllt ihre definierende Eigenschaft bis $4{,}44\cdot10^{-16}$; die reduzierte Phasenraumdimension ist $10-2\cdot2-2=4$, also zwei kanonische Paare. Neu gegenüber der ersten Zählung ist die Präzisierung der Klassifikation: $\bar{\mathcal H}$ selbst ist nicht erstklassig, erstklassig ist $\mathcal H_{\rm T}=\bar{\mathcal H}+u\,C_c$ mit $u=\{D,\bar{\mathcal H}\}/W$, und zwar exakt in dem Sinne, dass alle vier Klammern mit den Constraints im Constraintideal liegen. Der zugehörige geführte Fluss ist für kleines $\tau$ ($\tau\le0{,}05$) auf dem wohlbedingten Teilensemble quantitativ brauchbar, darüber hinaus nicht.

**Grenzen und offene Rechnungsschritte.** Erstens ist die Klassifikation an die geschlossene $\alpha_c=0$-Scheibe gebunden; über das vollständige Funktional ist schon $\{C_s,D\}\ne0$ (§19.7 Punkt (iv)), und eine $\alpha_c$-vollständige geschlossene Form von $\bar{\mathcal H}$ liegt nicht vor, bleibt also offen. Zweitens ist keine Eichfixierung der erstklassigen Menge $\{C_s,\mathcal H_{\rm T}\}$ getroffen; ohne sie ist die Parametrisierung der reduzierten Fläche und damit eine Führung im reduzierten Raum nicht eindeutig. Drittens bleibt ein die Constraintfläche tatsächlich projizierendes Integrationsverfahren (implizit oder mit Projektionsschritt nach jedem Schritt) offen; die hier verwendete explizite Runge-Kutta-Integration liefert die Konvergenzordnung des Schemas, aber keine Bahngenauigkeit auf der Fläche. Viertens sind Polarzerlegung, gemeinsames Quantenpotential, Strom und Führung auf dem reduzierten Raum der zwei kanonischen Paare noch nicht durchgeführt; die in §19.2 und §19.7 für den unreduzierten getrimmten Raum angegebenen Größen sind nicht mit der reduzierten Rechnung identisch zu setzen.

**Rechenhilfen und Laufzeiten.** Die Rechnung liegt in `work/phase4_model.py` (Modell, symbolische Klammern, vektorisierte Auswertung, mpmath-Route), `work/phase4_dirac_v5.py` (Stufen A bis I mit Fortschrittsausgabe, Zwischenstand `work/phase4_dirac_v5.json` nach jeder Stufe), `work/phase4_dirac_klaerung.py` (unabhängige Nachrechnung von Randdaten, Rängen und den Klammern von $\mathcal H_{\rm T}$) sowie `work/phase4_dirac_klassifikation_lean.lean` (formalisierte Algebra). Die gegenprüften Tabellenwerte werden von `work/pruefe_198.py` bis `work/pruefe_198d.py` unmittelbar aus dem Ergebnispapier erzeugt, damit Manuskript und Rechnung nicht auseinanderlaufen. Die unabhängigen Fluss- und Referenzläufe sind auf mehrere Rechnerprozesse verteilt ($6$ für das Ensemble, $4$ für die Hochpräzision). Die Stufenzeiten sind $0{,}41\,\mathrm{s}$ für die symbolischen Klammern, $8{,}74\,\mathrm{s}$ bis einschließlich Multiplikatorstufe, $59{,}7\,\mathrm{s}$ für die $48$ Ensembleflüsse und $312\,\mathrm{s}$ Gesamlaufzeit inklusive der Hochpräzisionsstufe und ihrer Zeitgrenzen; die Unabhängigkeits- und Rangnachrechnung braucht $5{,}8\,\mathrm{s}$. Eine frühere, vollsymbolisch mit `simplify` aufgebaute und seriell integrierte Fassung brauchte mehr als elf Minuten ohne jede Ausgabe auf einem Kern und wurde abgebrochen; die jetzige Fassung rechnet auf demselben Modell und mit denselben Zahlen.

"""


def main() -> None:
    pfad = Path("Publikation.md")
    text = pfad.read_text(encoding="utf-8", newline="")
    if MARKE in text:
        print("Sektion 19.8 ist bereits vorhanden, keine Aenderung.")
        return
    if KAPITEL not in text:
        raise SystemExit("Anker fuer Kapitel 20 nicht gefunden.")
    neuer = ABSCHNITT
    if "\r\n" in text:
        neuer = neuer.replace("\n", "\r\n")
    text = text.replace(KAPITEL, neuer + KAPITEL, 1)
    pfad.write_text(text, encoding="utf-8", newline="")
    print("Sektion 19.8 eingefuegt, neue Groesse", pfad.stat().st_size)


if __name__ == "__main__":
    main()
