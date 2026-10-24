import Mathlib

/-!
# Lokale Identität für die feinkörnige relative Entropie

Der Integral- und Randsatz aus Abschnitt 5.2 der Publikation beruht auf einer
lokalen Identität: Wenn die Dichte q die Kontinuitätsgleichung erfüllt und die
Funktion f die Materialableitung D_t f = 0 besitzt, dann ist die lokale Größe
q * Phi(f) für jede differenzierbare Funktion Phi ebenfalls erhalten. Dieser
Satz formalisiert die algebraische Kernaussage an einem Punkt. Voraussetzungen
und Bedeutung werden getrennt genannt; Existenz, Regularität und Randterme des
Integralsatzes sind hier nicht enthalten.
-/

namespace BohmTransport

/-- Lokale Erhaltung von q * Phi(f) entlang des gemeinsamen Geschwindigkeitsfelds.

Wenn q die Kontinuitätsgleichung dq + v * dxq = 0 erfüllt und f die Gleichung
D_t f = 0, also df + v * dxf = 0, dann ist die lokale Ableitung
d(q * Phi(f)) + v * dx(q * Phi(f)) = 0.

Hier sind Phi und Phi' die Werte der differenzierbaren Funktion beziehungsweise
ihrer Ableitung an der Stelle f. Die Formel folgt aus der Produktregel und der
Kettenregel; sie beweist keine Integrabilität und keinen Randterm. -/
theorem entropy_density_transport (q v dq df dxq dxf Phi Phi' : ℝ)
    (hq : dq + v * dxq = 0) (hf : df + v * dxf = 0) :
    (dq * Phi + q * Phi' * df) + v * (dxq * Phi + q * Phi' * dxf) = 0 := by
  calc
    (dq * Phi + q * Phi' * df) + v * (dxq * Phi + q * Phi' * dxf)
        = Phi * (dq + v * dxq) + q * Phi' * (df + v * dxf) := by
          ring
    _ = Phi * 0 + q * Phi' * 0 := by
      rw [hq, hf]
    _ = 0 := by
      simp

/-- Spezialfall für die feinkörnige relative Entropie Phi(f) = f * log f.

Setzt man Phi(f) = f * log f und Phi'(f) = log f + 1, so folgt aus der
vorstehenden Identität die lokale Erhaltung von q * f * log f. Der Beweis ist
rein algebraisch und setzt keine Differenzierbarkeit von log im reellen Sinn
voraus; er übernimmt die angegebenen Werte von Phi und Phi'. -/
theorem entropy_density_relative (q f v dq df dxq dxf : ℝ)
    (hq : dq + v * dxq = 0) (hf : df + v * dxf = 0) :
    (dq * (f * Real.log f) + q * (Real.log f + 1) * df)
        + v * (dxq * (f * Real.log f) + q * (Real.log f + 1) * dxf) = 0 := by
  exact entropy_density_transport q v dq df dxq dxf (f * Real.log f) (Real.log f + 1) hq hf

end BohmTransport
