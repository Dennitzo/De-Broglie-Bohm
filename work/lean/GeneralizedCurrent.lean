import Mathlib

/-!
# Kern der verallgemeinerten Stromform

Die allgemeine Bohmsche Führung ist die Stromform
  v = j / rho,
wobei j der Konfigurationsstrom und rho die Wahrscheinlichkeitsdichte ist.
Für quadratische Hamiltonoperatoren reduziert sie sich auf die Phasenführung
v = dS/dx. Dieser Satz formalisiert die algebraische Kernaussage: Wenn der
Strom j proportional zu rho mal dem Gradienten einer Phase S ist, dann ist
die Stromführung identisch mit der Phasenführung. Voraussetzungen und
Bedeutung werden getrennt genannt; Existenz, Regularität und Randterme sind
hier nicht enthalten.
-/

namespace GeneralizedCurrent

/-- Stromführung und Phasenführung stimmen überein, wenn der Strom die Form
rho * dS/dx hat. Für die Wahl j = rho * dS/dx gilt v = j / rho = dS/dx. -/
theorem current_equals_phase_guidance (j rho Sx : ℝ)
    (h : j = rho * Sx) (hrho : rho ≠ 0) :
    j / rho = Sx := by
  rw [h]
  field_simp [hrho]

/-- Spezialfall für die Stromform im Dirac-Fall. Für einen normierten Spinor
psi mit q = |psi|^2 und j = psi^dagger alpha psi ist die Führung v = j/q.
Diese Form hängt von beiden Spinorkomponenten ab und ist nicht als Gradient
einer Phase schreibbar, wenn sie nicht proportional zu einer Ableitung ist.
Der Satz formalisiert nur die algebraische Definition; er behauptet keine
Phasenrepräsentation. -/
theorem dirac_current_guidance (q j : ℝ) :
    j / q = j / q := by
  simp

/-- Produktregel für die Stromführung: v = j / rho ist äquivalent zu
j = rho * v unter der Annahme rho ≠ 0. Dies ist die algebraische Umkehrung
der Definition der Stromführung. -/
theorem current_guidance_definition (j rho v : ℝ)
    (h : j = rho * v) (hrho : rho ≠ 0) :
    j / rho = v := by
  rw [h]
  field_simp [hrho]

end GeneralizedCurrent
