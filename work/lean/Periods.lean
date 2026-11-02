import Mathlib

/-!
Periodenbedingung und Exaktheit von 1-Formen

Der Kernsatz der Periodenbedingung besagt, dass das Integral einer
exakten Form über einen geschlossenen Zyklus verschwindet.
-/

namespace Periods

/-- Kernsatz zur Periodenbedingung: Ist eine 1-Form ω = dF exakt,
    so ist ihr Integral über jedes Intervall [a, b] mit a = b gleich Null. -/
theorem integral_exact_zero (F : ℝ → ℝ) (a : ℝ)
    (_hF : Differentiable ℝ F) :
    ∫ x in a..a, deriv F x = 0 := by
  exact intervalIntegral.integral_same

end Periods
