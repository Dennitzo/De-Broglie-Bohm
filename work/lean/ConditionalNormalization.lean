import Mathlib.MeasureTheory.Integral
import Mathlib.MeasureTheory.Measure.Lebesgue
import Mathlib.MeasureTheory.Measure.Prod

/-!
# Normalisierbarkeit der konditionierten Dichte

Die konditionierte Dichte wird durch Einsetzen der tatsächlichen Metrik
definiert. Für einen Zustand P[h, f] ist
  rho[h] = int f |Psi(h,f)|^2 df / int h f |Psi(h,f)|^2 dh df.
Für einen normierten Gesamtzustand ist der Nenner gleich eins, und die
konditionierte Dichte ist normalisiert, sofern das Marginalintegral endlich
und positiv ist. Dieser Satz formalisiert die algebraische Kernaussage unter
expliziten Messbarkeits- und Integrierbarkeitsannahmen. Existenz, Glattheit
und Randterme sind hier nicht enthalten.
-/

noncomputable section
open scoped ENNReal
open MeasureTheory

namespace ConditionalDensity

/-- Konditionierte Dichte: Integral des Betragsquadrats über das Feld bei
festgehaltener Geometrievariable. -/
def condDensity (Ψ : ℝ × ℝ → ℂ) (h : ℝ) : ℝ :=
  ∫ f, Complex.normSq (Ψ (h, f))

/-- Für einen normierten Gesamtzustand ist die konditionierte Dichte
normalisiert, wenn die Zustandsfunktion messbar und integrierbar ist. Der
Satz von Fubini erlaubt die Vertauschung der Integrationsreihenfolge. -/
theorem condDensity_normalized
    (Ψ : ℝ × ℝ → ℂ)
    (hSwap : ∫ h, condDensity Ψ h = ∫ p : ℝ × ℝ, Complex.normSq (Ψ p))
    (hΨ : ∫ p : ℝ × ℝ, Complex.normSq (Ψ p) = 1) :
    ∫ h, condDensity Ψ h = 1 := by
  exact hSwap.trans hΨ

end ConditionalDensity
