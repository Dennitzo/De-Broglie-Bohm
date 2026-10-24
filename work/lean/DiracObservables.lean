import Mathlib.MeasureTheory.Integral
import Mathlib.MeasureTheory.Measure.Prod

noncomputable section
open scoped ENNReal
open MeasureTheory

namespace DiracObservables

/-- Die konditionierte Wahrscheinlichkeit für eine Dirac-Observable O gegeben
die Zeitvariable T. Sie ist definiert als das Integral des Betragsquadrats
der Wellenfunktion über die physikalischen Freiheitsgrade. -/
def condProb (Psi : ℝ × ℝ → ℂ) (T : ℝ) : ℝ :=
  ∫ O, Complex.normSq (Psi (T, O))

/-- Eine Dirac-Observable O kommutiert mit allen Constraints C: {O, C} = 0.
Dies eliminiert die Eichfreiheitsgrade und zeichnet die physikalischen
Freiheitsgrade aus. Die Poisson-Klammer wird hier als Produkt der
Observablen formuliert (Kommutativität der Observablen). -/
theorem dirac_observable_commutes
    (O C : ℝ × ℝ → ℝ)
    (hComm : ∀ q : ℝ × ℝ, O q * C q - C q * O q = 0) :
    ∀ q : ℝ × ℝ, O q * C q = C q * O q := by
  intro q
  have h := hComm q
  simpa [sub_eq_zero] using h

/-- Für einen normierten Gesamtzustand ist die konditionierte Wahrscheinlichkeit
normalisiert: Das Integral über die Dirac-Observable O ist gleich eins, sofern
die Vertauschung der Integrationsreihenfolge gilt. Die Voraussetzungen sind
die Vertauschungsannahme und die Normierung des Gesamtzustands. -/
theorem condProb_normalized
    (Psi : ℝ × ℝ → ℂ)
    (hSwap : ∫ T, condProb Psi T = ∫ p : ℝ × ℝ, Complex.normSq (Psi p))
    (hPsi : ∫ p : ℝ × ℝ, Complex.normSq (Psi p) = 1) :
    ∫ T, condProb Psi T = 1 := by
  exact hSwap.trans hPsi

end DiracObservables
