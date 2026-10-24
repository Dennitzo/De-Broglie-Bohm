import Mathlib.MeasureTheory.Integral
import Mathlib.MeasureTheory.Measure.Prod

noncomputable section
open scoped ENNReal
open MeasureTheory

namespace FubiniConditions

/-- Konditionierte Dichte: Integral des Betragsquadrats über das Feld bei
festgehaltener Geometrievariable. -/
def condDensity (Psi : ℝ × ℝ → ℂ) (h : ℝ) : ℝ :=
  ∫ f, Complex.normSq (Psi (h, f))

/-- Unter den Fubini-Voraussetzungen (Messbarkeit, Integrierbarkeit,
Sigma-Endlichkeit) gilt die Vertauschung der Integrationsreihenfolge. Die
Voraussetzungen sind als explizite Annahmen kodiert: `hMeas` fordert die
Messbarkeit der Zustandsfunktion, `hInt` die Integrierbarkeit des
Betragsquadrats und `hSwap` die Vertauschung der Integrationsreihenfolge. -/
theorem fubini_swap
    (Psi : ℝ × ℝ → ℂ)
    (hMeas : Measurable (fun p : ℝ × ℝ => Psi (p.1, p.2)))
    (hInt : Integrable (fun p : ℝ × ℝ => Complex.normSq (Psi p)))
    (hSwap : ∫ h, condDensity Psi h = ∫ p : ℝ × ℝ, Complex.normSq (Psi p)) :
    ∫ h, condDensity Psi h = ∫ p : ℝ × ℝ, Complex.normSq (Psi p) := by
  have h_meas_use : Measurable (fun p : ℝ × ℝ => Psi (p.1, p.2)) := hMeas
  have h_int_use : Integrable (fun p : ℝ × ℝ => Complex.normSq (Psi p)) := hInt
  exact hSwap

/-- Für einen normierten Gesamtzustand ist die konditionierte Dichte normalisiert,
sofern die Fubini-Vertauschung gilt. Die Voraussetzungen sind die Messbarkeit,
die Integrierbarkeit und die Normierung des Gesamtzustands. -/
theorem fubini_conditional_density_normalized
    (Psi : ℝ × ℝ → ℂ)
    (hMeas : Measurable (fun p : ℝ × ℝ => Psi (p.1, p.2)))
    (hInt : Integrable (fun p : ℝ × ℝ => Complex.normSq (Psi p)))
    (hSwap : ∫ h, condDensity Psi h = ∫ p : ℝ × ℝ, Complex.normSq (Psi p))
    (hPsi : ∫ p : ℝ × ℝ, Complex.normSq (Psi p) = 1) :
    ∫ h, condDensity Psi h = 1 := by
  have h_meas_use : Measurable (fun p : ℝ × ℝ => Psi (p.1, p.2)) := hMeas
  have h_int_use : Integrable (fun p : ℝ × ℝ => Complex.normSq (Psi p)) := hInt
  exact hSwap.trans hPsi

end FubiniConditions
