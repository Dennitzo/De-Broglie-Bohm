import Mathlib.MeasureTheory.Integral
import Mathlib.MeasureTheory.Measure.Prod

noncomputable section
open scoped ENNReal
open MeasureTheory

namespace GeneralNormalization

/-- Konditionierte Dichte: Integral des Betragsquadrats über das Feld bei
festgehaltener Geometrievariable. -/
def condDensity (Psi : ℝ × ℝ → ℂ) (h : ℝ) : ℝ :=
  ∫ f, Complex.normSq (Psi (h, f))

/-- Die Zustandsklasse wird durch die Messbarkeit und die Integrierbarkeit
des Betragsquadrats der Wellenfunktion beschränkt. -/
structure ZustandsKlasse where
  Psi : ℝ × ℝ → ℂ
  meas : Measurable (fun p : ℝ × ℝ => Psi (p.1, p.2))
  integrable : Integrable (fun p : ℝ × ℝ => Complex.normSq (Psi p))

/-- Für jede Zustandsklasse, die normiert ist und die Fubini-Vertauschung erfüllt,
ist die konditionierte Dichte normalisiert. Die Voraussetzungen sind die
Messbarkeit und die Integrierbarkeit der Zustandsklasse sowie die Normierung
des Gesamtzustands und die Vertauschung der Integrationsreihenfolge. -/
theorem general_conditional_density_normalized
    (Z : ZustandsKlasse)
    (hSwap : ∫ h, condDensity Z.Psi h = ∫ p : ℝ × ℝ, Complex.normSq (Z.Psi p))
    (hNorm : ∫ p : ℝ × ℝ, Complex.normSq (Z.Psi p) = 1) :
    ∫ h, condDensity Z.Psi h = 1 := by
  have h_meas_use : Measurable (fun p : ℝ × ℝ => Z.Psi (p.1, p.2)) := Z.meas
  have h_int_use : Integrable (fun p : ℝ × ℝ => Complex.normSq (Z.Psi p)) := Z.integrable
  exact hSwap.trans hNorm

end GeneralNormalization
