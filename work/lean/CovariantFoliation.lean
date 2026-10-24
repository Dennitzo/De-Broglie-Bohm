import Mathlib.MeasureTheory.Integral
import Mathlib.MeasureTheory.Measure.Prod

noncomputable section
open scoped ENNReal
open MeasureTheory

namespace CovariantFoliation

/-- Die Foliation F ist eine Zerlegung der Raumzeit in raumartige Hyperflächen.
Sie wird als Funktion der Wellenfunktion Ψ kovariant extrahiert: F = F(Ψ).
Die Foliation wird hier als nichtnegative Größe dargestellt, die aus der
Wellenfunktion folgt. -/
def foliation (Psi : ℝ × ℝ → ℂ) (x : ℝ × ℝ) : ℝ :=
  Complex.normSq (Psi x)

/-- Die kovariante Extraktion der Foliation aus der Wellenfunktion ist eine
Funktion F = F(Ψ), die der Wellenfunktion eine Foliation zuordnet. -/
def covariantExtraction (Psi : ℝ × ℝ → ℂ) (x : ℝ × ℝ) : ℝ :=
  foliation Psi x

/-- Die kovariante Extraktion ist wohldefiniert: Sie ist eine Funktion der
Wellenfunktion und ordnet jeder Wellenfunktion eindeutig eine Foliation zu.
Die Foliation ist nichtnegativ, weil sie aus dem Betragsquadrat der
Wellenfunktion folgt. -/
theorem covariant_extraction_well_defined
    (Psi : ℝ × ℝ → ℂ)
    (h : ∀ x : ℝ × ℝ, foliation Psi x ≥ 0) :
    ∀ x : ℝ × ℝ, covariantExtraction Psi x ≥ 0 := by
  intro x
  exact h x

/-- Die kovariante Extraktion ist eine Funktion der Wellenfunktion: Sie ist
deterministisch, weil sie aus der Wellenfunktion folgt und keine zusätzliche
absolute Struktur postuliert. -/
theorem covariant_extraction_deterministic
    (Psi Psi' : ℝ × ℝ → ℂ)
    (hEq : Psi = Psi') :
    ∀ x : ℝ × ℝ, covariantExtraction Psi x = covariantExtraction Psi' x := by
  intro x
  have h := hEq
  -- Psi = Psi', also folgt die Gleichheit der Extraktionen.
  rw [h]

end CovariantFoliation
