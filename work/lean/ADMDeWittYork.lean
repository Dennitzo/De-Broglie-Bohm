import Mathlib

namespace ADMDeWittYork
noncomputable section
open Matrix

/-- Die spurumgekehrte Abbildung, wie sie in der ADM-Beziehung zwischen
kanonischem Impuls und Extrinischer Krümmung in einem orthonormalen
drei-dimensionalen Rahmen auftritt. -/
def traceReverseWith (c : ℝ) (A : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix (Fin 3) (Fin 3) ℝ :=
  A - c • (Matrix.trace A) • (1 : Matrix (Fin 3) (Fin 3) ℝ)

theorem trace_reverse_with_trace (c : ℝ) (A : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix.trace (traceReverseWith c A) = (1 - 3 * c) * Matrix.trace A := by
  unfold traceReverseWith
  rw [Matrix.trace_sub, Matrix.trace_smul]
  have h : ((Matrix.trace A : ℝ) • (1 : Matrix (Fin 3) (Fin 3) ℝ)).trace =
      Matrix.trace A * 3 := by
    rw [Matrix.trace_smul]
    have h1 : Matrix.trace (1 : Matrix (Fin 3) (Fin 3) ℝ) = (3 : ℝ) := by
      simp
    rw [h1]
    simp [smul_eq_mul]
  rw [h]
  simp only [smul_eq_mul]
  ring

theorem trace_reverse_with_invert {c : ℝ} (hc : 1 - 3 * c ≠ 0)
    (A : Matrix (Fin 3) (Fin 3) ℝ) :
    traceReverseWith c A + (c / (1 - 3 * c)) •
      (Matrix.trace (traceReverseWith c A)) • (1 : Matrix (Fin 3) (Fin 3) ℝ) = A := by
  ext i j
  rw [trace_reverse_with_trace]
  have hc' : (1 : ℝ) - c * 3 ≠ 0 := by
    simpa [mul_comm] using hc
  simp only [traceReverseWith, Matrix.add_apply, Matrix.sub_apply,
    Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
  split <;> field_simp [hc'] <;> ring

theorem adm_trace_reverse_inverse (A : Matrix (Fin 3) (Fin 3) ℝ) :
    traceReverseWith (1 / 2) A +
      ((1 / 2 : ℝ) / (1 - 3 * (1 / 2 : ℝ))) •
        Matrix.trace (traceReverseWith (1 / 2) A) •
      (1 : Matrix (Fin 3) (Fin 3) ℝ) = A :=
  trace_reverse_with_invert (by norm_num : 1 - 3 * (1 / 2 : ℝ) ≠ 0) A

/-- Der Cardano-Radius für eine reduzierte kubische Gleichung mit negativem
linearen Koeffizienten, wie er in der trigonometrischen Darstellung der
York-Masse auftritt. -/
def yorkRadius (p : ℝ) : ℝ := Real.sqrt (-p / 3)

/-- Die trigonometrische Cardano-Wurzel. -/
def yorkAngleRoot (p θ : ℝ) : ℝ := 2 * yorkRadius p * Real.cos θ

theorem york_radius_sq {p : ℝ} (hp : p < 0) : yorkRadius p ^ 2 = -p / 3 :=
  Real.sq_sqrt (by linarith)

theorem york_cubic_trig {p θ : ℝ} (hp : p < 0) :
    yorkAngleRoot p θ ^ 3 + p * yorkAngleRoot p θ =
      2 * yorkRadius p ^ 3 * Real.cos (3 * θ) := by
  have hr := york_radius_sq hp
  have h : p + 3 * yorkRadius p ^ 2 = 0 := by
    rw [hr]
    ring
  unfold yorkAngleRoot
  calc
    (2 * yorkRadius p * Real.cos θ) ^ 3 + p * (2 * yorkRadius p * Real.cos θ)
        = 2 * yorkRadius p ^ 3 *
            (4 * Real.cos θ ^ 3 - 3 * Real.cos θ) +
            2 * yorkRadius p * Real.cos θ * (p + 3 * yorkRadius p ^ 2) := by
          ring
    _ = 2 * yorkRadius p ^ 3 *
            (4 * Real.cos θ ^ 3 - 3 * Real.cos θ) := by
          rw [h]
          ring
    _ = 2 * yorkRadius p ^ 3 * Real.cos (3 * θ) := by
          rw [← Real.cos_three_mul θ]

theorem york_cubic_trig_add {p q θ : ℝ} (hp : p < 0) :
    yorkAngleRoot p θ ^ 3 + p * yorkAngleRoot p θ + q =
      2 * yorkRadius p ^ 3 * Real.cos (3 * θ) + q := by
  rw [york_cubic_trig hp]

end
end ADMDeWittYork
