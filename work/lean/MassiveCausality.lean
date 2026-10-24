import Mathlib

/-! Pointwise PG-current geometry. Strict time-likeness requires both chiral densities. -/
namespace MassiveCausality
noncomputable section

def q (A B : ℝ) : ℝ := A+B
def j (v A B : ℝ) : ℝ := (1-v)*A+(-1-v)*B
def velocity (v A B : ℝ) : ℝ := j v A B/q A B

theorem current_metric_identity (v A B : ℝ) :
    (q A B)^2-(j v A B+v*q A B)^2=4*A*B := by
  unfold q j
  ring

theorem causal_velocity_bounds (v A B : ℝ) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hq : 0 < q A B) : -1-v ≤ velocity v A B ∧ velocity v A B ≤ 1-v := by
  unfold velocity
  constructor
  · apply (le_div_iff₀ hq).2
    dsimp [q,j]
    nlinarith
  · apply (div_le_iff₀ hq).2
    dsimp [q,j]
    nlinarith

theorem normalized_metric_identity (v A B : ℝ) (hq : q A B ≠ 0) :
    1-(velocity v A B+v)^2=4*A*B/(q A B)^2 := by
  unfold velocity
  field_simp
  unfold q j
  ring

theorem both_components_strictly_timelike (v A B : ℝ) (hA : 0 < A) (hB : 0 < B) :
    0 < 1-(velocity v A B+v)^2 := by
  have hq : 0 < q A B := by dsimp [q]; positivity
  rw [normalized_metric_identity v A B (ne_of_gt hq)]
  positivity

theorem timelike_iff_both_components (v A B : ℝ) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hq : 0 < q A B) :
    0 < 1-(velocity v A B+v)^2 ↔ 0 < A ∧ 0 < B := by
  constructor
  · intro ht
    rw [normalized_metric_identity v A B (ne_of_gt hq)] at ht
    have hn : 0 < 4*A*B := (div_pos_iff.mp ht).resolve_right (by
      intro hh
      exact (not_lt_of_ge (sq_nonneg (q A B))) hh.2)
      |>.1
    constructor <;> nlinarith
  · rintro ⟨ha,hb⟩
    exact both_components_strictly_timelike v A B ha hb

theorem upper_only_current_is_null (v A : ℝ) (hA : A ≠ 0) :
    velocity v A 0=1-v ∧ 1-(velocity v A 0+v)^2=0 := by
  have he : velocity v A 0=1-v := by simp [velocity,j,q,hA]
  rw [he]
  constructor <;> ring

theorem lower_only_current_is_null (v B : ℝ) (hB : B ≠ 0) :
    velocity v 0 B= -1-v ∧ 1-(velocity v 0 B+v)^2=0 := by
  have he : velocity v 0 B= -1-v := by simp [velocity,j,q,hB]
  rw [he]
  constructor <;> ring

theorem strictly_inward_inside (v A B : ℝ) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hq : 0 < q A B) (hv : 1 < v) : velocity v A B < 0 := by
  have hh := (causal_velocity_bounds v A B hA hB hq).2
  linarith

theorem horizon_velocity (A B : ℝ) : velocity 1 A B= -2*B/q A B := by
  unfold velocity j
  ring

theorem horizon_strictly_inward_with_lower_component (A B : ℝ)
    (hA : 0 ≤ A) (hB : 0 < B) : velocity 1 A B < 0 := by
  rw [horizon_velocity]
  apply div_neg_of_neg_of_pos
  · linarith
  · dsimp [q]; positivity

theorem real_spinor_density_positive (a b c e : ℝ)
    (hn : a ≠ 0 ∨ b ≠ 0 ∨ c ≠ 0 ∨ e ≠ 0) :
    0 < q (a^2+b^2) (c^2+e^2) := by
  unfold q
  rcases hn with ha | hb | hc | he
  · nlinarith [sq_pos_of_ne_zero ha,sq_nonneg b,sq_nonneg c,sq_nonneg e]
  · nlinarith [sq_pos_of_ne_zero hb,sq_nonneg a,sq_nonneg c,sq_nonneg e]
  · nlinarith [sq_pos_of_ne_zero hc,sq_nonneg a,sq_nonneg b,sq_nonneg e]
  · nlinarith [sq_pos_of_ne_zero he,sq_nonneg a,sq_nonneg b,sq_nonneg c]

theorem schwarzschild_interior_shift (M r : ℝ) (hM : 0 < M) (hr : 0 < r)
    (hinside : r < 2*M) : 1 < Real.sqrt (2*M/r) := by
  have hp : 0 < 2*M/r := by positivity
  have hratio : 1 < 2*M/r := (lt_div_iff₀ hr).2 (by simpa using hinside)
  have hs := Real.sq_sqrt (le_of_lt hp)
  have hn := Real.sqrt_nonneg (2*M/r)
  nlinarith

theorem schwarzschild_interior_guidance (M r A B : ℝ)
    (hM : 0 < M) (hr : 0 < r) (hinside : r < 2*M)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hq : 0 < q A B) :
    velocity (Real.sqrt (2*M/r)) A B < 0 :=
  strictly_inward_inside _ A B hA hB hq (schwarzschild_interior_shift M r hM hr hinside)

end
end MassiveCausality
