import Mathlib

/-!
# The deterministic CHSH bound and finite probability ensembles

For each hidden state the same four local response values A0, A1, B0, B1
occur in all four products. A single nonnegative normalized weight function
is used for all setting pairs. These are mathematical hypotheses, not
conclusions about any particular physical experiment.
-/

open scoped BigOperators

namespace BohmBell

def chsh (a0 a1 b0 b1 : ℝ) : ℝ :=
  a0 * b0 + a0 * b1 + a1 * b0 - a1 * b1

theorem chsh_factorization (a0 a1 b0 b1 : ℝ) :
    chsh a0 a1 b0 b1 = a0 * (b0 + b1) + a1 * (b0 - b1) := by
  unfold chsh
  ring

theorem deterministic_chsh_abs_eq_two (a0 a1 b0 b1 : ℝ)
    (ha0 : a0 = 1 ∨ a0 = -1) (ha1 : a1 = 1 ∨ a1 = -1)
    (hb0 : b0 = 1 ∨ b0 = -1) (hb1 : b1 = 1 ∨ b1 = -1) :
    |chsh a0 a1 b0 b1| = 2 := by
  rcases ha0 with rfl | rfl <;> rcases ha1 with rfl | rfl <;>
    rcases hb0 with rfl | rfl <;> rcases hb1 with rfl | rfl <;>
    norm_num [chsh]

theorem deterministic_chsh_bound (a0 a1 b0 b1 : ℝ)
    (ha0 : a0 = 1 ∨ a0 = -1) (ha1 : a1 = 1 ∨ a1 = -1)
    (hb0 : b0 = 1 ∨ b0 = -1) (hb1 : b1 = 1 ∨ b1 = -1) :
    |chsh a0 a1 b0 b1| ≤ 2 := by
  exact le_of_eq (deterministic_chsh_abs_eq_two a0 a1 b0 b1 ha0 ha1 hb0 hb1)

theorem weighted_chsh_bound {ι : Type*} (s : Finset ι)
    (w a0 a1 b0 b1 : ι → ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (hnorm : ∑ i ∈ s, w i = 1)
    (ha0 : ∀ i ∈ s, a0 i = 1 ∨ a0 i = -1)
    (ha1 : ∀ i ∈ s, a1 i = 1 ∨ a1 i = -1)
    (hb0 : ∀ i ∈ s, b0 i = 1 ∨ b0 i = -1)
    (hb1 : ∀ i ∈ s, b1 i = 1 ∨ b1 i = -1) :
    |∑ i ∈ s, w i * chsh (a0 i) (a1 i) (b0 i) (b1 i)| ≤ 2 := by
  have hpoint (i : ι) (hi : i ∈ s) :
      -2 ≤ chsh (a0 i) (a1 i) (b0 i) (b1 i) ∧
      chsh (a0 i) (a1 i) (b0 i) (b1 i) ≤ 2 :=
    abs_le.mp (deterministic_chsh_bound _ _ _ _
      (ha0 i hi) (ha1 i hi) (hb0 i hi) (hb1 i hi))
  apply abs_le.mpr
  constructor
  · calc
      -2 = ∑ i ∈ s, w i * (-2) := by rw [← Finset.sum_mul, hnorm]; norm_num
      _ ≤ ∑ i ∈ s, w i * chsh (a0 i) (a1 i) (b0 i) (b1 i) := by
        exact Finset.sum_le_sum fun i hi =>
          mul_le_mul_of_nonneg_left (hpoint i hi).1 (hw i hi)
  · calc
      (∑ i ∈ s, w i * chsh (a0 i) (a1 i) (b0 i) (b1 i))
          ≤ ∑ i ∈ s, w i * 2 := by
        exact Finset.sum_le_sum fun i hi =>
          mul_le_mul_of_nonneg_left (hpoint i hi).2 (hw i hi)
      _ = 2 := by rw [← Finset.sum_mul, hnorm]; norm_num

def correlation {ι : Type*} (s : Finset ι) (w a b : ι → ℝ) : ℝ :=
  ∑ i ∈ s, w i * (a i * b i)

theorem chsh_correlation_bound {ι : Type*} (s : Finset ι)
    (w a0 a1 b0 b1 : ι → ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (hnorm : ∑ i ∈ s, w i = 1)
    (ha0 : ∀ i ∈ s, a0 i = 1 ∨ a0 i = -1)
    (ha1 : ∀ i ∈ s, a1 i = 1 ∨ a1 i = -1)
    (hb0 : ∀ i ∈ s, b0 i = 1 ∨ b0 i = -1)
    (hb1 : ∀ i ∈ s, b1 i = 1 ∨ b1 i = -1) :
    |correlation s w a0 b0 + correlation s w a0 b1 +
      correlation s w a1 b0 - correlation s w a1 b1| ≤ 2 := by
  have h := weighted_chsh_bound s w a0 a1 b0 b1 hw hnorm ha0 ha1 hb0 hb1
  simpa only [chsh, correlation, mul_sub, mul_add,
    Finset.sum_add_distrib, Finset.sum_sub_distrib] using h

end BohmBell

#print axioms BohmBell.chsh_factorization
#print axioms BohmBell.deterministic_chsh_abs_eq_two
#print axioms BohmBell.deterministic_chsh_bound
#print axioms BohmBell.weighted_chsh_bound
#print axioms BohmBell.chsh_correlation_bound
