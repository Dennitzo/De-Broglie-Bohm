import Mathlib

/-! Infinite frequency and angular convergence of an explicit joint majorant.
The spectral comparison that bounds a particular Green kernel is an input. -/
namespace FullRSETUV
noncomputable section
open MeasureTheory

def majorant (C a t : ℝ) (l : ℕ) (k : ℝ) : ℝ :=
  C*(2*(l : ℝ)+1)*t^l*Real.exp (-a*k)

theorem angular_weighted_geometric_summable (t : ℝ) (ht : |t| < 1) :
    Summable (fun l : ℕ => (2*(l : ℝ)+1)*t^l) := by
  have hp := summable_pow_mul_geometric_of_norm_lt_one 1 (r := t) (by simpa using ht)
  have hg := summable_geometric_of_norm_lt_one (x := t) (by simpa using ht)
  convert (hp.mul_left 2).add hg using 1
  funext l
  simp only [pow_one]
  ring

theorem frequency_majorant_integrable (C a t : ℝ) (l : ℕ) (ha : 0 < a) :
    IntegrableOn (majorant C a t l) (Set.Ioi 0) := by
  exact (integrableOn_exp_mul_Ioi (by linarith : -a < 0) 0).const_mul (C*(2*(l : ℝ)+1)*t^l)

theorem frequency_majorant_integral (C a t : ℝ) (l : ℕ) (ha : 0 < a) :
    (∫ k in Set.Ioi (0 : ℝ), majorant C a t l k) = C*(2*(l : ℝ)+1)*t^l/a := by
  unfold majorant
  rw [integral_const_mul, integral_exp_mul_Ioi (by linarith : -a < 0)]
  simp
  ring

theorem joint_majorant_summable (C a t : ℝ) (ha : 0 < a) (ht : |t| < 1) :
    Summable (fun l : ℕ => ∫ k in Set.Ioi (0 : ℝ), majorant C a t l k) := by
  simp_rw [frequency_majorant_integral C a t _ ha]
  convert (angular_weighted_geometric_summable t ht).mul_left (C/a) using 1
  funext l
  ring

theorem exponentially_bounded_mode_integrable (F : ℕ → ℝ → ℝ) (C a t : ℝ)
    (ha : 0 < a) (l : ℕ)
    (hm : AEStronglyMeasurable (F l) (volume.restrict (Set.Ioi 0)))
    (hb : ∀ k ∈ Set.Ioi (0 : ℝ), |F l k| ≤ majorant C a t l k) :
    IntegrableOn (F l) (Set.Ioi 0) := by
  apply (frequency_majorant_integrable C a t l ha).mono' hm
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with k hk
  simpa [Real.norm_eq_abs] using hb k hk

theorem exponentially_bounded_nonnegative_joint_sum (F : ℕ → ℝ → ℝ) (C a t : ℝ)
    (ha : 0 < a) (ht : |t| < 1)
    (hm : ∀ l, AEStronglyMeasurable (F l) (volume.restrict (Set.Ioi 0)))
    (hb : ∀ l k, k ∈ Set.Ioi (0 : ℝ) → 0 ≤ F l k ∧ F l k ≤ majorant C a t l k) :
    Summable (fun l : ℕ => ∫ k in Set.Ioi (0 : ℝ), F l k) := by
  apply (joint_majorant_summable C a t ha ht).of_nonneg_of_le
  · intro l
    apply integral_nonneg_of_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with k hk
    exact (hb l k hk).1
  · intro l
    have hF := exponentially_bounded_mode_integrable F C a t ha l (hm l) (by
      intro k hk
      rw [abs_of_nonneg (hb l k hk).1]
      exact (hb l k hk).2)
    apply integral_mono_ae hF (frequency_majorant_integrable C a t l ha)
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with k hk
    exact (hb l k hk).2

end
end FullRSETUV
