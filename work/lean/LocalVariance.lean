import Mathlib

/-! Algebra of Gaussian moments and divergence of the cubic-cutoff variance.
The connection of the moment values to a Gaussian measure is analytic input;
the polynomial lower bound and divergence are checked here. -/
namespace LocalVariance
noncomputable section
open Filter
open scoped Topology BigOperators

def cutoffVariance (m c n : ℝ) : ℝ :=
  (2*n+1)^3*(m^2+c^2*n*(n+1))/8

theorem centered_second_moment_shift (first second shift : ℝ) :
    (second-2*shift*first+shift^2)-(first-shift)^2=second-first^2 := by ring

theorem gaussian_mode_potential_mean (w : ℝ) (hw : w ≠ 0) :
    w^2/2*(1/(2*w))=w/4 := by field_simp; ring

theorem gaussian_mode_potential_variance (w : ℝ) (hw : w ≠ 0) :
    w^4/4*(3/(4*w^2)-(1/(2*w))^2)=w^2/8 := by field_simp; ring

theorem positive_integer_square_sum (n : ℕ) :
    6*(∑ k ∈ Finset.range (n+1), (k : ℝ)^2)=
      (n : ℝ)*(n+1)*(2*n+1) := by
  induction n with
  | zero => norm_num
  | succ n ih =>
    rw [Finset.sum_range_succ]
    push_cast
    nlinarith [ih]

theorem cubic_mode_frequency_trace (m c : ℝ) (n : ℕ) :
    (2*(n : ℝ)+1)^3*m^2+
      3*(2*(n : ℝ)+1)^2*c^2*(2*∑ k ∈ Finset.range (n+1), (k : ℝ)^2)=
      (2*(n : ℝ)+1)^3*(m^2+c^2*(n : ℝ)*(n+1)) := by
  have hs := positive_integer_square_sum n
  linear_combination (2*(n : ℝ)+1)^2*c^2*hs

theorem cutoff_variance_remainder (m c n : ℝ) :
    cutoffVariance m c n-c^2*n^5=
      m^2*(2*n+1)^3/8+c^2*(20*n^4+18*n^3+7*n^2+n)/8 := by
  dsimp [cutoffVariance]
  ring

theorem cutoff_variance_lower_bound (m c n : ℝ) (hn : 0 ≤ n) :
    c^2*n^5 ≤ cutoffVariance m c n := by
  have hrem := cutoff_variance_remainder m c n
  have hnonneg : 0 ≤ m^2*(2*n+1)^3/8+c^2*(20*n^4+18*n^3+7*n^2+n)/8 := by
    positivity
  linarith

theorem cutoff_variance_diverges (m c : ℝ) (hc : c ≠ 0) :
    Tendsto (cutoffVariance m c) atTop atTop := by
  have hp : Tendsto (fun n : ℝ => n^5) atTop atTop :=
    tendsto_pow_atTop (by norm_num : (5 : ℕ) ≠ 0)
  have hscaled : Tendsto (fun n : ℝ => c^2*n^5) atTop atTop :=
    Tendsto.const_mul_atTop (sq_pos_of_ne_zero hc) hp
  apply tendsto_atTop_mono' atTop ?_ hscaled
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with n hn
  exact cutoff_variance_lower_bound m c n hn

theorem natural_cutoff_variance_diverges (m c : ℝ) (hc : c ≠ 0) :
    Tendsto (fun n : ℕ => cutoffVariance m c n) atTop atTop := by
  exact (cutoff_variance_diverges m c hc).comp tendsto_natCast_atTop_atTop

theorem centered_sequence_norm_unbounded (m c : ℝ) (hc : c ≠ 0)
    (normSquared : ℕ → ℝ)
    (hnorm : ∀ n : ℕ, normSquared n=cutoffVariance m c n) :
    Tendsto normSquared atTop atTop := by
  have heq : normSquared=(fun n : ℕ => cutoffVariance m c n) := funext hnorm
  rw [heq]
  exact natural_cutoff_variance_diverges m c hc

theorem no_cauchy_from_variance_lower_bound {E : Type*}
    [NormedAddCommGroup E] [CompleteSpace E]
    (f : ℕ → E) (m c : ℝ) (hc : c ≠ 0)
    (hbound : ∀ n : ℕ, cutoffVariance m c n ≤ ‖f n‖^2) :
    ¬CauchySeq f := by
  intro hcauchy
  obtain ⟨limit,hlimit⟩ := cauchySeq_tendsto_of_complete hcauchy
  have hfinite : Tendsto (fun n => ‖f n‖^2) atTop (𝓝 (‖limit‖^2)) :=
    hlimit.norm.pow 2
  have hinfinite : Tendsto (fun n => ‖f n‖^2) atTop atTop :=
    tendsto_atTop_mono hbound (natural_cutoff_variance_diverges m c hc)
  exact not_tendsto_nhds_of_tendsto_atTop hinfinite _ hfinite

end
end LocalVariance
