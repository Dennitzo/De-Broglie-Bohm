import Mathlib

/-! A genuine differential residence-time bound, conditional on a differentiable guided path. -/
namespace NonSingCrossing
noncomputable section

theorem positive_weights_radial_bounds (u d v : ℝ)
    (hu : 0 ≤ u) (hd : 0 ≤ d) (hq : 0 < u+d) :
    -1-v ≤ ((1-v)*u+(-1-v)*d)/(u+d) ∧
      ((1-v)*u+(-1-v)*d)/(u+d) ≤ 1-v := by
  constructor
  · apply (le_div_iff₀ hq).2
    nlinarith
  · apply (div_le_iff₀ hq).2
    nlinarith

theorem core_radial_velocity_bound (u d v vc : ℝ)
    (hu : 0 ≤ u) (hd : 0 ≤ d) (hq : 0 < u+d) (hv : vc ≤ v) :
    ((1-v)*u+(-1-v)*d)/(u+d) ≤ -(vc-1) := by
  have hh := (positive_weights_radial_bounds u d v hu hd hq).2
  linarith

theorem radial_differential_comparison (r dr : ℝ → ℝ) (t0 t1 c : ℝ)
    (ht : t0 ≤ t1)
    (hr : ∀ t ∈ Set.Icc t0 t1, HasDerivAt r (dr t) t)
    (hv : ∀ t ∈ Set.Icc t0 t1, dr t ≤ -c) :
    r t1 ≤ r t0-c*(t1-t0) := by
  have hder : ∀ t ∈ Set.Icc t0 t1,
      HasDerivAt (fun s => -r s-c*s) (-dr t-c) t := by
    intro t hmem
    convert (hr t hmem).neg.sub ((hasDerivAt_id t).const_mul c) using 1 <;> dsimp <;> ring
  have hm : MonotoneOn (fun s => -r s-c*s) (Set.Icc t0 t1) := by
    apply monotoneOn_of_deriv_nonneg (convex_Icc t0 t1)
    · intro t hmem
      exact (hder t hmem).continuousAt.continuousWithinAt
    · intro t hmem
      exact (hder t (interior_subset hmem)).differentiableAt.differentiableWithinAt
    · intro t hmem
      rw [(hder t (interior_subset hmem)).deriv]
      linarith [hv t (interior_subset hmem)]
  have hh := hm ⟨le_rfl,ht⟩ ⟨ht,le_rfl⟩ ht
  nlinarith

theorem core_residence_time_bound (r dr : ℝ → ℝ) (t0 t1 rc vc : ℝ)
    (ht : t0 ≤ t1) (hvc : 1 < vc)
    (hr : ∀ t ∈ Set.Icc t0 t1, HasDerivAt r (dr t) t)
    (hv : ∀ t ∈ Set.Icc t0 t1, dr t ≤ -(vc-1))
    (hstart : r t0 ≤ rc) (hend : -rc ≤ r t1) :
    t1-t0 ≤ 2*rc/(vc-1) := by
  have hh := radial_differential_comparison r dr t0 t1 (vc-1) ht hr hv
  apply (le_div_iff₀ (by linarith : 0 < vc-1)).2
  nlinarith

theorem reaches_opposite_boundary_by_deadline (r dr : ℝ → ℝ) (t0 rc vc : ℝ)
    (hrc : 0 ≤ rc) (hvc : 1 < vc)
    (hr : ∀ t ∈ Set.Icc t0 (t0+2*rc/(vc-1)), HasDerivAt r (dr t) t)
    (hv : ∀ t ∈ Set.Icc t0 (t0+2*rc/(vc-1)), dr t ≤ -(vc-1))
    (hstart : r t0 ≤ rc) : r (t0+2*rc/(vc-1)) ≤ -rc := by
  have hvp : 0 < vc-1 := by linarith
  have ht : t0 ≤ t0+2*rc/(vc-1) := by
    have hh : 0 ≤ 2*rc/(vc-1) := by positivity
    linarith
  have hh := radial_differential_comparison r dr t0 (t0+2*rc/(vc-1)) (vc-1) ht hr hv
  have he : (vc-1)*((t0+2*rc/(vc-1))-t0)=2*rc := by field_simp; ring
  rw [he] at hh
  linarith

end
end NonSingCrossing
