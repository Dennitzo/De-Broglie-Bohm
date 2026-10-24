import Mathlib

/-! Regular square-radius chart and a derivative comparison ensuring a positive scale factor. -/
namespace DiracBackRegularity
noncomputable section
open Filter Topology

def rootFactor (C x : ℝ) : ℝ := Real.sqrt (2*C-x^2)
def vacuumScale (C x : ℝ) : ℝ := rootFactor C x/x
def polynomialControl (M mu kap x : ℝ) : ℝ := M*(kap*x+mu*x^3/3)

theorem rootFactor_hasDerivAt (C x : ℝ) (h : 0 < 2*C-x^2) :
    HasDerivAt (rootFactor C) (-x/rootFactor C x) x := by
  have hh := (((hasDerivAt_id x).pow 2).const_sub (2*C)).sqrt (ne_of_gt h)
  convert hh using 1 <;> dsimp [rootFactor] <;> field_simp <;> ring

theorem vacuumScale_hasDerivAt (C x : ℝ) (h : 0 < 2*C-x^2) (hx : x ≠ 0) :
    HasDerivAt (vacuumScale C) (-2*C/(rootFactor C x*x^2)) x := by
  have hw : rootFactor C x ≠ 0 := ne_of_gt (Real.sqrt_pos.2 h)
  have hs : (rootFactor C x)^2=2*C-x^2 := Real.sq_sqrt (le_of_lt h)
  convert (rootFactor_hasDerivAt C x h).div (hasDerivAt_id x) hx using 1
  dsimp
  field_simp
  nlinarith

theorem regular_chart_constraint (C x f fp gam mu kap S1 S3 : ℝ)
    (h : 0 < 2*C-x^2) (hx : x ≠ 0)
    (hf : fp=gam*(kap*S1+mu*x^2*S3)/(rootFactor C x)^3) :
    let w := rootFactor C x
    let A := w*f/x
    let b := -w/x
    let a := (-2*C*f/(w*x^2)+(w/x)*fp)/(-2*x^2/w)
    2*x^2*b*a+A*(b^2+1)=gam*(mu*S3+kap*S1/x^2) := by
  have hw : rootFactor C x ≠ 0 := ne_of_gt (Real.sqrt_pos.2 h)
  have hs : (rootFactor C x)^2=2*C-x^2 := Real.sq_sqrt (le_of_lt h)
  dsimp
  rw [hf]
  field_simp
  rw [hs]
  ring

theorem chart_time_radius_ratio (x w : ℝ) (hx : x ≠ 0) (hw : w ≠ 0) :
    (2*x)/(-2*x^2/w) = -w/x := by field_simp <;> ring

theorem bloch_component_bounds (S1 S2 S3 : ℝ) (h : S1^2+S2^2+S3^2=1) :
    |S1| ≤ 1 ∧ |S2| ≤ 1 ∧ |S3| ≤ 1 := by
  refine ⟨abs_le.mpr ⟨?_,?_⟩,abs_le.mpr ⟨?_,?_⟩,abs_le.mpr ⟨?_,?_⟩⟩ <;>
    nlinarith [sq_nonneg S1,sq_nonneg S2,sq_nonneg S3]

theorem regular_chart_rate_upper (gam mu kap x S1 S3 w wmin : ℝ)
    (hg : 0 ≤ gam) (hm : 0 ≤ mu) (hk : 0 ≤ kap)
    (h1 : S1 ≤ 1) (h3 : S3 ≤ 1) (hwmin : 0 < wmin) (hw : wmin ≤ w) :
    gam*(kap*S1+mu*x^2*S3)/w^3 ≤ (gam/wmin^3)*(kap+mu*x^2) := by
  have hwp : 0 < w := lt_of_lt_of_le hwmin hw
  have hpow : wmin^3 ≤ w^3 := pow_le_pow_left₀ (le_of_lt hwmin) hw 3
  have hnum : gam*(kap*S1+mu*x^2*S3) ≤ gam*(kap+mu*x^2) := by
    apply mul_le_mul_of_nonneg_left _ hg
    have hmx : 0 ≤ mu*x^2 := mul_nonneg hm (sq_nonneg x)
    nlinarith [mul_nonneg hk (sub_nonneg.mpr h1),mul_nonneg hmx (sub_nonneg.mpr h3)]
  calc
    gam*(kap*S1+mu*x^2*S3)/w^3 ≤ gam*(kap+mu*x^2)/w^3 :=
      div_le_div_of_nonneg_right hnum (le_of_lt (pow_pos hwp 3))
    _ ≤ gam*(kap+mu*x^2)/wmin^3 :=
      div_le_div_of_nonneg_left (by positivity) (pow_pos hwmin 3) hpow
    _ = (gam/wmin^3)*(kap+mu*x^2) := by ring

theorem polynomialControl_hasDerivAt (M mu kap x : ℝ) :
    HasDerivAt (polynomialControl M mu kap) (M*(kap+mu*x^2)) x := by
  convert (((hasDerivAt_id x).const_mul kap).add
    ((((hasDerivAt_id x).pow 3).const_mul mu).div_const 3)).const_mul M using 1 <;>
      dsimp [polynomialControl] <;> ring

theorem scale_factor_comparison (f fp : ℝ → ℝ) (M mu kap x0 x : ℝ)
    (hM : 0 ≤ M) (hm : 0 ≤ mu) (hk : 0 ≤ kap)
    (hx : x ∈ Set.Icc 0 x0)
    (hf : ∀ y ∈ Set.Icc 0 x0, HasDerivAt f (fp y) y)
    (hbound : ∀ y ∈ Set.Icc 0 x0, fp y ≤ M*(kap+mu*y^2)) :
    f x0-polynomialControl M mu kap x0 ≤ f x := by
  have hd : ∀ y ∈ Set.Icc 0 x0,
      HasDerivAt (fun z => polynomialControl M mu kap z-f z)
        (M*(kap+mu*y^2)-fp y) y := by
    intro y hy
    exact (polynomialControl_hasDerivAt M mu kap y).sub (hf y hy)
  have hmono : MonotoneOn (fun z => polynomialControl M mu kap z-f z) (Set.Icc 0 x0) := by
    apply monotoneOn_of_deriv_nonneg (convex_Icc 0 x0)
    · intro y hy
      exact (hd y hy).continuousAt.continuousWithinAt
    · intro y hy
      exact (hd y (interior_subset hy)).differentiableAt.differentiableWithinAt
    · intro y hy
      rw [(hd y (interior_subset hy)).deriv]
      exact sub_nonneg.mpr (hbound y (interior_subset hy))
  have hh := hmono hx ⟨le_trans hx.1 hx.2,le_rfl⟩ hx.2
  have hPx : 0 ≤ polynomialControl M mu kap x := by
    unfold polynomialControl
    have := hx.1
    positivity
  linarith

theorem scale_factor_positive (f fp : ℝ → ℝ) (M mu kap x0 x : ℝ)
    (hM : 0 ≤ M) (hm : 0 ≤ mu) (hk : 0 ≤ kap)
    (hx : x ∈ Set.Icc 0 x0)
    (hf : ∀ y ∈ Set.Icc 0 x0, HasDerivAt f (fp y) y)
    (hbound : ∀ y ∈ Set.Icc 0 x0, fp y ≤ M*(kap+mu*y^2))
    (hmargin : polynomialControl M mu kap x0 < f x0) : 0 < f x := by
  have hh := scale_factor_comparison f fp M mu kap x0 x hM hm hk hx hf hbound
  linarith

theorem scale_factor_endpoint_tendsto (f : ℝ → ℝ) (C : ℝ)
    (hf : ContinuousAt f 0) (hC : 0 < C) :
    Tendsto (fun x => rootFactor C x*f x) (𝓝[>] (0 : ℝ)) (𝓝 (Real.sqrt (2*C)*f 0)) := by
  have hw : ContinuousAt (rootFactor C) 0 := by unfold rootFactor; fun_prop
  simpa [rootFactor] using (hw.mul hf).tendsto.mono_left inf_le_left

end
end DiracBackRegularity
