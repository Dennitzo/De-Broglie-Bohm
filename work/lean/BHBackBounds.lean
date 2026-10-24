import Mathlib

/-! Conditional comparison and curvature bounds. Existence and invariant
domain arguments for the model ODE are not assumed to be kernel results. -/

namespace BHBackBounds

noncomputable section

open Filter
open scoped Topology

theorem weighted_D_hasDerivAt (D : ℝ → ℝ) (s dD : ℝ)
    (hD : HasDerivAt D dD s) :
    HasDerivAt (fun x => D x*Real.exp (4*x/3))
      ((dD+4*D s/3)*Real.exp (4*s/3)) s := by
  have he := (Real.hasDerivAt_exp (4*s/3)).comp s
    (((hasDerivAt_id s).const_mul 4).div_const 3)
  convert hD.mul he using 1 <;> dsimp <;> ring

theorem weighted_D_derivative_nonpositive (D : ℝ → ℝ) (s dD : ℝ)
    (hD : HasDerivAt D dD s) (hdecay : dD ≤ -4*D s/3) :
    deriv (fun x => D x*Real.exp (4*x/3)) s ≤ 0 := by
  rw [(weighted_D_hasDerivAt D s dD hD).deriv]
  exact mul_nonpos_of_nonpos_of_nonneg (by linarith) (le_of_lt (Real.exp_pos _))

theorem weighted_D_comparison (D dD : ℝ → ℝ) (s : ℝ) (hs : 0 ≤ s)
    (hD : ∀ x ∈ Set.Icc 0 s, HasDerivAt D (dD x) x)
    (hdecay : ∀ x ∈ Set.Icc 0 s, dD x ≤ -4*D x/3) :
    D s*Real.exp (4*s/3) ≤ D 0 := by
  have hw : ∀ x ∈ Set.Icc 0 s, HasDerivAt (fun z => D z*Real.exp (4*z/3))
      ((dD x+4*D x/3)*Real.exp (4*x/3)) x := by
    intro x hx
    exact weighted_D_hasDerivAt D x (dD x) (hD x hx)
  have hm : AntitoneOn (fun x => D x*Real.exp (4*x/3)) (Set.Icc 0 s) := by
    apply antitoneOn_of_deriv_nonpos (convex_Icc 0 s)
    · intro x hx
      exact (hw x hx).continuousAt.continuousWithinAt
    · intro x hx
      exact (hw x (interior_subset hx)).differentiableAt.differentiableWithinAt
    · intro x hx
      exact weighted_D_derivative_nonpositive D x (dD x) (hD x (interior_subset hx))
        (hdecay x (interior_subset hx))
  have hh := hm (by exact ⟨le_rfl, hs⟩) (by exact ⟨hs, le_rfl⟩) hs
  simpa using hh

theorem D_exponential_upper_bound (D dD : ℝ → ℝ) (s alpha : ℝ) (hs : 0 ≤ s)
    (hinit : D 0=alpha)
    (hD : ∀ x ∈ Set.Icc 0 s, HasDerivAt D (dD x) x)
    (hdecay : ∀ x ∈ Set.Icc 0 s, dD x ≤ -4*D x/3) :
    D s ≤ alpha*Real.exp (-4*s/3) := by
  have hh := weighted_D_comparison D dD s hs hD hdecay
  rw [hinit] at hh
  have he : Real.exp (-4*s/3)=1/Real.exp (4*s/3) := by
    rw [show -4*s/3= -(4*s/3) by ring, Real.exp_neg]
    simp only [one_div]
  rw [he, ← div_eq_mul_one_div]
  exact (le_div_iff₀ (Real.exp_pos _)).2 hh

theorem radius_square_hasDerivAt (r : ℝ → ℝ) (s q : ℝ)
    (hr : HasDerivAt r q s) :
    HasDerivAt (fun x => (r x)^2) (2*r s*q) s := by
  convert hr.pow 2 using 1 <;> dsimp <;> ring

theorem radius_square_rate_upper_bound (r : ℝ → ℝ) (s q alpha : ℝ)
    (hr : HasDerivAt r q s) (hr0 : 0 < r s) (hq : q ≤ alpha/r s) :
    deriv (fun x => (r x)^2) s ≤ 2*alpha := by
  rw [(radius_square_hasDerivAt r s q hr).deriv]
  have hh := (le_div_iff₀ hr0).mp hq
  nlinarith

theorem curvature_bound_transfer (alpha r D q A R E K : ℝ)
    (hr : 0 < r) (hD : 0 < D) (hE : 0 < E) (hA : 0 ≤ A)
    (hq : 0 ≤ q) (hqlower : alpha*A ≤ 2*r*q) (hDupper : D*E ≤ alpha)
    (hR : r^2 ≤ R) (hK : 8*q^2/(r^2*D^2) ≤ K) :
    2*A^2*E^2/R^2 ≤ K := by
  have hR0 : 0 < R := lt_of_lt_of_le (sq_pos_of_pos hr) hR
  have hp := mul_le_mul_of_nonneg_right hDupper hA
  have hlin : D*E*A ≤ 2*r*q := by linarith
  have hleft : 0 ≤ D*E*A := by positivity
  have hright : 0 ≤ 2*r*q := by positivity
  have hsq : (D*E*A)^2 ≤ (2*r*q)^2 := by nlinarith [sq_nonneg (2*r*q-D*E*A)]
  have hprod := mul_le_mul_of_nonneg_left hsq (by positivity : 0 ≤ 2*r^2)
  have hmid : 2*A^2*E^2/r^4 ≤ 8*q^2/(r^2*D^2) := by
    apply (div_le_div_iff₀ (by positivity : 0 < r^4) (by positivity : 0 < r^2*D^2)).2
    nlinarith [hprod]
  have hr4 : r^4 ≤ R^2 := by nlinarith [sq_nonneg (R-r^2)]
  have hfirst : 2*A^2*E^2/R^2 ≤ 2*A^2*E^2/r^4 := by
    apply (div_le_div_iff₀ (by positivity : 0 < R^2) (by positivity : 0 < r^4)).2
    exact mul_le_mul_of_nonneg_left hr4 (by positivity)
  exact le_trans hfirst (le_trans hmid hK)

def curvatureLowerBound (alpha r0 s : ℝ) : ℝ :=
  2*(1-Real.exp (-2*s))^2*Real.exp (8*s/3)/(r0^2+2*alpha*s)^2

theorem exponential_curvature_bound (alpha r0 s r D q K : ℝ)
    (hr : 0 < r) (hD : 0 < D) (hs : 0 ≤ s) (hq : 0 ≤ q)
    (hqlower : alpha*(1-Real.exp (-2*s)) ≤ 2*r*q)
    (hDupper : D*Real.exp (4*s/3) ≤ alpha)
    (hR : r^2 ≤ r0^2+2*alpha*s) (hK : 8*q^2/(r^2*D^2) ≤ K) :
    curvatureLowerBound alpha r0 s ≤ K := by
  have hA : 0 ≤ 1-Real.exp (-2*s) := by
    have he : Real.exp (-2*s) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
    linarith
  have hh := curvature_bound_transfer alpha r D q (1-Real.exp (-2*s))
    (r0^2+2*alpha*s) (Real.exp (4*s/3)) K hr hD (Real.exp_pos _) hA hq hqlower hDupper hR hK
  have he : (Real.exp (4*s/3))^2=Real.exp (8*s/3) := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  simpa only [he, curvatureLowerBound] using hh

theorem radius_ratio_tendsto (alpha r0 : ℝ) (ha : 0 < alpha) :
    Tendsto (fun s : ℝ => s/(r0^2+2*alpha*s)) atTop (𝓝 (1/(2*alpha))) := by
  have hq : Tendsto (fun s : ℝ => r0^2/s) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop tendsto_id
  have hd : Tendsto (fun s : ℝ => r0^2/s+2*alpha) atTop (𝓝 (2*alpha)) := by
    simpa using hq.add_const (2*alpha)
  have hh : Tendsto (fun s : ℝ => 1/(r0^2/s+2*alpha)) atTop (𝓝 (1/(2*alpha))) :=
    tendsto_const_nhds.div hd (by positivity)
  apply hh.congr'
  filter_upwards [eventually_gt_atTop (0:ℝ)] with s hs
  have hs0 : s ≠ 0 := ne_of_gt hs
  have hn : r0^2+2*alpha*s ≠ 0 := by positivity
  field_simp
  <;> ring

theorem curvature_lower_bound_tendsto (alpha r0 : ℝ) (ha : 0 < alpha) :
    Tendsto (curvatureLowerBound alpha r0) atTop atTop := by
  have hneg : Tendsto (fun s : ℝ => (-2)*s) atTop atBot :=
    Tendsto.const_mul_atTop_of_neg (by norm_num) tendsto_id
  have he0 : Tendsto (fun s : ℝ => Real.exp (-2*s)) atTop (𝓝 0) :=
    Real.tendsto_exp_atBot.comp hneg
  have hA : Tendsto (fun s : ℝ => 1-Real.exp (-2*s)) atTop (𝓝 1) := by
    simpa using tendsto_const_nhds.sub he0
  have hR := radius_ratio_tendsto alpha r0 ha
  have hB : Tendsto (fun s : ℝ => (8/3:ℝ)*(s/(r0^2+2*alpha*s))) atTop
      (𝓝 ((8/3:ℝ)*(1/(2*alpha)))) := hR.const_mul _
  have hfac : Tendsto (fun s : ℝ => 2*(1-Real.exp (-2*s))^2*
      ((8/3:ℝ)*(s/(r0^2+2*alpha*s)))^2) atTop
      (𝓝 (2*((8/3:ℝ)*(1/(2*alpha)))^2)) := by
    simpa using ((hA.pow 2).const_mul 2).mul (hB.pow 2)
  have hx : Tendsto (fun s : ℝ => (8/3:ℝ)*s) atTop atTop :=
    Tendsto.const_mul_atTop (by norm_num) tendsto_id
  have he : Tendsto (fun s : ℝ => Real.exp ((8/3:ℝ)*s)/((8/3:ℝ)*s)^2) atTop atTop :=
    (Real.tendsto_exp_div_pow_atTop 2).comp hx
  have hpos : 0 < 2*((8/3:ℝ)*(1/(2*alpha)))^2 := by positivity
  have hh := Tendsto.pos_mul_atTop hpos hfac he
  apply hh.congr'
  filter_upwards [eventually_gt_atTop (0:ℝ)] with s hs
  have hs0 : s ≠ 0 := ne_of_gt hs
  have hn : r0^2+2*alpha*s ≠ 0 := by positivity
  unfold curvatureLowerBound
  rw [show 8*s/3=(8/3:ℝ)*s by ring]
  field_simp
  <;> ring

theorem conditional_curvature_diverges (K : ℝ → ℝ) (alpha r0 : ℝ)
    (ha : 0 < alpha)
    (hK : ∀ᶠ s in atTop, curvatureLowerBound alpha r0 s ≤ K s) :
    Tendsto K atTop atTop := by
  exact tendsto_atTop_mono' atTop hK (curvature_lower_bound_tendsto alpha r0 ha)

theorem relative_decay_lower_bound (alpha r q s : ℝ)
    (ha : 0 < alpha) (hr2 : 4*alpha ≤ r^2)
    (hq : alpha*(1-Real.exp (-2*s)) ≤ 2*r*q) :
    -(2:ℝ)/3*Real.exp (-2*s) ≤ (4*r*q-2*alpha)/(r^2-alpha) := by
  have hd : 0 < r^2-alpha := by nlinarith
  apply (le_div_iff₀ hd).2
  have hp := mul_nonneg (le_of_lt (Real.exp_pos (-2*s)))
    (by linarith : 0 ≤ r^2-4*alpha)
  nlinarith [hp]

theorem corrected_log_D_hasDerivAt (D : ℝ → ℝ) (s dD : ℝ)
    (hD : HasDerivAt D dD s) (hD0 : D s ≠ 0) :
    HasDerivAt (fun x => Real.log (D x)+2*x+(1-Real.exp (-2*x))/3)
      (dD/D s+2+(2:ℝ)/3*Real.exp (-2*s)) s := by
  have hlog := (Real.hasDerivAt_log hD0).comp s hD
  have he := (Real.hasDerivAt_exp (-2*s)).comp s ((hasDerivAt_id s).const_mul (-2))
  have hall := (hlog.add ((hasDerivAt_id s).const_mul 2)).add
    (((hasDerivAt_const s 1).sub he).div_const 3)
  convert hall using 1 <;> dsimp <;> (try funext x) <;> (try dsimp) <;> ring

theorem weighted_D_lower_comparison (D dD : ℝ → ℝ) (s alpha : ℝ)
    (hs : 0 ≤ s) (ha : 0 < alpha) (hinit : D 0=alpha)
    (hpos : ∀ x ∈ Set.Icc 0 s, 0 < D x)
    (hD : ∀ x ∈ Set.Icc 0 s, HasDerivAt D (dD x) x)
    (hrate : ∀ x ∈ Set.Icc 0 s,
      0 ≤ dD x/D x+2+(2:ℝ)/3*Real.exp (-2*x)) :
    alpha*Real.exp (-(1:ℝ)/3) ≤ D s*Real.exp (2*s) := by
  have hw : ∀ x ∈ Set.Icc 0 s,
      HasDerivAt (fun z => Real.log (D z)+2*z+(1-Real.exp (-2*z))/3)
        (dD x/D x+2+(2:ℝ)/3*Real.exp (-2*x)) x := by
    intro x hx
    exact corrected_log_D_hasDerivAt D x (dD x) (hD x hx) (ne_of_gt (hpos x hx))
  have hm : MonotoneOn (fun x => Real.log (D x)+2*x+(1-Real.exp (-2*x))/3) (Set.Icc 0 s) := by
    apply monotoneOn_of_deriv_nonneg (convex_Icc 0 s)
    · intro x hx
      exact (hw x hx).continuousAt.continuousWithinAt
    · intro x hx
      exact (hw x (interior_subset hx)).differentiableAt.differentiableWithinAt
    · intro x hx
      rw [(hw x (interior_subset hx)).deriv]
      exact hrate x (interior_subset hx)
  have hh := hm (by exact ⟨le_rfl, hs⟩) (by exact ⟨hs, le_rfl⟩) hs
  simp only [hinit, mul_zero, Real.exp_zero, sub_self, zero_div, add_zero] at hh
  have hlog : Real.log alpha-(1:ℝ)/3 ≤ Real.log (D s)+2*s := by
    have he := Real.exp_pos (-2*s)
    linarith
  have he := Real.exp_le_exp.mpr hlog
  rw [Real.exp_sub, Real.exp_log ha, Real.exp_add, Real.exp_log (hpos s ⟨hs, le_rfl⟩)] at he
  simpa only [show -(1:ℝ)/3= -((1:ℝ)/3) by ring, Real.exp_neg, div_eq_mul_inv] using he

theorem optical_factor_lower_bound (alpha s D : ℝ) (ha : 0 < alpha) (hD : 0 < D)
    (hbound : alpha*Real.exp (-(1:ℝ)/3) ≤ D*Real.exp (2*s)) :
    Real.sqrt alpha*Real.exp (-(1:ℝ)/6) ≤ Real.exp s*Real.sqrt D := by
  have he : (Real.exp s)^2=Real.exp (2*s) := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  have hec : (Real.exp (-(1:ℝ)/6))^2=Real.exp (-(1:ℝ)/3) := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  have hsq : (Real.sqrt alpha*Real.exp (-(1:ℝ)/6))^2 ≤ (Real.exp s*Real.sqrt D)^2 := by
    rw [mul_pow, mul_pow, Real.sq_sqrt (le_of_lt ha), Real.sq_sqrt (le_of_lt hD), he, hec]
    nlinarith [hbound]
  have hl : 0 ≤ Real.sqrt alpha*Real.exp (-(1:ℝ)/6) := by positivity
  have hr : 0 ≤ Real.exp s*Real.sqrt D := by positivity
  nlinarith [hsq]

end

end BHBackBounds
