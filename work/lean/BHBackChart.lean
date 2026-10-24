import Mathlib

/-! Redshift chart and reduced curvature identities for the specified
static model. The solution-domain inequalities are stated explicitly. -/

namespace BHBackChart

noncomputable section

def chartD (alpha r q : ℝ) : ℝ := alpha-2*r*q+q^2
def qRate (alpha r q : ℝ) : ℝ := (r-q)*(alpha-2*r*q)/(r^2-alpha)
def properURate (alpha r u v : ℝ) : ℝ := v*((r^2+alpha)*u+alpha*r*v)/(r^2-alpha)
def properVRate (alpha r u v : ℝ) : ℝ := -r*v*(2*u+r*v)/(r^2-alpha)
def properU (alpha r q : ℝ) : ℝ := -q/Real.sqrt (chartD alpha r q)
def properV (alpha r q : ℝ) : ℝ := 1/Real.sqrt (chartD alpha r q)
def reducedCurvature (r u v du dv : ℝ) : ℝ :=
  4*((dv+v^2)^2+2*(u*v/r)^2+2*(du/r)^2+((1-u^2)/r^2)^2)

theorem redshift_radial_rate (u v : ℝ) (hv : v ≠ 0) : u/(-v)= -u/v := by ring

theorem redshift_q_hasDerivAt (u v : ℝ → ℝ) (t alpha r : ℝ)
    (hv0 : v t ≠ 0) (hd : r^2-alpha ≠ 0)
    (hu : HasDerivAt u (properURate alpha r (u t) (v t)) t)
    (hv : HasDerivAt v (properVRate alpha r (u t) (v t)) t) :
    HasDerivAt (fun x => -u x/v x) (-v t*qRate alpha r (-u t/v t)) t := by
  convert hu.neg.div hv hv0 using 1 <;>
    dsimp [properURate, properVRate, qRate] <;> field_simp <;> ring

theorem redshift_q_rate_ratio (alpha r q v : ℝ) (hv : v ≠ 0) :
    (-v*qRate alpha r q)/(-v)=qRate alpha r q := by
  field_simp

theorem chart_D_hasDerivAt (r q : ℝ → ℝ) (s alpha : ℝ)
    (hd : (r s)^2-alpha ≠ 0)
    (hr : HasDerivAt r (q s) s)
    (hq : HasDerivAt q (qRate alpha (r s) (q s)) s) :
    HasDerivAt (fun x => chartD alpha (r x) (q x))
      (-2*r s*(r s-2*q s)*chartD alpha (r s) (q s)/((r s)^2-alpha)) s := by
  have hall := ((hasDerivAt_const s alpha).sub ((hr.mul hq).const_mul 2)).add (hq.pow 2)
  convert hall using 1 <;> dsimp [chartD, qRate] <;>
    (try funext x) <;> (try dsimp) <;> field_simp <;> ring

theorem chart_log_weight_hasDerivAt (r q : ℝ → ℝ) (s alpha : ℝ)
    (hd : (r s)^2-alpha ≠ 0) (hD : 0 < chartD alpha (r s) (q s))
    (hr : HasDerivAt r (q s) s)
    (hq : HasDerivAt q (qRate alpha (r s) (q s)) s) :
    HasDerivAt (fun x => Real.log (chartD alpha (r x) (q x))+2*x)
      ((4*r s*q s-2*alpha)/((r s)^2-alpha)) s := by
  have hD0 : chartD alpha (r s) (q s) ≠ 0 := ne_of_gt hD
  have hh := (Real.hasDerivAt_log hD0).comp s (chart_D_hasDerivAt r q s alpha hd hr hq)
  convert hh.add ((hasDerivAt_id s).const_mul 2) using 1 <;>
    dsimp <;> field_simp <;> ring

theorem normalized_constraint (alpha r q : ℝ) (hd : 0 < chartD alpha r q) :
    (properU alpha r q)^2+2*r*properU alpha r q*properV alpha r q+
      alpha*(properV alpha r q)^2=1 := by
  have hs0 : Real.sqrt (chartD alpha r q) ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hd)
  have hs2 := Real.sq_sqrt (le_of_lt hd)
  unfold properU properV
  field_simp
  dsimp [chartD] at *
  nlinarith [hs2]

theorem normalized_signs (alpha r q : ℝ) (hd : 0 < chartD alpha r q) (hq : 0 < q) :
    properU alpha r q < 0 ∧ 0 < properV alpha r q := by
  have hs : 0 < Real.sqrt (chartD alpha r q) := Real.sqrt_pos.2 hd
  constructor
  · unfold properU
    exact div_neg_of_neg_of_pos (neg_neg_of_pos hq) hs
  · unfold properV
    exact div_pos (by norm_num) hs

theorem curvature_component_lower_bound (r u v du dv : ℝ) :
    8*(u*v/r)^2 ≤ reducedCurvature r u v du dv := by
  unfold reducedCurvature
  nlinarith [sq_nonneg (dv+v^2), sq_nonneg (du/r), sq_nonneg ((1-u^2)/r^2)]

theorem chart_curvature_lower_bound (alpha r q du dv : ℝ)
    (hr : r ≠ 0) (hd : 0 < chartD alpha r q) :
    8*q^2/(r^2*(chartD alpha r q)^2) ≤
      reducedCurvature r (properU alpha r q) (properV alpha r q) du dv := by
  have hs0 : Real.sqrt (chartD alpha r q) ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hd)
  have hd0 : chartD alpha r q ≠ 0 := ne_of_gt hd
  have hs2 := Real.sq_sqrt (le_of_lt hd)
  have he : 8*(properU alpha r q*properV alpha r q/r)^2=
      8*q^2/(r^2*(chartD alpha r q)^2) := by
    unfold properU properV
    have hs4 : (Real.sqrt (chartD alpha r q))^4=(chartD alpha r q)^2 := by
      calc
        _ = ((Real.sqrt (chartD alpha r q))^2)^2 := by ring
        _ = _ := by rw [hs2]
    field_simp
    nlinarith [congrArg (fun x : ℝ => q^2*x) hs4]
  rw [← he]
  exact curvature_component_lower_bound r (properU alpha r q) (properV alpha r q) du dv

theorem inner_root_identity (alpha r : ℝ) (hr : 0 < r)
    (ha : 0 < alpha) (hd : alpha < r^2) :
    r-Real.sqrt (r^2-alpha)=alpha/(r+Real.sqrt (r^2-alpha)) := by
  have hs2 := Real.sq_sqrt (le_of_lt (sub_pos.mpr hd))
  have hs0 := Real.sqrt_nonneg (r^2-alpha)
  have hden : r+Real.sqrt (r^2-alpha) ≠ 0 := by positivity
  apply (eq_div_iff hden).2
  nlinarith [hs2]

theorem inner_root_bounds (alpha r : ℝ) (hr : 0 < r)
    (ha : 0 < alpha) (hd : alpha < r^2) :
    0 < r-Real.sqrt (r^2-alpha) ∧ r-Real.sqrt (r^2-alpha) < alpha/r := by
  rw [inner_root_identity alpha r hr ha hd]
  have hs : 0 < Real.sqrt (r^2-alpha) := Real.sqrt_pos.2 (sub_pos.mpr hd)
  constructor
  · exact div_pos ha (by positivity)
  · exact div_lt_div_of_pos_left ha hr (by linarith)

theorem chart_D_upper_bound (alpha r q : ℝ) (hq : 0 ≤ q) (hqr : q ≤ 2*r) :
    chartD alpha r q ≤ alpha := by
  have hp := mul_nonneg hq (sub_nonneg.mpr hqr)
  unfold chartD
  nlinarith

theorem chart_decay_coefficient_lower_bound (alpha r q : ℝ)
    (ha : 0 < alpha) (hr : 0 < r) (hr2 : 4*alpha ≤ r^2)
    (hq : q ≤ alpha/r) : (4:ℝ)/3 ≤ 2*r*(r-2*q)/(r^2-alpha) := by
  have hd : 0 < r^2-alpha := by nlinarith
  have hqr : q*r ≤ alpha := (le_div_iff₀ hr).mp hq
  apply (le_div_iff₀ hd).2
  nlinarith

end

end BHBackChart
