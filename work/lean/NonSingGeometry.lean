import Mathlib

/-! Regular areal radius, curvature bounds and the necessary source change. -/
namespace NonSingGeometry
noncomputable section

def radius (a r : ℝ) : ℝ := Real.sqrt (r^2+a^2)
def radiusPrime (a r : ℝ) : ℝ := r/radius a r
def normalizedCurvature (u v x : ℝ) : ℝ :=
  4*((u*(3*x-1))^2+2*(u*x)^2+2*(u*(2-3*x)-v*(1-x))^2+(2*u*x+v*(1-x))^2)
def curvature (M a r : ℝ) : ℝ :=
  normalizedCurvature (M/(radius a r)^3) (1/(radius a r)^2) (r^2/(radius a r)^2)

theorem radius_lower (a r : ℝ) (ha : 0 < a) : a ≤ radius a r := by
  unfold radius
  have hs := Real.sq_sqrt (show 0 ≤ r^2+a^2 by positivity)
  have hp := Real.sqrt_nonneg (r^2+a^2)
  nlinarith [sq_nonneg r]

theorem radius_hasDerivAt (a r : ℝ) (ha : 0 < a) :
    HasDerivAt (radius a) (radiusPrime a r) r := by
  have h : r^2+a^2 ≠ 0 := by positivity
  convert (((hasDerivAt_id r).pow 2).add_const (a^2)).sqrt h using 1 <;>
    dsimp [radius,radiusPrime] <;> field_simp <;> ring

theorem radiusPrime_hasDerivAt (a r : ℝ) (ha : 0 < a) :
    HasDerivAt (radiusPrime a) (a^2/(radius a r)^3) r := by
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  have hs : (radius a r)^2=r^2+a^2 := Real.sq_sqrt (by positivity)
  convert (hasDerivAt_id r).div (radius_hasDerivAt a r ha) (ne_of_gt hR) using 1
  dsimp [radiusPrime]
  field_simp
  nlinarith

theorem radius_throat_strict_minimum (a r : ℝ) (ha : 0 < a) (hr : r ≠ 0) :
    radius a 0=a ∧ a < radius a r := by
  constructor
  · simp [radius,Real.sqrt_sq (le_of_lt ha)]
  · unfold radius
    exact (Real.lt_sqrt (le_of_lt ha)).2 (by nlinarith [sq_pos_of_ne_zero hr])

theorem core_shift_lower_bound (M a r rc : ℝ) (hM : 0 ≤ M) (ha : 0 < a)
    (hr : r^2 ≤ rc^2) :
    Real.sqrt (2*M/radius a rc) ≤ Real.sqrt (2*M/radius a r) := by
  apply Real.sqrt_le_sqrt
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  have hrr : radius a r ≤ radius a rc := by
    unfold radius
    exact Real.sqrt_le_sqrt (by linarith)
  exact div_le_div_of_nonneg_left (by positivity) hR hrr

theorem core_shift_strictly_superluminal (M a rc : ℝ) (ha : 0 < a)
    (hcore : radius a rc < 2*M) : 1 < Real.sqrt (2*M/radius a rc) := by
  have hR : 0 < radius a rc := lt_of_lt_of_le ha (radius_lower a rc ha)
  have hrat : 1 < 2*M/radius a rc := (lt_div_iff₀ hR).2 (by linarith)
  have hs := Real.sq_sqrt (show 0 ≤ 2*M/radius a rc by linarith)
  have hp := Real.sqrt_nonneg (2*M/radius a rc)
  nlinarith

theorem normalized_curvature_bound (u v U V x : ℝ)
    (hu : 0 ≤ u) (hv : 0 ≤ v) (huU : u ≤ U) (hvV : v ≤ V)
    (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    normalizedCurvature u v x ≤ 4*(6*U^2+3*(2*U+V)^2) := by
  have hU : 0 ≤ U := le_trans hu huU
  have hV : 0 ≤ V := le_trans hv hvV
  have hux0 : 0 ≤ u*x := mul_nonneg hu hx0
  have hux1 : u*x ≤ u := by nlinarith [mul_nonneg hu (sub_nonneg.mpr hx1)]
  have hvx0 : 0 ≤ v*(1-x) := mul_nonneg hv (sub_nonneg.mpr hx1)
  have hvx1 : v*(1-x) ≤ v := by nlinarith [mul_nonneg hv hx0]
  have hA : |u*(3*x-1)| ≤ 2*U := abs_le.mpr ⟨by nlinarith,by nlinarith⟩
  have hB : |u*x| ≤ U := abs_le.mpr ⟨by nlinarith,by nlinarith⟩
  have hC : |u*(2-3*x)-v*(1-x)| ≤ 2*U+V := abs_le.mpr ⟨by nlinarith,by nlinarith⟩
  have hD : |2*u*x+v*(1-x)| ≤ 2*U+V := abs_le.mpr ⟨by nlinarith,by nlinarith⟩
  have hAs := (sq_le_sq₀ (abs_nonneg _) (by positivity)).2 hA
  have hBs := (sq_le_sq₀ (abs_nonneg _) hU).2 hB
  have hCs := (sq_le_sq₀ (abs_nonneg _) (by positivity)).2 hC
  have hDs := (sq_le_sq₀ (abs_nonneg _) (by positivity)).2 hD
  simp only [sq_abs] at hAs hBs hCs hDs
  unfold normalizedCurvature
  nlinarith

theorem curvature_uniform_bound (M a r : ℝ) (hM : 0 ≤ M) (ha : 0 < a) :
    curvature M a r ≤ 4*(6*(M/a^3)^2+3*(2*M/a^3+1/a^2)^2) := by
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  have hRa := radius_lower a r ha
  have hs : (radius a r)^2=r^2+a^2 := Real.sq_sqrt (by positivity)
  have hx : r^2/(radius a r)^2 ≤ 1 := by
    apply (div_le_one (sq_pos_of_pos hR)).2
    nlinarith [sq_nonneg a]
  have hu : M/(radius a r)^3 ≤ M/a^3 := by gcongr
  have hv : 1/(radius a r)^2 ≤ 1/a^2 := by gcongr
  have hh := normalized_curvature_bound (M/(radius a r)^3) (1/(radius a r)^2)
    (M/a^3) (1/a^2) (r^2/(radius a r)^2) (by positivity) (by positivity) hu hv (by positivity) hx
  simpa [curvature,mul_div_assoc] using hh

theorem ingoing_null_source_strictly_negative (a R v : ℝ)
    (ha : 0 < a) (hR : 0 < R) (hv : 0 ≤ v) :
    -2*a^2*(1+v)^2/R^4 < 0 := by
  apply div_neg_of_neg_of_pos
  · have h1 : 0 < (1+v)^2 := sq_pos_of_pos (by linarith)
    have h2 : 0 < a^2 := sq_pos_of_pos ha
    nlinarith [mul_pos h1 h2]
  · positivity

theorem old_zero_pressure_no_throat (B ddB : ℝ) (hB : 0 < B)
    (hEinstein : 2*B*ddB+1=0) : ddB < 0 := by
  by_contra hh
  have hdd : 0 ≤ ddB := le_of_not_gt hh
  nlinarith [mul_nonneg (le_of_lt hB) hdd]

theorem black_bounce_throat_requires_negative_pressure (a M : ℝ)
    (ha : 0 < a) (hblack : a < 2*M) : (a-4*M)/a^3 < 0 := by
  exact div_neg_of_neg_of_pos (by linarith) (by positivity)

end
end NonSingGeometry
