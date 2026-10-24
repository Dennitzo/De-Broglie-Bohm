import Mathlib

/-! Two-dimensional conformal stress tensor on the explicitly specified
static metric. These are reduced 2D identities, not a four-dimensional
stress tensor or a dimensional reduction of its quantum state. -/

namespace BlackHoleStress

noncomputable section

def lapse (m r : ℝ) : ℝ := 1-2*m/r
def lapseFirst (m r : ℝ) : ℝ := 2*m/r^2
def lapseSecond (m r : ℝ) : ℝ := -4*m/r^3
def lapseThird (m r : ℝ) : ℝ := 12*m/r^4

def geometricB (f f1 f2 : ℝ) : ℝ := (2*f*f2-f1^2)/(192*Real.pi)
def mixedStress (f f2 : ℝ) : ℝ := f*f2/(96*Real.pi)
def radialB (m r : ℝ) : ℝ := geometricB (lapse m r) (lapseFirst m r) (lapseSecond m r)
def radialMixed (m r : ℝ) : ℝ := mixedStress (lapse m r) (lapseSecond m r)
def stressUU (m r tu : ℝ) : ℝ := radialB m r+tu
def stressVV (m r tv : ℝ) : ℝ := radialB m r+tv
def surfaceGravity (m : ℝ) : ℝ := 1/(4*m)
def thermalFlux (m : ℝ) : ℝ := (surfaceGravity m)^2/(48*Real.pi)
def horizonRegularUU (m r : ℝ) : ℝ :=
  (r^2+4*m*r+12*m^2)/(768*Real.pi*m^2*r^2)

theorem lapse_hasDerivAt (m r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (lapse m) (lapseFirst m r) r := by
  convert (hasDerivAt_const r 1).sub
    ((hasDerivAt_const r (2*m)).div (hasDerivAt_id r) hr) using 1 <;>
    dsimp [lapse, lapseFirst] <;> field_simp <;> ring

theorem lapse_first_hasDerivAt (m r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (lapseFirst m) (lapseSecond m r) r := by
  have hr2 : r^2 ≠ 0 := pow_ne_zero 2 hr
  convert (hasDerivAt_const r (2*m)).div ((hasDerivAt_id r).pow 2) hr2 using 1 <;>
    dsimp [lapseFirst, lapseSecond] <;> field_simp <;> ring

theorem lapse_second_hasDerivAt (m r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (lapseSecond m) (lapseThird m r) r := by
  have hr3 : r^3 ≠ 0 := pow_ne_zero 3 hr
  convert (hasDerivAt_const r (-4*m)).div ((hasDerivAt_id r).pow 3) hr3 using 1 <;>
    dsimp [lapseSecond, lapseThird] <;> field_simp <;> ring

theorem geometric_B_hasDerivAt (f f1 f2 : ℝ → ℝ) (x f3 : ℝ)
    (hf : HasDerivAt f (f1 x) x) (hf1 : HasDerivAt f1 (f2 x) x)
    (hf2 : HasDerivAt f2 f3 x) :
    HasDerivAt (fun u => geometricB (f u) (f1 u) (f2 u))
      (f x*f3/(96*Real.pi)) x := by
  convert (((hf.mul hf2).const_mul 2).sub (hf1.pow 2)).div_const (192*Real.pi) using 1 <;>
    dsimp [geometricB] <;> (try funext u) <;> ring

theorem mixed_stress_hasDerivAt (f f2 : ℝ → ℝ) (x f1 f3 : ℝ)
    (hf : HasDerivAt f f1 x) (hf2 : HasDerivAt f2 f3 x) :
    HasDerivAt (fun u => mixedStress (f u) (f2 u))
      ((f1*f2 x+f x*f3)/(96*Real.pi)) x := by
  convert (hf.mul hf2).div_const (96*Real.pi) using 1 <;>
    dsimp [mixedStress] <;> ring

theorem reduced_ward_algebra (f f1 f2 f3 : ℝ) :
    f*(f*f3/(96*Real.pi)-(f1*f2+f*f3)/(96*Real.pi))+
      f1*mixedStress f f2=0 := by
  unfold mixedStress
  ring

theorem reduced_ward_identity (f f1 f2 : ℝ → ℝ) (x f3 : ℝ)
    (hf : HasDerivAt f (f1 x) x) (hf1 : HasDerivAt f1 (f2 x) x)
    (hf2 : HasDerivAt f2 f3 x) :
    f x*(deriv (fun u => geometricB (f u) (f1 u) (f2 u)) x-
      deriv (fun u => mixedStress (f u) (f2 u)) x)+f1 x*mixedStress (f x) (f2 x)=0 := by
  rw [(geometric_B_hasDerivAt f f1 f2 x f3 hf hf1 hf2).deriv,
    (mixed_stress_hasDerivAt f f2 x (f1 x) f3 hf hf2).deriv]
  exact reduced_ward_algebra (f x) (f1 x) (f2 x) f3

theorem radial_ward_identity (m r : ℝ) (hr : r ≠ 0) :
    lapse m r*(deriv (radialB m) r-deriv (radialMixed m) r)+
      lapseFirst m r*radialMixed m r=0 := by
  exact reduced_ward_identity (lapse m) (lapseFirst m) (lapseSecond m) r
    (lapseThird m r) (lapse_hasDerivAt m r hr) (lapse_first_hasDerivAt m r hr)
    (lapse_second_hasDerivAt m r hr)

theorem conformal_trace_identity (f f2 : ℝ) (hf : f ≠ 0) :
    -4*mixedStress f f2/f= -f2/(24*Real.pi) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  unfold mixedStress
  field_simp
  <;> ring

theorem radial_B_reduced (m r : ℝ) (hr : r ≠ 0) :
    radialB m r=m*(3*m-2*r)/(48*Real.pi*r^4) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  unfold radialB geometricB lapse lapseFirst lapseSecond
  field_simp
  <;> ring

theorem stationary_net_flux (m r tu tv : ℝ) :
    stressUU m r tu-stressVV m r tv=tu-tv := by
  unfold stressUU stressVV
  ring

theorem thermal_flux_reduced (m : ℝ) : thermalFlux m=1/(768*Real.pi*m^2) := by
  unfold thermalFlux surfaceGravity
  ring

theorem thermal_flux_positive (m : ℝ) (hm : 0 < m) : 0 < thermalFlux m := by
  rw [thermal_flux_reduced]
  positivity

theorem horizon_B_cancellation (m : ℝ) (hm : m ≠ 0) :
    radialB m (2*m)= -thermalFlux m := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  rw [radial_B_reduced m (2*m) (mul_ne_zero (by norm_num) hm), thermal_flux_reduced]
  field_simp
  <;> ring

theorem boulware_horizon_negative (m : ℝ) (hm : 0 < m) :
    stressUU m (2*m) 0 < 0 := by
  unfold stressUU
  rw [horizon_B_cancellation m (ne_of_gt hm)]
  have hh := thermal_flux_positive m hm
  linarith

theorem unruh_outgoing_horizon_zero (m : ℝ) (hm : m ≠ 0) :
    stressUU m (2*m) (thermalFlux m)=0 := by
  unfold stressUU
  rw [horizon_B_cancellation m hm]
  ring

theorem unruh_outgoing_factorization (m r : ℝ) (hm : m ≠ 0) (hr : r ≠ 0) :
    stressUU m r (thermalFlux m)=
      (r-2*m)^2*(r^2+4*m*r+12*m^2)/(768*Real.pi*m^2*r^4) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  unfold stressUU
  rw [radial_B_reduced m r hr, thermal_flux_reduced]
  field_simp
  <;> ring

theorem unruh_regular_quotient (m r : ℝ) (hm : m ≠ 0) (hr : r ≠ 0)
    (hh : r-2*m ≠ 0) :
    stressUU m r (thermalFlux m)/(lapse m r)^2=horizonRegularUU m r := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  rw [unruh_outgoing_factorization m r hm hr]
  unfold lapse horizonRegularUU
  field_simp
  <;> ring

theorem regular_quotient_horizon_value (m : ℝ) (hm : m ≠ 0) :
    horizonRegularUU m (2*m)=1/(128*Real.pi*m^2) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  unfold horizonRegularUU
  field_simp
  <;> ring

theorem regular_quotient_continuous_at_horizon (m : ℝ) (hm : m ≠ 0) :
    ContinuousAt (horizonRegularUU m) (2*m) := by
  unfold horizonRegularUU
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  have hd : 768*Real.pi*m^2*(2*m)^2 ≠ 0 := by positivity
  fun_prop

theorem hartle_hawking_flux_zero (m r : ℝ) :
    stressUU m r (thermalFlux m)-stressVV m r (thermalFlux m)=0 := by
  rw [stationary_net_flux]
  ring

theorem unruh_flux_positive (m r : ℝ) (hm : 0 < m) :
    0 < stressUU m r (thermalFlux m)-stressVV m r 0 := by
  rw [stationary_net_flux, sub_zero]
  exact thermal_flux_positive m hm

end

end BlackHoleStress
