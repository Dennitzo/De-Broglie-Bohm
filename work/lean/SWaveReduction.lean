import Mathlib

/-! Local spherical reduction for a real scalar. The spacetime tensor
construction and quantization are not part of these statements. -/
namespace SWaveReduction
noncomputable section
open Filter
open scoped Topology

def potential (m r : ℝ) : ℝ := (1-2*m/r)*(2*m/r^3)

theorem rescaled_field_hasDerivAt (r f : ℝ → ℝ) (x u v : ℝ)
    (hr : HasDerivAt r u x) (hf : HasDerivAt f v x) :
    HasDerivAt (fun y => r y*f y) (u*f x+r x*v) x := by
  convert hr.mul hf using 1 <;> ring

theorem spherical_flux_identity (r U : ℝ → ℝ) (x u v : ℝ)
    (hr : HasDerivAt r u x) (hU : HasDerivAt U v x) (hr0 : r x ≠ 0) :
    (r x)^2*deriv (fun y => U y/r y) x = r x*v-u*U x := by
  change (r x)^2*deriv (U/r) x = r x*v-u*U x
  rw [(hU.div hr hr0).deriv]
  field_simp
  <;> ring

theorem spherical_flux_hasDerivAt (r rp U Up : ℝ → ℝ) (x rpp Upp : ℝ)
    (hr : HasDerivAt r (rp x) x) (hrp : HasDerivAt rp rpp x)
    (hU : HasDerivAt U (Up x) x) (hUp : HasDerivAt Up Upp x) :
    HasDerivAt (fun y => r y*Up y-rp y*U y)
      (r x*Upp-rpp*U x) x := by
  convert (hr.mul hUp).sub (hrp.mul hU) using 1 <;> ring

theorem reduced_wave_operator (r rx rxx f fx fxx ftt : ℝ) (hr : r ≠ 0) :
    -(r*ftt)+(rxx*f+2*rx*fx+r*fxx)-(rxx/r)*(r*f) =
      r*(-ftt+fxx+2*rx/r*fx) := by
  field_simp
  <;> ring

theorem schwarzschild_radial_acceleration (r : ℝ → ℝ) (x m : ℝ)
    (hr : HasDerivAt r (1-2*m/r x) x) (hr0 : r x ≠ 0) :
    HasDerivAt (fun y => 1-2*m/r y)
      ((2*m/(r x)^2)*(1-2*m/r x)) x := by
  convert (hasDerivAt_const x (1 : ℝ)).sub
    ((hasDerivAt_const x (2*m)).div hr hr0) using 1 <;>
    dsimp <;> field_simp <;> ring

theorem potential_from_radial_acceleration (m r : ℝ) (hr : r ≠ 0) :
    ((2*m/r^2)*(1-2*m/r))/r=potential m r := by
  unfold potential
  field_simp
  <;> ring

theorem potential_positive (m r : ℝ) (hm : 0 < m) (hr : 2*m < r) :
    0 < potential m r := by
  have hr0 : 0 < r := by linarith
  unfold potential
  apply mul_pos
  · exact sub_pos.mpr ((div_lt_one hr0).2 hr)
  · positivity

theorem potential_hasDerivAt (m r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (potential m) (2*m*(8*m-3*r)/r^5) r := by
  have hi := hasDerivAt_id r
  have h1 := (hasDerivAt_const r (1 : ℝ)).sub
    ((hasDerivAt_const r (2*m)).div hi hr)
  have h2 := (hasDerivAt_const r (2*m)).div (hi.pow 3) (pow_ne_zero 3 hr)
  convert h1.mul h2 using 1 <;> dsimp [potential] <;>
    (try funext y) <;> field_simp <;> ring

theorem potential_peak_stationary (m : ℝ) (hm : 0 < m) :
    HasDerivAt (potential m) 0 (8*m/3) := by
  convert potential_hasDerivAt m (8*m/3) (by positivity) using 1 <;> ring

theorem potential_horizon_tendsto (m : ℝ) (hm : 0 < m) :
    Tendsto (potential m) (𝓝 (2*m)) (𝓝 0) := by
  have h := (potential_hasDerivAt m (2*m) (by positivity)).continuousAt.tendsto
  have hp : potential m (2*m)=0 := by
    unfold potential
    field_simp
    <;> ring
  simpa only [hp] using h


theorem second_zero_mode_hasDerivAt (r rp I : ℝ → ℝ) (x rpp : ℝ)
    (hr : HasDerivAt r (rp x) x) (hrp : HasDerivAt rp rpp x)
    (hI : HasDerivAt I (1/(r x)^2) x) (hr0 : r x ≠ 0) :
    HasDerivAt (fun y => rp y*I y+1/r y) ((rpp/r x)*(r x*I x)) x := by
  convert (hrp.mul hI).add ((hasDerivAt_const x (1 : ℝ)).div hr hr0) using 1 <;>
    dsimp <;> (try funext y) <;> (try dsimp) <;> field_simp <;> ring

theorem zero_mode_wronskian (r I : ℝ → ℝ) (x rp : ℝ)
    (hr : HasDerivAt r rp x) (hI : HasDerivAt I (1/(r x)^2) x) (hr0 : r x ≠ 0) :
    r x*deriv (fun y => r y*I y) x-rp*(r x*I x)=1 := by
  have hd : deriv (fun y => r y*I y) x=rp*I x+r x*(1/(r x)^2) := by
    simpa using (hr.mul hI).deriv
  rw [hd]
  field_simp
  <;> ring

end
end SWaveReduction
