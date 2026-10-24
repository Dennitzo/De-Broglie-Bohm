import Mathlib

/-! Exact globally differentiable eigenmode trajectories for a free Gaussian.
The density-transport statement is an exact scalar Jacobian identity. -/
namespace FoundSynchronous
noncomputable section

def denominator (l t : ℝ) : ℝ := 1+l^2*t^2
def scale (l t : ℝ) : ℝ := Real.sqrt (denominator l t)
def trajectory (l z₀ t : ℝ) : ℝ := z₀*scale l t
def velocity (l t z : ℝ) : ℝ := l^2*t/denominator l t*z
def density (l t z : ℝ) : ℝ := Real.exp (-l*z^2/denominator l t)/scale l t

theorem denominator_ge_one (l t : ℝ) : 1 ≤ denominator l t := by
  dsimp [denominator]
  nlinarith [mul_nonneg (sq_nonneg l) (sq_nonneg t)]

theorem denominator_positive (l t : ℝ) : 0 < denominator l t := by
  have := denominator_ge_one l t
  linarith

theorem scale_positive (l t : ℝ) : 0 < scale l t := by
  exact Real.sqrt_pos.2 (denominator_positive l t)

theorem scale_squared (l t : ℝ) : (scale l t)^2=denominator l t := by
  exact Real.sq_sqrt (le_of_lt (denominator_positive l t))

theorem scale_at_zero (l : ℝ) : scale l 0=1 := by
  simp [scale,denominator]

theorem denominator_derivative (l t : ℝ) :
    HasDerivAt (denominator l) (2*l^2*t) t := by
  convert (((hasDerivAt_id t).pow 2).const_mul (l^2)).const_add 1 using 1 <;>
    dsimp [denominator] <;> ring

theorem scale_derivative (l t : ℝ) :
    HasDerivAt (scale l) (l^2*t/scale l t) t := by
  convert (denominator_derivative l t).sqrt (ne_of_gt (denominator_positive l t)) using 1
  dsimp [scale]
  ring

theorem trajectory_derivative (l z₀ t : ℝ) :
    HasDerivAt (trajectory l z₀) (z₀*(l^2*t/scale l t)) t := by
  exact (scale_derivative l t).const_mul z₀

theorem trajectory_initial (l z₀ : ℝ) : trajectory l z₀ 0=z₀ := by
  simp [trajectory,scale_at_zero]

theorem trajectory_solves_guidance (l z₀ t : ℝ) :
    deriv (trajectory l z₀) t=velocity l t (trajectory l z₀ t) := by
  rw [(trajectory_derivative l z₀ t).deriv]
  have hs := ne_of_gt (scale_positive l t)
  have hd := ne_of_gt (denominator_positive l t)
  have hsq := scale_squared l t
  dsimp [velocity,trajectory]
  field_simp
  rw [hsq]

theorem trajectory_differentiable (l z₀ : ℝ) : Differentiable ℝ (trajectory l z₀) := by
  intro t
  exact (trajectory_derivative l z₀ t).differentiableAt

theorem trajectory_continuous (l z₀ : ℝ) : Continuous (trajectory l z₀) :=
  (trajectory_differentiable l z₀).continuous

theorem trajectory_nonzero (l z₀ t : ℝ) (hz : z₀ ≠ 0) : trajectory l z₀ t ≠ 0 := by
  exact mul_ne_zero hz (ne_of_gt (scale_positive l t))

theorem trajectory_normalized_square (l z₀ t : ℝ) :
    (trajectory l z₀ t)^2/denominator l t=z₀^2 := by
  have hd := ne_of_gt (denominator_positive l t)
  dsimp [trajectory]
  rw [mul_pow,scale_squared]
  field_simp

theorem density_transport (l z₀ t : ℝ) :
    density l t (trajectory l z₀ t)*scale l t=Real.exp (-l*z₀^2) := by
  have hs := ne_of_gt (scale_positive l t)
  have hexp : -l*(trajectory l z₀ t)^2/denominator l t = -l*z₀^2 := by
    rw [mul_div_assoc,trajectory_normalized_square]
  dsimp [density]
  rw [hexp]
  exact div_mul_cancel₀ _ hs

theorem two_mode_density_transport (l μ x₀ y₀ t : ℝ) :
    (density l t (trajectory l x₀ t)*density μ t (trajectory μ y₀ t))*
    (scale l t*scale μ t)=Real.exp (-l*x₀^2)*Real.exp (-μ*y₀^2) := by
  calc
    _ = (density l t (trajectory l x₀ t)*scale l t)*
        (density μ t (trajectory μ y₀ t)*scale μ t) := by ring
    _ = _ := by rw [density_transport,density_transport]

theorem two_mode_jacobian_positive (l μ t : ℝ) : 0 < scale l t*scale μ t :=
  mul_pos (scale_positive l t) (scale_positive μ t)

theorem scale_covariance (l t : ℝ) :
    (scale l t)^2*(1/(2*l))=denominator l t/(2*l) := by
  rw [scale_squared]
  ring

theorem common_time_reparametrization (l z₀ s f' : ℝ) (f : ℝ → ℝ)
    (hf : HasDerivAt f f' s) :
    HasDerivAt (fun r => trajectory l z₀ (f r))
      (f'*velocity l (f s) (trajectory l z₀ (f s))) s := by
  have heq : z₀*(l^2*(f s)/scale l (f s)) =
      velocity l (f s) (trajectory l z₀ (f s)) := by
    rw [←(trajectory_derivative l z₀ (f s)).deriv]
    exact trajectory_solves_guidance l z₀ (f s)
  convert (trajectory_derivative l z₀ (f s)).comp s hf using 1
  rw [heq]
  ring

end
end FoundSynchronous
