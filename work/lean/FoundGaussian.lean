import Mathlib

/-! A genuine two-mode Gaussian, its derivatives and finite quantum energy.
No convergence or infinite-dimensional measure is asserted. -/
namespace FoundGaussian
noncomputable section

def exponent (a b c x y : ℝ) : ℝ := -(a*x^2+2*b*x*y+c*y^2)/2
def amplitude (a b c x y : ℝ) : ℝ := Real.exp (exponent a b c x y)
def quantum1 (a b x y : ℝ) : ℝ := a/2-(a*x+b*y)^2/2
def quantum2 (b c x y : ℝ) : ℝ := c/2-(b*x+c*y)^2/2
def potential (a b c x y : ℝ) : ℝ :=
  ((a^2+b^2)*x^2+2*b*(a+c)*x*y+(b^2+c^2)*y^2)/2

theorem exponent_dx (a b c x y : ℝ) :
    HasDerivAt (fun z => exponent a b c z y) (-(a*x+b*y)) x := by
  convert ((((hasDerivAt_id x).pow 2).const_mul a).add
    (((hasDerivAt_id x).const_mul (2*b)).mul_const y)).add_const
      (c*y^2) |>.neg.div_const 2 using 1 <;>
    dsimp [exponent] <;> ring

theorem exponent_dy (a b c x y : ℝ) :
    HasDerivAt (fun z => exponent a b c x z) (-(b*x+c*y)) y := by
  convert (((((hasDerivAt_id y).const_mul (2*b*x)).add
    (((hasDerivAt_id y).pow 2).const_mul c)).const_add (a*x^2)).neg.div_const 2)
      using 1 <;> dsimp [exponent] <;> try funext z
  all_goals ring

theorem amplitude_positive (a b c x y : ℝ) : 0 < amplitude a b c x y :=
  Real.exp_pos _

theorem amplitude_dx (a b c x y : ℝ) :
    HasDerivAt (fun z => amplitude a b c z y)
      (-(a*x+b*y)*amplitude a b c x y) x := by
  convert (exponent_dx a b c x y).exp using 1 <;> dsimp [amplitude] <;> ring

theorem amplitude_dy (a b c x y : ℝ) :
    HasDerivAt (fun z => amplitude a b c x z)
      (-(b*x+c*y)*amplitude a b c x y) y := by
  convert (exponent_dy a b c x y).exp using 1 <;> dsimp [amplitude] <;> ring

theorem amplitude_dxx (a b c x y : ℝ) :
    HasDerivAt (fun z => deriv (fun w => amplitude a b c w y) z)
      (((a*x+b*y)^2-a)*amplitude a b c x y) x := by
  have heq : (fun z => deriv (fun w => amplitude a b c w y) z)=
      (fun z => -(a*z+b*y)*amplitude a b c z y) := by
    funext z
    exact (amplitude_dx a b c z y).deriv
  rw [heq]
  convert ((((hasDerivAt_id x).const_mul a).add_const (b*y)).neg).mul
    (amplitude_dx a b c x y) using 1 <;> dsimp only [Pi.neg_apply,id_eq] <;> ring

theorem amplitude_dyy (a b c x y : ℝ) :
    HasDerivAt (fun z => deriv (fun w => amplitude a b c x w) z)
      (((b*x+c*y)^2-c)*amplitude a b c x y) y := by
  have heq : (fun z => deriv (fun w => amplitude a b c x w) z)=
      (fun z => -(b*x+c*z)*amplitude a b c x z) := by
    funext z
    exact (amplitude_dy a b c x z).deriv
  rw [heq]
  convert ((((hasDerivAt_id y).const_mul c).const_add (b*x)).neg).mul
    (amplitude_dy a b c x y) using 1 <;> dsimp only [Pi.neg_apply,id_eq] <;> ring

theorem quantum1_from_amplitude (a b c x y : ℝ) :
    -deriv (fun z => deriv (fun w => amplitude a b c w y) z) x/
      (2*amplitude a b c x y)=quantum1 a b x y := by
  rw [(amplitude_dxx a b c x y).deriv]
  have hr := ne_of_gt (amplitude_positive a b c x y)
  dsimp [quantum1]
  field_simp
  <;> ring

theorem quantum2_from_amplitude (a b c x y : ℝ) :
    -deriv (fun z => deriv (fun w => amplitude a b c x w) z) y/
      (2*amplitude a b c x y)=quantum2 b c x y := by
  rw [(amplitude_dyy a b c x y).deriv]
  have hr := ne_of_gt (amplitude_positive a b c x y)
  dsimp [quantum2]
  field_simp
  <;> ring


theorem quantum1_dy (a b x y : ℝ) :
    HasDerivAt (fun z => quantum1 a b x z) (-b*(a*x+b*y)) y := by
  convert (((((hasDerivAt_id y).const_mul b).const_add (a*x)).pow 2).div_const 2).neg.const_add (a/2) using 1 <;> dsimp [quantum1] <;> ring

theorem quantum2_dx (b c x y : ℝ) :
    HasDerivAt (fun z => quantum2 b c z y) (-b*(b*x+c*y)) x := by
  convert (((((hasDerivAt_id x).const_mul b).add_const (c*y)).pow 2).div_const 2).neg.const_add (c/2) using 1 <;> dsimp [quantum2] <;> ring

def slice (k p t : ℝ) : ℝ := p*t/(1+k^2*t^2)

theorem slice_denominator_positive (k t : ℝ) : 0 < 1+k^2*t^2 := by
  positivity

theorem slice_deriv (k p t : ℝ) :
    HasDerivAt (slice k p) (p*(1-k^2*t^2)/(1+k^2*t^2)^2) t := by
  have hden : 1+k^2*t^2 ≠ 0 := ne_of_gt (slice_denominator_positive k t)
  convert ((hasDerivAt_id t).const_mul p).div
    ((((hasDerivAt_id t).pow 2).const_mul (k^2)).const_add 1) hden using 1 <;>
    dsimp [slice] <;> ring

theorem slice_at_zero (k p : ℝ) : slice k p 0 = 0 := by simp [slice]

theorem slice_deriv_at_zero (k p : ℝ) : deriv (slice k p) 0 = p := by
  rw [(slice_deriv k p 0).deriv]
  simp

theorem cross_acceleration1_from_quantum (a b x y : ℝ) :
    deriv (slice a (b*(a*x+b*y))) 0 = -deriv (fun z => quantum1 a b x z) y := by
  rw [slice_deriv_at_zero, (quantum1_dy a b x y).deriv]
  ring

theorem cross_acceleration2_from_quantum (b c x y : ℝ) :
    deriv (slice c (b*(b*x+c*y))) 0 = -deriv (fun z => quantum2 b c z y) x := by
  rw [slice_deriv_at_zero, (quantum2_dx b c x y).deriv]
  ring

theorem initial_commutator_from_slices (a b c x y : ℝ) :
    (-deriv (slice c (b*(b*x+c*y))) 0,
      deriv (slice a (b*(a*x+b*y))) 0) =
    (-b*(b*x+c*y),b*(a*x+b*y)) := by
  rw [slice_deriv_at_zero,slice_deriv_at_zero]
  congr 1
  ring

theorem weighted_curvature_dx (a b c x y : ℝ) :
    HasDerivAt (fun z => (-b*(b*z+c*y))*(amplitude a b c z y)^2)
      ((-b^2+2*b*(b*x+c*y)*(a*x+b*y))*(amplitude a b c x y)^2) x := by
  convert ((((hasDerivAt_id x).const_mul b).add_const (c*y)).const_mul (-b)).mul
    ((amplitude_dx a b c x y).pow 2) using 1 <;> dsimp only [id_eq,Pi.pow_apply] <;> ring

theorem weighted_curvature_dy (a b c x y : ℝ) :
    HasDerivAt (fun z => (b*(a*x+b*z))*(amplitude a b c x z)^2)
      ((b^2-2*b*(a*x+b*y)*(b*x+c*y))*(amplitude a b c x y)^2) y := by
  convert ((((hasDerivAt_id y).const_mul b).const_add (a*x)).const_mul b).mul
    ((amplitude_dy a b c x y).pow 2) using 1 <;> dsimp only [id_eq,Pi.pow_apply] <;> ring

theorem weighted_curvature_divergence_zero (a b c x y : ℝ) :
    deriv (fun z => (-b*(b*z+c*y))*(amplitude a b c z y)^2) x +
    deriv (fun z => (b*(a*x+b*z))*(amplitude a b c x z)^2) y=0 := by
  rw [(weighted_curvature_dx a b c x y).deriv,
    (weighted_curvature_dy a b c x y).deriv]
  ring

theorem amplitude_separation_iff (a b c : ℝ) :
    (∀ x y : ℝ, amplitude a b c x y*amplitude a b c 0 0 =
      amplitude a b c x 0*amplitude a b c 0 y) ↔ b=0 := by
  constructor
  · intro h
    have h1 := h 1 1
    simp only [amplitude,←Real.exp_add] at h1
    have hexp := Real.exp_injective h1
    dsimp [exponent] at hexp
    nlinarith
  · rintro rfl x y
    simp only [amplitude,←Real.exp_add]
    congr 1
    dsimp [exponent]
    ring

end
end FoundGaussian
