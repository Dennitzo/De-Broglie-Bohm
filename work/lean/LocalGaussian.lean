import Mathlib

/-! A genuine two-mode Gaussian, its derivatives and finite quantum energy.
No convergence or infinite-dimensional measure is asserted. -/
namespace LocalGaussian
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

theorem gaussian_total_cancellation (a b c x y : ℝ) :
    potential a b c x y+quantum1 a b x y+quantum2 b c x y=(a+c)/2 := by
  dsimp [potential,quantum1,quantum2]
  ring

theorem normal_ordered_total_cancellation (a b c x y : ℝ) :
    potential a b c x y+(quantum1 a b x y-a/2)+(quantum2 b c x y-c/2)=0 := by
  dsimp [potential,quantum1,quantum2]
  ring

theorem rational_kernel_positive (x y : ℝ) (h : x ≠ 0 ∨ y ≠ 0) :
    0 < 2*x^2+2*x*y+2*y^2 := by
  rcases h with hx | hy
  · nlinarith [sq_pos_of_ne_zero hx, sq_nonneg y, sq_nonneg (x+y)]
  · nlinarith [sq_pos_of_ne_zero hy, sq_nonneg x, sq_nonneg (x+y)]

theorem constant_energy_guidance (s : ℝ → ℝ) (x sx e t : ℝ)
    (hs : HasDerivAt s sx x) : deriv (fun y => s y-e*t) x=deriv s x := by
  rw [(hs.sub_const (e*t)).deriv,hs.deriv]

theorem constant_subtraction_force (q : ℝ → ℝ) (x qx c : ℝ)
    (hq : HasDerivAt q qx x) : deriv (fun y => q y-c) x=deriv q x := by
  rw [(hq.sub_const c).deriv,hq.deriv]

end
end LocalGaussian
