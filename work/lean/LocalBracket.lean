import Mathlib

/-! Exact finite Poisson correction of a nonlocal Gaussian quantum potential.
The two lapse weights select different sites. The classical spatial-gradient
Hamiltonian is omitted: this is its additive quantum correction. -/
namespace LocalBracket
noncomputable section

def q1 (a b x y : ℝ) := a/2-(a*x+b*y)^2/2
def q2 (b c x y : ℝ) := c/2-(b*x+c*y)^2/2
def h1 (a b x y p r : ℝ) := p^2/2+q1 a b x y
def h2 (b c x y p r : ℝ) := r^2/2+q2 b c x y

def poisson (f g : ℝ → ℝ → ℝ → ℝ → ℝ) (x y p r : ℝ) :=
  deriv (fun z => f z y p r) x*deriv (fun z => g x y z r) p+
  deriv (fun z => f x z p r) y*deriv (fun z => g x y p z) r-
  deriv (fun z => f x y z r) p*deriv (fun z => g z y p r) x-
  deriv (fun z => f x y p z) r*deriv (fun z => g x z p r) y

theorem q1_dy (a b x y : ℝ) :
    HasDerivAt (fun z => q1 a b x z) (-b*(a*x+b*y)) y := by
  convert (((((hasDerivAt_id y).const_mul b).const_add (a*x)).pow 2).div_const 2).const_sub
    (a/2) using 1 <;> dsimp [q1] <;> ring

theorem q2_dx (b c x y : ℝ) :
    HasDerivAt (fun z => q2 b c z y) (-b*(b*x+c*y)) x := by
  convert (((((hasDerivAt_id x).const_mul b).add_const (c*y)).pow 2).div_const 2).const_sub
    (c/2) using 1 <;> dsimp [q2] <;> ring

theorem h1_dy (a b x y p r : ℝ) :
    HasDerivAt (fun z => h1 a b x z p r) (-b*(a*x+b*y)) y := by
  exact (q1_dy a b x y).const_add (p^2/2)

theorem h2_dx (b c x y p r : ℝ) :
    HasDerivAt (fun z => h2 b c z y p r) (-b*(b*x+c*y)) x := by
  exact (q2_dx b c x y).const_add (r^2/2)

theorem h1_dp (a b x y p r : ℝ) :
    HasDerivAt (fun z => h1 a b x y z r) p p := by
  convert (((hasDerivAt_id p).pow 2).div_const 2).add_const (q1 a b x y) using 1 <;>
    dsimp [h1] <;> ring

theorem h2_dr (b c x y p r : ℝ) :
    HasDerivAt (fun z => h2 b c x y p z) r r := by
  convert (((hasDerivAt_id r).pow 2).div_const 2).add_const (q2 b c x y) using 1 <;>
    dsimp [h2] <;> ring

theorem gaussian_quantum_poisson (a b c x y p r : ℝ) :
    poisson (h1 a b) (h2 b c) x y p r=
      b*(p*(b*x+c*y)-r*(a*x+b*y)) := by
  unfold poisson
  have hconst1 : deriv (fun z => h2 b c x y z r) p=0 := by
    exact (hasDerivAt_const p (r^2/2+q2 b c x y)).deriv
  have hconst2 : deriv (fun z => h1 a b x y p z) r=0 := by
    exact (hasDerivAt_const r (p^2/2+q1 a b x y)).deriv
  rw [hconst1,hconst2,(h1_dy a b x y p r).deriv,
    (h2_dx b c x y p r).deriv,(h1_dp a b x y p r).deriv,
    (h2_dr b c x y p r).deriv]
  ring

theorem rational_nonlocal_witness : poisson (h1 2 1) (h2 1 2) 1 0 1 0=1 := by
  rw [gaussian_quantum_poisson]
  norm_num

theorem rational_nonlocal_nonzero : poisson (h1 2 1) (h2 1 2) 1 0 1 0 ≠ 0 := by
  rw [rational_nonlocal_witness]
  norm_num

theorem gaussian_guidance_weak_zero (a b c x y : ℝ) :
    poisson (h1 a b) (h2 b c) x y 0 0=0 := by
  rw [gaussian_quantum_poisson]
  ring

theorem local_kernel_zero (a c x y p r : ℝ) :
    poisson (h1 a 0) (h2 0 c) x y p r=0 := by
  rw [gaussian_quantum_poisson]
  ring

theorem ultralocal_configuration_zero (f g : ℝ → ℝ) (x y p r : ℝ) :
    poisson (fun x _ p _ => p^2/2+f x) (fun _ y _ r => r^2/2+g y) x y p r=0 := by
  unfold poisson
  rw [(hasDerivAt_const p (r^2/2+g y)).deriv,
    (hasDerivAt_const y (p^2/2+f x)).deriv,
    (hasDerivAt_const x (r^2/2+g y)).deriv,
    (hasDerivAt_const r (p^2/2+f x)).deriv]
  ring

theorem lapse_antisymmetric_factor (n1 n2 m1 m2 u v : ℝ) :
    n1*m2*u+n2*m1*v-m1*n2*u-m2*n1*v=
      (n1*m2-n2*m1)*(u-v) := by ring

theorem subtract_constants_poisson (f g : ℝ → ℝ → ℝ → ℝ → ℝ)
    (x y p r cf cg : ℝ)
    (hfx : DifferentiableAt ℝ (fun z => f z y p r) x)
    (hfy : DifferentiableAt ℝ (fun z => f x z p r) y)
    (hfp : DifferentiableAt ℝ (fun z => f x y z r) p)
    (hfr : DifferentiableAt ℝ (fun z => f x y p z) r)
    (hgx : DifferentiableAt ℝ (fun z => g z y p r) x)
    (hgy : DifferentiableAt ℝ (fun z => g x z p r) y)
    (hgp : DifferentiableAt ℝ (fun z => g x y z r) p)
    (hgr : DifferentiableAt ℝ (fun z => g x y p z) r) :
    poisson (fun x y p r => f x y p r-cf) (fun x y p r => g x y p r-cg) x y p r=
      poisson f g x y p r := by
  unfold poisson
  rw [(hfx.hasDerivAt.sub_const cf).deriv,(hfy.hasDerivAt.sub_const cf).deriv,
    (hfp.hasDerivAt.sub_const cf).deriv,(hfr.hasDerivAt.sub_const cf).deriv,
    (hgx.hasDerivAt.sub_const cg).deriv,(hgy.hasDerivAt.sub_const cg).deriv,
    (hgp.hasDerivAt.sub_const cg).deriv,(hgr.hasDerivAt.sub_const cg).deriv]

end
end LocalBracket
