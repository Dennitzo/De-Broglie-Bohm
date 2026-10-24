import Mathlib

/-! One real normal-mode coordinate of a free scalar field. The oscillator
mass is one and hbar is one. All calculus claims below are actual real
derivative statements. Functional infinite-dimensional limits are outside
their scope. -/

namespace QuantumFieldGaussian

noncomputable section

def gaussian (C a mu x : ℝ) : ℝ := C * Real.exp (-a*(x-mu)^2)

theorem gaussian_positive (C a mu x : ℝ) (hC : 0 < C) :
    0 < gaussian C a mu x := by
  unfold gaussian
  exact mul_pos hC (Real.exp_pos _)

theorem gaussian_hasDerivAt (C a mu x : ℝ) :
    HasDerivAt (gaussian C a mu) (-2*a*(x-mu)*gaussian C a mu x) x := by
  have hi := (((hasDerivAt_id x).sub_const mu).pow 2).const_mul (-a)
  convert hi.exp.const_mul C using 1 <;> dsimp [gaussian] <;> ring

theorem gaussian_deriv (C a mu x : ℝ) :
    deriv (gaussian C a mu) x = -2*a*(x-mu)*gaussian C a mu x :=
  (gaussian_hasDerivAt C a mu x).deriv

theorem gaussian_first_derivative_hasDerivAt (C a mu x : ℝ) :
    HasDerivAt (fun y => -2*a*(y-mu)*gaussian C a mu y)
      ((4*a^2*(x-mu)^2-2*a)*gaussian C a mu x) x := by
  have hl := ((hasDerivAt_id x).sub_const mu).const_mul (-2*a)
  convert hl.mul (gaussian_hasDerivAt C a mu x) using 1 <;> dsimp <;> ring

theorem gaussian_second_deriv (C a mu x : ℝ) :
    deriv (fun y => deriv (gaussian C a mu) y) x =
      (4*a^2*(x-mu)^2-2*a)*gaussian C a mu x := by
  simp_rw [gaussian_deriv]
  exact (gaussian_first_derivative_hasDerivAt C a mu x).deriv

def quantumPotential (C mu sigma x : ℝ) : ℝ :=
  -(1/2:ℝ) * deriv (fun y => deriv (gaussian C (1/(4*sigma^2)) mu) y) x /
    gaussian C (1/(4*sigma^2)) mu x

theorem quantum_potential_formula (C mu sigma x : ℝ)
    (hC : 0 < C) (hs : sigma ≠ 0) :
    quantumPotential C mu sigma x = 1/(4*sigma^2)-(x-mu)^2/(8*sigma^4) := by
  have hR := ne_of_gt (gaussian_positive C (1/(4*sigma^2)) mu x hC)
  unfold quantumPotential
  rw [gaussian_second_deriv]
  field_simp
  <;> ring

theorem quantum_potential_hasDerivAt (C mu sigma x : ℝ)
    (hC : 0 < C) (hs : sigma ≠ 0) :
    HasDerivAt (quantumPotential C mu sigma) (-(x-mu)/(4*sigma^4)) x := by
  have hf : quantumPotential C mu sigma =
      fun y => 1/(4*sigma^2)-(y-mu)^2/(8*sigma^4) := by
    funext y
    exact quantum_potential_formula C mu sigma y hC hs
  rw [hf]
  have h := (hasDerivAt_const x (1/(4*sigma^2))).sub
    ((((hasDerivAt_id x).sub_const mu).pow 2).div_const (8*sigma^4))
  convert h using 1 <;> dsimp <;> field_simp <;> ring

theorem quantum_force_formula (C mu sigma x : ℝ)
    (hC : 0 < C) (hs : sigma ≠ 0) :
    -deriv (quantumPotential C mu sigma) x = (x-mu)/(4*sigma^4) := by
  rw [(quantum_potential_hasDerivAt C mu sigma x hC hs).deriv]
  ring

def logDensity (mu sigma x : ℝ) : ℝ :=
  -Real.log sigma - (x-mu)^2/(2*sigma^2)

def density (C mu sigma x : ℝ) : ℝ := C * Real.exp (logDensity mu sigma x)

def guidance (mu sigma p sd x : ℝ) : ℝ := p + sd/sigma*(x-mu)

theorem log_density_time_hasDerivAt (mu sigma : ℝ → ℝ) (t p sd x : ℝ)
    (hmu : HasDerivAt mu p t) (hsigma : HasDerivAt sigma sd t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => logDensity (mu u) (sigma u) x)
      (-sd/sigma t + (x-mu t)*p/(sigma t)^2 + (x-mu t)^2*sd/(sigma t)^3) t := by
  have hnum := ((hasDerivAt_const t x).sub hmu).pow 2
  have hden := (hsigma.pow 2).const_mul 2
  have hdenne : 2*(sigma t)^2 ≠ 0 := by positivity
  have h := (hsigma.log hs).neg.sub (hnum.div hden hdenne)
  convert h using 1 <;> dsimp [logDensity] <;> field_simp <;> ring

theorem log_density_space_hasDerivAt (mu sigma x : ℝ) (hs : sigma ≠ 0) :
    HasDerivAt (logDensity mu sigma) (-(x-mu)/sigma^2) x := by
  have h := (hasDerivAt_const x (-Real.log sigma)).sub
    ((((hasDerivAt_id x).sub_const mu).pow 2).div_const (2*sigma^2))
  convert h using 1 <;> dsimp [logDensity] <;> field_simp <;> ring

theorem guidance_space_hasDerivAt (mu sigma p sd x : ℝ) :
    HasDerivAt (guidance mu sigma p sd) (sd/sigma) x := by
  have h := (((hasDerivAt_id x).sub_const mu).const_mul (sd/sigma)).const_add p
  convert h using 1 <;> dsimp [guidance] <;> ring

theorem gaussian_log_continuity (mu sigma p sd x : ℝ) (hs : sigma ≠ 0) :
    (-sd/sigma+(x-mu)*p/sigma^2+(x-mu)^2*sd/sigma^3) +
      guidance mu sigma p sd x * (-(x-mu)/sigma^2) + sd/sigma = 0 := by
  unfold guidance
  field_simp
  <;> ring

theorem density_time_hasDerivAt (mu sigma : ℝ → ℝ) (C t p sd x : ℝ)
    (hmu : HasDerivAt mu p t) (hsigma : HasDerivAt sigma sd t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => density C (mu u) (sigma u) x)
      (density C (mu t) (sigma t) x *
        (-sd/sigma t+(x-mu t)*p/(sigma t)^2+(x-mu t)^2*sd/(sigma t)^3)) t := by
  have h := (log_density_time_hasDerivAt mu sigma t p sd x hmu hsigma hs).exp.const_mul C
  convert h using 1 <;> dsimp [density] <;> ring

theorem density_space_hasDerivAt (C mu sigma x : ℝ) (hs : sigma ≠ 0) :
    HasDerivAt (density C mu sigma) (density C mu sigma x * (-(x-mu)/sigma^2)) x := by
  have h := (log_density_space_hasDerivAt mu sigma x hs).exp.const_mul C
  convert h using 1 <;> dsimp [density] <;> ring

theorem density_continuity (mu sigma : ℝ → ℝ) (C t p sd x : ℝ)
    (hmu : HasDerivAt mu p t) (hsigma : HasDerivAt sigma sd t)
    (hs : sigma t ≠ 0) :
    deriv (fun u => density C (mu u) (sigma u) x) t +
      deriv (fun y => density C (mu t) (sigma t) y * guidance (mu t) (sigma t) p sd y) x = 0 := by
  rw [(density_time_hasDerivAt mu sigma C t p sd x hmu hsigma hs).deriv]
  have hflux :
      deriv (fun y => density C (mu t) (sigma t) y * guidance (mu t) (sigma t) p sd y) x =
      (density C (mu t) (sigma t) x * (-(x-mu t)/(sigma t)^2)) *
          guidance (mu t) (sigma t) p sd x + density C (mu t) (sigma t) x * (sd/sigma t) := by
    exact ((density_space_hasDerivAt C (mu t) (sigma t) x hs).mul
      (guidance_space_hasDerivAt (mu t) (sigma t) p sd x)).deriv
  rw [hflux]
  have h := gaussian_log_continuity (mu t) (sigma t) p sd x hs
  linear_combination density C (mu t) (sigma t) x * h

theorem density_explicit_formula (C mu sigma x : ℝ) (hs : 0 < sigma) :
    density C mu sigma x = C/sigma*Real.exp (-((x-mu)/sigma)^2/2) := by
  have hs0 : sigma ≠ 0 := ne_of_gt hs
  have he : -((x-mu)^2/(2*sigma^2)) = -((x-mu)/sigma)^2/2 := by
    field_simp
  unfold density logDensity
  rw [sub_eq_add_neg, Real.exp_add, Real.exp_neg, Real.exp_log hs, he]
  ring

def phase (mu p b theta x : ℝ) : ℝ :=
  p*(x-mu)+(b/2)*(x-mu)^2+theta

theorem phase_space_hasDerivAt (mu p b theta x : ℝ) :
    HasDerivAt (phase mu p b theta) (p+b*(x-mu)) x := by
  have hy := (hasDerivAt_id x).sub_const mu
  have h := ((hy.const_mul p).add ((hy.pow 2).const_mul (b/2))).add_const theta
  convert h using 1 <;> dsimp [phase] <;> ring

theorem phase_time_hasDerivAt (mu p b theta : ℝ → ℝ) (t pd bd td x : ℝ)
    (hmu : HasDerivAt mu (p t) t) (hp : HasDerivAt p pd t)
    (hb : HasDerivAt b bd t) (htheta : HasDerivAt theta td t) :
    HasDerivAt (fun u => phase (mu u) (p u) (b u) (theta u) x)
      (pd*(x-mu t)-(p t)^2+(bd/2)*(x-mu t)^2-b t*p t*(x-mu t)+td) t := by
  have hy := (hasDerivAt_const t x).sub hmu
  have h := ((hp.mul hy).add ((hb.div_const 2).mul (hy.pow 2))).add htheta
  convert h using 1 <;> dsimp [phase] <;> ring

theorem phase_hamilton_jacobi (mu p b theta : ℝ → ℝ) (C t omega sigma x : ℝ)
    (hC : 0 < C) (hs : sigma ≠ 0)
    (hmu : HasDerivAt mu (p t) t)
    (hp : HasDerivAt p (-omega^2*mu t) t)
    (hb : HasDerivAt b (1/(4*sigma^4)-omega^2-(b t)^2) t)
    (htheta : HasDerivAt theta (((p t)^2-omega^2*(mu t)^2)/2-1/(4*sigma^2)) t) :
    deriv (fun u => phase (mu u) (p u) (b u) (theta u) x) t +
      (deriv (phase (mu t) (p t) (b t) (theta t)) x)^2/2 +
      omega^2*x^2/2+quantumPotential C (mu t) sigma x = 0 := by
  rw [(phase_time_hasDerivAt mu p b theta t _ _ _ x hmu hp hb htheta).deriv]
  rw [(phase_space_hasDerivAt (mu t) (p t) (b t) (theta t) x).deriv]
  rw [quantum_potential_formula C (mu t) sigma x hC hs]
  field_simp
  <;> ring

end

end QuantumFieldGaussian
