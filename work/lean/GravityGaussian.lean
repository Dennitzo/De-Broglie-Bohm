import Mathlib

/-! Finite real FLRW mode: volume M=a^3, positive physical lapse, and
Gaussian transport. These are real derivative statements on a specified
background; they do not identify a generally covariant stress tensor. -/

namespace GravityGaussian

noncomputable section

def logDensity (mu sigma x : ℝ) : ℝ := -Real.log sigma-(x-mu)^2/(2*sigma^2)
def density (C mu sigma x : ℝ) : ℝ := C*Real.exp (logDensity mu sigma x)
def guidance (N M mu sigma p sd x : ℝ) : ℝ := N*(p/M+sd/sigma*(x-mu))

theorem volume_hasDerivAt (a : ℝ → ℝ) (t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) :
    HasDerivAt (fun u => (a u)^3) (3*H*(a t)^3) t := by
  convert ha.pow 3 using 1 <;> dsimp <;> ring

theorem log_density_time_hasDerivAt (mu sigma : ℝ → ℝ) (t N M p sd x : ℝ)
    (hmu : HasDerivAt mu (N*p/M) t) (hsigma : HasDerivAt sigma (N*sd) t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => logDensity (mu u) (sigma u) x)
      (-N*sd/sigma t+(x-mu t)*(N*p/M)/(sigma t)^2+
        (x-mu t)^2*(N*sd)/(sigma t)^3) t := by
  have hn := ((hasDerivAt_const t x).sub hmu).pow 2
  have hd := (hsigma.pow 2).const_mul 2
  have hdn : 2*(sigma t)^2 ≠ 0 := by positivity
  convert (hsigma.log hs).neg.sub (hn.div hd hdn) using 1 <;>
    dsimp [logDensity] <;> field_simp <;> ring

theorem log_density_space_hasDerivAt (mu sigma x : ℝ) (hs : sigma ≠ 0) :
    HasDerivAt (logDensity mu sigma) (-(x-mu)/sigma^2) x := by
  have h := (hasDerivAt_const x (-Real.log sigma)).sub
    ((((hasDerivAt_id x).sub_const mu).pow 2).div_const (2*sigma^2))
  convert h using 1 <;> dsimp [logDensity] <;> field_simp <;> ring

theorem guidance_space_hasDerivAt (N M mu sigma p sd x : ℝ) :
    HasDerivAt (guidance N M mu sigma p sd) (N*sd/sigma) x := by
  have h := ((((hasDerivAt_id x).sub_const mu).const_mul (sd/sigma)).const_add (p/M)).const_mul N
  convert h using 1 <;> dsimp [guidance] <;> ring

theorem lapse_log_continuity (N M mu sigma p sd x : ℝ)
    (hs : sigma ≠ 0) (hM : M ≠ 0) :
    (-N*sd/sigma+(x-mu)*(N*p/M)/sigma^2+(x-mu)^2*(N*sd)/sigma^3)+
      guidance N M mu sigma p sd x * (-(x-mu)/sigma^2)+N*sd/sigma=0 := by
  unfold guidance
  field_simp
  <;> ring

theorem density_time_hasDerivAt (mu sigma : ℝ → ℝ) (C t N M p sd x : ℝ)
    (hmu : HasDerivAt mu (N*p/M) t) (hsigma : HasDerivAt sigma (N*sd) t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => density C (mu u) (sigma u) x)
      (density C (mu t) (sigma t) x *
        (-N*sd/sigma t+(x-mu t)*(N*p/M)/(sigma t)^2+
          (x-mu t)^2*(N*sd)/(sigma t)^3)) t := by
  have h := (log_density_time_hasDerivAt mu sigma t N M p sd x hmu hsigma hs).exp.const_mul C
  convert h using 1 <;> dsimp [density] <;> ring

theorem density_space_hasDerivAt (C mu sigma x : ℝ) (hs : sigma ≠ 0) :
    HasDerivAt (density C mu sigma) (density C mu sigma x * (-(x-mu)/sigma^2)) x := by
  have h := (log_density_space_hasDerivAt mu sigma x hs).exp.const_mul C
  convert h using 1 <;> dsimp [density] <;> ring

theorem lapse_density_continuity (mu sigma : ℝ → ℝ) (C t N M p sd x : ℝ)
    (hmu : HasDerivAt mu (N*p/M) t) (hsigma : HasDerivAt sigma (N*sd) t)
    (hs : sigma t ≠ 0) (hM : M ≠ 0) :
    deriv (fun u => density C (mu u) (sigma u) x) t+
      deriv (fun y => density C (mu t) (sigma t) y * guidance N M (mu t) (sigma t) p sd y) x=0 := by
  rw [(density_time_hasDerivAt mu sigma C t N M p sd x hmu hsigma hs).deriv]
  have hf : deriv (fun y => density C (mu t) (sigma t) y * guidance N M (mu t) (sigma t) p sd y) x =
      (density C (mu t) (sigma t) x * (-(x-mu t)/(sigma t)^2))*guidance N M (mu t) (sigma t) p sd x+
      density C (mu t) (sigma t) x*(N*sd/sigma t) := by
    exact ((density_space_hasDerivAt C (mu t) (sigma t) x hs).mul
      (guidance_space_hasDerivAt N M (mu t) (sigma t) p sd x)).deriv
  rw [hf]
  have h := lapse_log_continuity N M (mu t) (sigma t) p sd x hs hM
  linear_combination density C (mu t) (sigma t) x * h

theorem centered_affine_guidance (sigma : ℝ → ℝ) (t N sd q0 sigma0 : ℝ)
    (hsigma : HasDerivAt sigma (N*sd) t) (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => q0*sigma u/sigma0)
      (N*(sd/sigma t)*(q0*sigma t/sigma0)) t := by
  convert (hsigma.const_mul q0).div_const sigma0 using 1 <;> dsimp <;> field_simp <;> ring

theorem centered_standardized_coordinate (q sigma : ℝ → ℝ) (t N sd : ℝ)
    (hs : sigma t ≠ 0) (hsigma : HasDerivAt sigma (N*sd) t)
    (hq : HasDerivAt q (N*(sd/sigma t)*q t) t) :
    HasDerivAt (fun u => q u/sigma u) 0 t := by
  convert hq.div hsigma hs using 1 <;> dsimp <;> field_simp <;> ring

theorem lapse_chain_rule (q clock : ℝ → ℝ) (t N v : ℝ)
    (hq : HasDerivAt q v (clock t)) (ht : HasDerivAt clock N t) :
    HasDerivAt (fun u => q (clock u)) (N*v) t := by
  convert hq.comp t ht using 1 <;> ring

theorem lapse_rate_cancellation (N v : ℝ) (hN : 0 < N) : (N*v)/N=v := by
  field_simp

theorem positive_lapse_strictly_monotone (clock N : ℝ → ℝ)
    (hc : ∀ t, HasDerivAt clock (N t) t) (hN : ∀ t, 0 < N t) : StrictMono clock := by
  apply strictMono_of_deriv_pos
  intro t
  rw [(hc t).deriv]
  exact hN t

theorem affine_gaussian_jacobian (C sigma sigma0 q0 : ℝ)
    (hs : sigma ≠ 0) (hs0 : sigma0 ≠ 0) :
    (C/sigma*Real.exp (-((q0*sigma/sigma0)/sigma)^2/2))*(sigma/sigma0)=
      C/sigma0*Real.exp (-(q0/sigma0)^2/2) := by
  have h : (q0*sigma/sigma0)/sigma=q0/sigma0 := by field_simp
  rw [h]
  field_simp

end

end GravityGaussian
