import Mathlib

/-! Exact affine Gaussian flow for a finite collection of real oscillator
modes. Width dynamics use the Ermakov equation; no claim is made about
an unregularized functional measure on infinitely many field modes. -/

namespace QuantumFieldDynamics

noncomputable section

theorem affine_trajectory_guidance (mu sigma : ℝ → ℝ) (t p sd xi : ℝ)
    (hmu : HasDerivAt mu p t) (hsigma : HasDerivAt sigma sd t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => mu u + sigma u * xi)
      (p + sd/sigma t*((mu t+sigma t*xi)-mu t)) t := by
  convert hmu.add (hsigma.mul_const xi) using 1 <;> field_simp <;> ring

theorem standardized_coordinate_hasDerivAt_zero (q mu sigma : ℝ → ℝ)
    (t p sd : ℝ) (hs : sigma t ≠ 0)
    (hq : HasDerivAt q (p+sd/sigma t*(q t-mu t)) t)
    (hmu : HasDerivAt mu p t) (hsigma : HasDerivAt sigma sd t) :
    HasDerivAt (fun u => (q u-mu u)/sigma u) 0 t := by
  convert (hq.sub hmu).div hsigma hs using 1 <;> dsimp <;> field_simp <;> ring

theorem squeezing_ratio_riccati_hasDerivAt (sigma sd : ℝ → ℝ) (t omega : ℝ)
    (hs : sigma t ≠ 0) (hsigma : HasDerivAt sigma (sd t) t)
    (hsd : HasDerivAt sd (1/(4*(sigma t)^3)-omega^2*sigma t) t) :
    HasDerivAt (fun u => sd u/sigma u)
      (1/(4*(sigma t)^4)-omega^2-(sd t/sigma t)^2) t := by
  convert hsd.div hsigma hs using 1 <;> dsimp <;> field_simp <;> ring

theorem ermakov_effective_acceleration (omega mu sigma xi mu2 sigma2 : ℝ)
    (hs : sigma ≠ 0) (hmu : mu2 = -omega^2*mu)
    (he : sigma2+omega^2*sigma=1/(4*sigma^3)) :
    mu2+sigma2*xi = -omega^2*(mu+sigma*xi)+((mu+sigma*xi)-mu)/(4*sigma^4) := by
  have hsig : sigma2=1/(4*sigma^3)-omega^2*sigma := by linarith
  rw [hmu, hsig]
  field_simp
  <;> ring

theorem affine_trajectory_acceleration (mu sigma p sd : ℝ → ℝ)
    (t omega xi : ℝ)
    (hp : HasDerivAt p (-omega^2*mu t) t)
    (hsd : HasDerivAt sd (1/(4*(sigma t)^3)-omega^2*sigma t) t)
    (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => p u+sd u*xi)
      (-omega^2*(mu t+sigma t*xi)+((mu t+sigma t*xi)-mu t)/(4*(sigma t)^4)) t := by
  convert hp.add (hsd.mul_const xi) using 1 <;> field_simp <;> ring

theorem vacuum_quantum_force_cancellation (omega sigma q : ℝ)
    (hs : sigma ≠ 0) (hwidth : sigma^2=1/(2*omega)) (hw : omega ≠ 0) :
    -omega^2*q+q/(4*sigma^4)=0 := by
  have hp : 4*sigma^4*omega^2=1 := by
    calc
      _ = (2*omega*sigma^2)^2 := by ring
      _ = 1 := by rw [hwidth]; field_simp
  field_simp
  nlinarith [congrArg (fun z : ℝ => z*q) hp]

theorem vacuum_width_solves_ermakov (omega sigma : ℝ)
    (hs : sigma ≠ 0) (hwidth : sigma^2=1/(2*omega)) (hw : omega ≠ 0) :
    omega^2*sigma=1/(4*sigma^3) := by
  have hp : 4*sigma^4*omega^2=1 := by
    calc
      _ = (2*omega*sigma^2)^2 := by ring
      _ = 1 := by rw [hwidth]; field_simp
  field_simp
  nlinarith [hp]

def affineDensity (C mu sigma q : ℝ) : ℝ :=
  C/sigma*Real.exp (-((q-mu)/sigma)^2/2)

theorem affine_density_scaled (C mu sigma xi : ℝ) (hs : sigma ≠ 0) :
    affineDensity C mu sigma (mu+sigma*xi)*sigma = C*Real.exp (-xi^2/2) := by
  unfold affineDensity
  have harg : (mu+sigma*xi-mu)/sigma=xi := by field_simp; ring
  rw [harg]
  field_simp

theorem affine_density_jacobian (C mu mu0 sigma sigma0 xi : ℝ)
    (hs : sigma ≠ 0) (hs0 : sigma0 ≠ 0) :
    affineDensity C mu sigma (mu+sigma*xi)*(sigma/sigma0) =
      affineDensity C mu0 sigma0 (mu0+sigma0*xi) := by
  have h := affine_density_scaled C mu sigma xi hs
  have h0 := affine_density_scaled C mu0 sigma0 xi hs0
  apply (mul_right_cancel₀ hs0)
  calc
    _ = affineDensity C mu sigma (mu+sigma*xi)*sigma := by field_simp
    _ = _ := h.trans h0.symm

theorem gaussian_width_square_positive (omega sigma0 c s : ℝ)
    (hw : 0 < omega) (hs : 0 < sigma0) (hcs : c^2+s^2=1) :
    0 < sigma0^2*c^2+s^2/(4*omega^2*sigma0^2) := by
  have hc0 : 0 ≤ c^2 := sq_nonneg c
  have hs0 : 0 ≤ s^2 := sq_nonneg s
  have hden : 0 < 4*omega^2*sigma0^2 := by positivity
  have ha : 0 ≤ sigma0^2*c^2 := by positivity
  have hb : 0 ≤ s^2/(4*omega^2*sigma0^2) := by positivity
  by_cases hc : c=0
  · have hs1 : s^2=1 := by rw [hc] at hcs; nlinarith
    rw [hc, hs1]
    positivity
  · have hapos : 0 < sigma0^2*c^2 := by positivity
    linarith

theorem squeezed_width_square_positive (omega sigma0 t : ℝ)
    (hw : 0 < omega) (hs : 0 < sigma0) :
    0 < sigma0^2*(Real.cos (omega*t))^2+
      (Real.sin (omega*t))^2/(4*omega^2*sigma0^2) := by
  apply gaussian_width_square_positive omega sigma0 _ _ hw hs
  nlinarith [Real.sin_sq_add_cos_sq (omega*t)]

theorem squeezed_width_positive (omega sigma0 t : ℝ)
    (hw : 0 < omega) (hs : 0 < sigma0) :
    0 < Real.sqrt (sigma0^2*(Real.cos (omega*t))^2+
      (Real.sin (omega*t))^2/(4*omega^2*sigma0^2)) := by
  exact Real.sqrt_pos.mpr (squeezed_width_square_positive omega sigma0 t hw hs)

theorem oscillator_energy_nonnegative (omega q p : ℝ) :
    0 ≤ (p^2+omega^2*q^2)/2 := by positivity

theorem finite_mode_energy_nonnegative {ι : Type*} [Fintype ι]
    (omega q p : ι → ℝ) :
    0 ≤ ∑ k : ι, ((p k)^2+(omega k)^2*(q k)^2)/2 := by
  apply Finset.sum_nonneg
  intro k hk
  exact oscillator_energy_nonnegative _ _ _

def gaussianEnergy (omega mu sigma p sd : ℝ) : ℝ :=
  (p^2+sd^2+omega^2*(mu^2+sigma^2))/2+1/(8*sigma^2)

theorem gaussian_energy_positive (omega mu sigma p sd : ℝ) (hs : sigma ≠ 0) :
    0 < gaussianEnergy omega mu sigma p sd := by
  unfold gaussianEnergy
  positivity

theorem gaussian_energy_hasDerivAt_zero (mu sigma p sd : ℝ → ℝ) (t omega : ℝ)
    (hs : sigma t ≠ 0)
    (hmu : HasDerivAt mu (p t) t) (hsigma : HasDerivAt sigma (sd t) t)
    (hp : HasDerivAt p (-omega^2*mu t) t)
    (hsd : HasDerivAt sd (1/(4*(sigma t)^3)-omega^2*sigma t) t) :
    HasDerivAt (fun u => gaussianEnergy omega (mu u) (sigma u) (p u) (sd u)) 0 t := by
  have hk := (hp.pow 2).add (hsd.pow 2)
  have hv := ((hmu.pow 2).add (hsigma.pow 2)).const_mul (omega^2)
  have hm := (hk.add hv).div_const 2
  have hden := (hsigma.pow 2).const_mul 8
  have hdenne : 8*(sigma t)^2 ≠ 0 := by positivity
  have hq := (hasDerivAt_const t 1).div hden hdenne
  convert hm.add hq using 1 <;> dsimp [gaussianEnergy] <;> field_simp <;> ring

theorem vacuum_gaussian_energy (omega sigma : ℝ) (hw : omega ≠ 0)
    (hwidth : sigma^2=1/(2*omega)) :
    gaussianEnergy omega 0 sigma 0 0 = omega/2 := by
  unfold gaussianEnergy
  rw [hwidth]
  field_simp
  <;> ring

theorem finite_gaussian_energy_nonnegative {ι : Type*} [Fintype ι]
    (omega mu sigma p sd : ι → ℝ) (hs : ∀ k, sigma k ≠ 0) :
    0 ≤ ∑ k : ι, gaussianEnergy (omega k) (mu k) (sigma k) (p k) (sd k) := by
  apply Finset.sum_nonneg
  intro k hk
  exact le_of_lt (gaussian_energy_positive _ _ _ _ _ (hs k))

end

end QuantumFieldDynamics
