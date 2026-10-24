import Mathlib

/-! Finite positive-weight Gaussian modes coupled to flat FLRW through
their complete operator expectation source. All evolution equations used
as premises are displayed. No existence theorem for the coupled ODE or
unregularized continuum field is asserted by these statements. -/

namespace BackreactionCoupled

noncomputable section

def modeDensity (k a s v : ℝ) : ℝ :=
  (v^2+(k^2/a^2)*s^2+1/(4*a^6*s^2))/2

def modePressure (k a s v : ℝ) : ℝ :=
  (v^2+1/(4*a^6*s^2)-(k^2/a^2)*s^2/3)/2

def totalDensity {ι : Type*} [Fintype ι] (w k s v : ι → ℝ) (a : ℝ) : ℝ :=
  ∑ i : ι, w i*modeDensity (k i) a (s i) (v i)

def totalPressure {ι : Type*} [Fintype ι] (w k s v : ι → ℝ) (a : ℝ) : ℝ :=
  ∑ i : ι, w i*modePressure (k i) a (s i) (v i)

def spectralWeight {ι : Type*} [Fintype ι] (w k : ι → ℝ) : ℝ :=
  ∑ i : ι, w i*k i/2

theorem completed_square_energy (k a s v : ℝ) (ha : a ≠ 0) (hs : s ≠ 0) :
    modeDensity k a s v-k/(2*a^4)=v^2/2+(k*s/a-1/(2*a^3*s))^2/2 := by
  unfold modeDensity
  field_simp
  <;> ring

theorem mode_uncertainty_lower_bound (k a s v : ℝ) (ha : a ≠ 0) (hs : s ≠ 0) :
    k/(2*a^4) ≤ modeDensity k a s v := by
  have h := completed_square_energy k a s v ha hs
  have h1 := sq_nonneg v
  have h2 := sq_nonneg (k*s/a-1/(2*a^3*s))
  linarith

theorem mode_density_positive (k a s v : ℝ) (hk : 0 < k) (ha : 0 < a) (hs : 0 < s) :
    0 < modeDensity k a s v := by
  have h := mode_uncertainty_lower_bound k a s v (ne_of_gt ha) (ne_of_gt hs)
  have hl : 0 < k/(2*a^4) := by positivity
  linarith

theorem mode_pressure_bounds (k a s v : ℝ) :
    -modeDensity k a s v/3 ≤ modePressure k a s v ∧
      modePressure k a s v ≤ modeDensity k a s v := by
  have h1 : 0 ≤ v^2 := sq_nonneg v
  have h2 : 0 ≤ (k^2/a^2)*s^2 := by positivity
  have h3 : 0 ≤ 1/(4*a^6*s^2) := by positivity
  unfold modeDensity modePressure
  constructor <;> nlinarith

theorem mode_null_energy_positive (k a s v : ℝ) (hk : 0 < k) (ha : 0 < a) (hs : 0 < s) :
    0 < modeDensity k a s v+modePressure k a s v := by
  have hr := mode_density_positive k a s v hk ha hs
  have hp := (mode_pressure_bounds k a s v).1
  linarith

theorem spectral_weight_positive {ι : Type*} [Fintype ι] [Nonempty ι]
    (w k : ι → ℝ) (hw : ∀ i, 0 < w i) (hk : ∀ i, 0 < k i) :
    0 < spectralWeight w k := by
  unfold spectralWeight
  apply Finset.sum_pos
  · intro i hi
    exact div_pos (mul_pos (hw i) (hk i)) (by norm_num)
  · exact Finset.univ_nonempty

theorem total_uncertainty_lower_bound {ι : Type*} [Fintype ι]
    (w k s v : ι → ℝ) (a : ℝ) (hw : ∀ i, 0 ≤ w i)
    (ha : a ≠ 0) (hs : ∀ i, s i ≠ 0) :
    spectralWeight w k/a^4 ≤ totalDensity w k s v a := by
  calc
    _ = ∑ i : ι, w i*((k i)/(2*a^4)) := by
      unfold spectralWeight
      rw [Finset.sum_div]
      apply Finset.sum_congr rfl
      intro i hi
      ring
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro i hi
      exact mul_le_mul_of_nonneg_left (mode_uncertainty_lower_bound _ _ _ _ ha (hs i)) (hw i)

theorem total_density_positive {ι : Type*} [Fintype ι] [Nonempty ι]
    (w k s v : ι → ℝ) (a : ℝ) (hw : ∀ i, 0 < w i) (hk : ∀ i, 0 < k i)
    (ha : 0 < a) (hs : ∀ i, 0 < s i) :
    0 < totalDensity w k s v a := by
  have h := total_uncertainty_lower_bound w k s v a (fun i => le_of_lt (hw i))
    (ne_of_gt ha) (fun i => ne_of_gt (hs i))
  have hL := spectral_weight_positive w k hw hk
  have hp : 0 < spectralWeight w k/a^4 := by positivity
  linarith

theorem total_pressure_bounds {ι : Type*} [Fintype ι]
    (w k s v : ι → ℝ) (a : ℝ) (hw : ∀ i, 0 ≤ w i) :
    -totalDensity w k s v a/3 ≤ totalPressure w k s v a ∧
      totalPressure w k s v a ≤ totalDensity w k s v a := by
  have h1 : 0 ≤ ∑ i : ι, w i*(modeDensity (k i) a (s i) (v i)-modePressure (k i) a (s i) (v i)) := by
    apply Finset.sum_nonneg
    intro i hi
    exact mul_nonneg (hw i) (sub_nonneg.mpr (mode_pressure_bounds _ _ _ _).2)
  have h2 : 0 ≤ ∑ i : ι, w i*(modeDensity (k i) a (s i) (v i)+3*modePressure (k i) a (s i) (v i)) := by
    apply Finset.sum_nonneg
    intro i hi
    have hp := (mode_pressure_bounds (k i) a (s i) (v i)).1
    apply mul_nonneg (hw i)
    linarith
  have he1 : (∑ i : ι, w i*(modeDensity (k i) a (s i) (v i)-modePressure (k i) a (s i) (v i)))=
      totalDensity w k s v a-totalPressure w k s v a := by
    simp only [mul_sub, Finset.sum_sub_distrib, totalDensity, totalPressure]
  have he2 : (∑ i : ι, w i*(modeDensity (k i) a (s i) (v i)+3*modePressure (k i) a (s i) (v i)))=
      totalDensity w k s v a+3*totalPressure w k s v a := by
    simp only [mul_add, Finset.sum_add_distrib, totalDensity, totalPressure]
    congr 1
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i hi
    ring
  rw [he1] at h1
  rw [he2] at h2
  constructor <;> linarith

theorem mode_density_hasDerivAt (a s v : ℝ → ℝ) (k t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hs : HasDerivAt s (v t) t)
    (hv : HasDerivAt v (-3*H*v t-(k^2/(a t)^2)*s t+1/(4*(a t)^6*(s t)^3)) t)
    (ha0 : a t ≠ 0) (hs0 : s t ≠ 0) :
    HasDerivAt (fun u => modeDensity k (a u) (s u) (v u))
      (-3*H*(modeDensity k (a t) (s t) (v t)+modePressure k (a t) (s t) (v t))) t := by
  have ha2 : (a t)^2 ≠ 0 := by positivity
  have hW := (hasDerivAt_const t (k^2)).div (ha.pow 2) ha2
  have hd := ((ha.pow 6).const_mul 4).mul (hs.pow 2)
  have hdn : 4*(a t)^6*(s t)^2 ≠ 0 := by positivity
  have hL := (hasDerivAt_const t 1).div hd hdn
  have h := (((hv.pow 2).add (hW.mul (hs.pow 2))).add hL).div_const 2
  convert h using 1 <;> dsimp [modeDensity, modePressure] <;> field_simp <;> ring

theorem mode_source_conserved (a s v : ℝ → ℝ) (k t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hs : HasDerivAt s (v t) t)
    (hv : HasDerivAt v (-3*H*v t-(k^2/(a t)^2)*s t+1/(4*(a t)^6*(s t)^3)) t)
    (ha0 : a t ≠ 0) (hs0 : s t ≠ 0) :
    deriv (fun u => modeDensity k (a u) (s u) (v u)) t+
      3*H*(modeDensity k (a t) (s t) (v t)+modePressure k (a t) (s t) (v t))=0 := by
  rw [(mode_density_hasDerivAt a s v k t H ha hs hv ha0 hs0).deriv]
  ring

theorem total_density_hasDerivAt {ι : Type*} [Fintype ι]
    (w k : ι → ℝ) (a : ℝ → ℝ) (s v : ι → ℝ → ℝ) (t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) (ha0 : a t ≠ 0)
    (hs : ∀ i, HasDerivAt (s i) (v i t) t) (hs0 : ∀ i, s i t ≠ 0)
    (hv : ∀ i, HasDerivAt (v i)
      (-3*H*v i t+(-(k i)^2/(a t)^2)*s i t+1/(4*(a t)^6*(s i t)^3)) t) :
    HasDerivAt (fun u => totalDensity w k (fun i => s i u) (fun i => v i u) (a u))
      (-3*H*(totalDensity w k (fun i => s i t) (fun i => v i t) (a t)+
        totalPressure w k (fun i => s i t) (fun i => v i t) (a t))) t := by
  have hi : ∀ i ∈ (Finset.univ : Finset ι),
      HasDerivAt (fun u => w i*modeDensity (k i) (a u) (s i u) (v i u))
        (w i*(-3*H*(modeDensity (k i) (a t) (s i t) (v i t)+modePressure (k i) (a t) (s i t) (v i t)))) t := by
    intro i _
    have hvi : HasDerivAt (v i)
        (-3*H*v i t-((k i)^2/(a t)^2)*s i t+1/(4*(a t)^6*(s i t)^3)) t := by
      convert hv i using 1 <;> ring
    exact (mode_density_hasDerivAt a (s i) (v i) (k i) t H ha (hs i) hvi ha0 (hs0 i)).const_mul (w i)
  have hh := HasDerivAt.fun_sum hi
  have he : (∑ i : ι, w i*(-3*H*(modeDensity (k i) (a t) (s i t) (v i t)+modePressure (k i) (a t) (s i t) (v i t))))=
      -3*H*(totalDensity w k (fun i => s i t) (fun i => v i t) (a t)+
        totalPressure w k (fun i => s i t) (fun i => v i t) (a t)) := by
    unfold totalDensity totalPressure
    rw [← Finset.sum_add_distrib, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i hi
    ring
  rw [he] at hh
  exact hh

theorem density_decay_rate (H rho pressure : ℝ) (hH : 0 ≤ H) (hp : -rho/3 ≤ pressure) :
    -3*H*(rho+pressure) ≤ -2*H*rho := by
  have h := mul_nonneg hH (by linarith : 0 ≤ rho+3*pressure)
  nlinarith

theorem weighted_density_scale_hasDerivAt (a rho : ℝ → ℝ) (t H pressure : ℝ)
    (ha : HasDerivAt a (H*a t) t)
    (hr : HasDerivAt rho (-3*H*(rho t+pressure)) t) :
    HasDerivAt (fun u => rho u*(a u)^2) (-H*(a t)^2*(rho t+3*pressure)) t := by
  convert hr.mul (ha.pow 2) using 1 <;> dsimp <;> ring

theorem weighted_density_scale_derivative_nonpositive (a rho : ℝ → ℝ) (t H pressure : ℝ)
    (ha : HasDerivAt a (H*a t) t)
    (hr : HasDerivAt rho (-3*H*(rho t+pressure)) t)
    (hH : 0 ≤ H) (hp : -rho t/3 ≤ pressure) :
    deriv (fun u => rho u*(a u)^2) t ≤ 0 := by
  rw [(weighted_density_scale_hasDerivAt a rho t H pressure ha hr).deriv]
  have h := mul_nonneg (mul_nonneg hH (sq_nonneg (a t))) (by linarith : 0 ≤ rho t+3*pressure)
  nlinarith

def constraint (c Lambda H rho : ℝ) : ℝ := H^2-c*rho-Lambda/3

theorem constraint_hasDerivAt (H rho : ℝ → ℝ) (t c Lambda pressure drho : ℝ)
    (hH : HasDerivAt H (-(3*c/2)*(rho t+pressure)) t)
    (hr : HasDerivAt rho drho t) :
    HasDerivAt (fun u => constraint c Lambda (H u) (rho u))
      (-c*(drho+3*H t*(rho t+pressure))) t := by
  have h := ((hH.pow 2).sub (hr.const_mul c)).sub_const (Lambda/3)
  convert h using 1 <;> dsimp [constraint] <;> ring

theorem constraint_zero_derivative (H rho : ℝ → ℝ) (t c Lambda pressure : ℝ)
    (hH : HasDerivAt H (-(3*c/2)*(rho t+pressure)) t)
    (hr : HasDerivAt rho (-3*H t*(rho t+pressure)) t) :
    HasDerivAt (fun u => constraint c Lambda (H u) (rho u)) 0 t := by
  convert constraint_hasDerivAt H rho t c Lambda pressure _ hH hr using 1 <;> ring

theorem constraint_preserved_on_interval (H rho pressure : ℝ → ℝ) (c Lambda l r t t0 : ℝ)
    (ht : t ∈ Set.Ioo l r) (ht0 : t0 ∈ Set.Ioo l r)
    (hH : ∀ u ∈ Set.Ioo l r, HasDerivAt H (-(3*c/2)*(rho u+pressure u)) u)
    (hr : ∀ u ∈ Set.Ioo l r, HasDerivAt rho (-3*H u*(rho u+pressure u)) u) :
    constraint c Lambda (H t) (rho t)=constraint c Lambda (H t0) (rho t0) := by
  have hz : ∀ u ∈ Set.Ioo l r,
      HasDerivAt (fun z => constraint c Lambda (H z) (rho z)) 0 u := by
    intro u hu
    exact constraint_zero_derivative H rho u c Lambda (pressure u) (hH u hu) (hr u hu)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo l r).isPreconnected
    (fun u hu => (hz u hu).differentiableAt.differentiableWithinAt)
    (fun u hu => (hz u hu).deriv) ht ht0

theorem coupled_constraint_hasDerivAt_zero {ι : Type*} [Fintype ι]
    (w k : ι → ℝ) (a H : ℝ → ℝ) (s v : ι → ℝ → ℝ) (t c Lambda : ℝ)
    (ha : HasDerivAt a (H t*a t) t) (ha0 : a t ≠ 0)
    (hs : ∀ i, HasDerivAt (s i) (v i t) t) (hs0 : ∀ i, s i t ≠ 0)
    (hv : ∀ i, HasDerivAt (v i)
      (-3*H t*v i t+(-(k i)^2/(a t)^2)*s i t+1/(4*(a t)^6*(s i t)^3)) t)
    (hH : HasDerivAt H (-(3*c/2)*
      (totalDensity w k (fun i => s i t) (fun i => v i t) (a t)+
       totalPressure w k (fun i => s i t) (fun i => v i t) (a t))) t) :
    HasDerivAt (fun u => constraint c Lambda (H u)
      (totalDensity w k (fun i => s i u) (fun i => v i u) (a u))) 0 t := by
  exact constraint_zero_derivative H
    (fun u => totalDensity w k (fun i => s i u) (fun i => v i u) (a u))
    t c Lambda (totalPressure w k (fun i => s i t) (fun i => v i t) (a t)) hH
    (total_density_hasDerivAt w k a s v t (H t) ha ha0 hs hs0 hv)

theorem constrained_hubble_nonzero (c Lambda H rho L a : ℝ)
    (hc : 0 < c) (hL : 0 < L) (ha : 0 < a) (hLambda : 0 ≤ Lambda)
    (hr : L/a^4 ≤ rho) (hC : H^2=c*rho+Lambda/3) : H ≠ 0 := by
  have hl : 0 < L/a^4 := by positivity
  have hp : 0 < c*rho := mul_pos hc (lt_of_lt_of_le hl hr)
  intro hzero
  rw [hzero] at hC
  nlinarith

theorem expanding_branch_preserved (H : ℝ → ℝ) (t0 t : ℝ)
    (htt : t0 ≤ t) (hcont : ContinuousOn H (Set.Icc t0 t))
    (hnever : ∀ u ∈ Set.Icc t0 t, H u ≠ 0) (hinit : 0 < H t0) : 0 < H t := by
  by_contra h
  have hle : H t ≤ 0 := le_of_not_gt h
  have hz : 0 ∈ H '' Set.Icc t0 t := intermediate_value_Icc' htt hcont ⟨hle, le_of_lt hinit⟩
  obtain ⟨u, hu, hz⟩ := hz
  exact hnever u hu hz

theorem hubble_scaled_lower_bound (c Lambda H rho L a : ℝ)
    (hc : 0 < c) (hL : 0 < L) (ha : 0 < a) (hLambda : 0 ≤ Lambda)
    (hH : 0 ≤ H) (hr : L/a^4 ≤ rho) (hC : H^2=c*rho+Lambda/3) :
    Real.sqrt (c*L) ≤ H*a^2 := by
  have ha4 : 0 < a^4 := by positivity
  have ha4n : a^4 ≠ 0 := ne_of_gt ha4
  have hprod := mul_le_mul_of_nonneg_right hr (le_of_lt ha4)
  have hlr : L ≤ rho*a^4 := by simpa only [div_mul_cancel₀ _ ha4n] using hprod
  have hcr := mul_le_mul_of_nonneg_left hlr (le_of_lt hc)
  have hlam : 0 ≤ (Lambda/3)*a^4 := by positivity
  have hscale : (H*a^2)^2=c*rho*a^4+(Lambda/3)*a^4 := by
    calc
      _ = H^2*a^4 := by ring
      _ = _ := by rw [hC]; ring
  have hs := Real.sq_sqrt (le_of_lt (mul_pos hc hL))
  have hs0 := Real.sqrt_nonneg (c*L)
  have hh : 0 ≤ H*a^2 := mul_nonneg hH (sq_nonneg a)
  nlinarith

theorem hubble_lower_bound (c Lambda H rho L a : ℝ)
    (hc : 0 < c) (hL : 0 < L) (ha : 0 < a) (hLambda : 0 ≤ Lambda)
    (hH : 0 ≤ H) (hr : L/a^4 ≤ rho) (hC : H^2=c*rho+Lambda/3) :
    Real.sqrt (c*L)/a^2 ≤ H := by
  have hs := hubble_scaled_lower_bound c Lambda H rho L a hc hL ha hLambda hH hr hC
  exact (div_le_iff₀ (by positivity : 0 < a^2)).mpr hs

theorem expanding_scale_square_growth (a : ℝ → ℝ) (t H c Lambda rho L : ℝ)
    (ha : HasDerivAt a (H*a t) t) (ha0 : 0 < a t)
    (hc : 0 < c) (hL : 0 < L) (hLambda : 0 ≤ Lambda) (hH : 0 ≤ H)
    (hr : L/(a t)^4 ≤ rho) (hC : H^2=c*rho+Lambda/3) :
    2*Real.sqrt (c*L) ≤ deriv (fun u => (a u)^2) t := by
  have hd : deriv (fun u => (a u)^2) t=2*H*(a t)^2 := by
    convert (ha.pow 2).deriv using 1 <;> dsimp <;> ring
  rw [hd]
  have hh := hubble_scaled_lower_bound c Lambda H rho L (a t) hc hL ha0 hLambda hH hr hC
  nlinarith

theorem hubble_derivative_nonpositive (H : ℝ → ℝ) (t c rho pressure : ℝ)
    (hc : 0 ≤ c) (hr : 0 ≤ rho) (hp : -rho/3 ≤ pressure)
    (hH : HasDerivAt H (-(3*c/2)*(rho+pressure)) t) : deriv H t ≤ 0 := by
  rw [hH.deriv]
  have hsum : 0 ≤ rho+pressure := by linarith
  have h := mul_nonneg hc hsum
  nlinarith

def covarianceCasimir (X Y Z : ℝ) : ℝ := X*Z-Y^2
def canonicalEnergy (k a X Z : ℝ) : ℝ := (Z/a^3+a*k^2*X)/2

theorem covariance_casimir_hasDerivAt (X Y Z : ℝ → ℝ) (t a k : ℝ)
    (hX : HasDerivAt X (2*Y t/a^3) t)
    (hY : HasDerivAt Y (Z t/a^3-a*k^2*X t) t)
    (hZ : HasDerivAt Z (-2*a*k^2*Y t) t) :
    HasDerivAt (fun u => covarianceCasimir (X u) (Y u) (Z u)) 0 t := by
  convert (hX.mul hZ).sub (hY.pow 2) using 1 <;> dsimp [covarianceCasimir] <;> ring

theorem covariance_casimir_preserved_on_interval (X Y Z a : ℝ → ℝ) (k l r t t0 : ℝ)
    (ht : t ∈ Set.Ioo l r) (ht0 : t0 ∈ Set.Ioo l r)
    (hX : ∀ u ∈ Set.Ioo l r, HasDerivAt X (2*Y u/(a u)^3) u)
    (hY : ∀ u ∈ Set.Ioo l r, HasDerivAt Y (Z u/(a u)^3-a u*k^2*X u) u)
    (hZ : ∀ u ∈ Set.Ioo l r, HasDerivAt Z (-2*a u*k^2*Y u) u) :
    covarianceCasimir (X t) (Y t) (Z t)=covarianceCasimir (X t0) (Y t0) (Z t0) := by
  have hz : ∀ u ∈ Set.Ioo l r,
      HasDerivAt (fun z => covarianceCasimir (X z) (Y z) (Z z)) 0 u := by
    intro u hu
    exact covariance_casimir_hasDerivAt X Y Z u (a u) k (hX u hu) (hY u hu) (hZ u hu)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo l r).isPreconnected
    (fun u hu => (hz u hu).differentiableAt.differentiableWithinAt)
    (fun u hu => (hz u hu).deriv) ht ht0

theorem pure_gaussian_covariance_casimir (a s v : ℝ) (hs : s ≠ 0) :
    covarianceCasimir (s^2) (a^3*s*v) (a^6*v^2+1/(4*s^2))=1/4 := by
  unfold covarianceCasimir
  field_simp
  <;> ring

theorem canonical_energy_density_identity (k a s v : ℝ) (ha : a ≠ 0) (hs : s ≠ 0) :
    canonicalEnergy k a (s^2) (a^6*v^2+1/(4*s^2))=a^3*modeDensity k a s v := by
  unfold canonicalEnergy modeDensity
  field_simp
  <;> ring

theorem canonical_gaussian_energy_lower_bound (k a s v : ℝ)
    (ha : 0 < a) (hs : s ≠ 0) :
    k/(2*a) ≤ canonicalEnergy k a (s^2) (a^6*v^2+1/(4*s^2)) := by
  have ha0 : a ≠ 0 := ne_of_gt ha
  calc
    _ = a^3*(k/(2*a^4)) := by field_simp <;> ring
    _ ≤ a^3*modeDensity k a s v :=
      mul_le_mul_of_nonneg_left (mode_uncertainty_lower_bound k a s v ha0 hs) (by positivity)
    _ = _ := (canonical_energy_density_identity k a s v ha0 hs).symm

theorem gaussian_energy_variance_expression_nonnegative (k a s v : ℝ)
    (hk : 0 ≤ k) (ha : 0 < a) (hs : s ≠ 0) :
    0 ≤ 2*(canonicalEnergy k a (s^2) (a^6*v^2+1/(4*s^2)))^2-k^2/(2*a^2) := by
  have hb := canonical_gaussian_energy_lower_bound k a s v ha hs
  have hg : 0 ≤ k/(2*a) := by positivity
  have he : 0 ≤ canonicalEnergy k a (s^2) (a^6*v^2+1/(4*s^2)) := le_trans hg hb
  have hp := mul_nonneg (sub_nonneg.mpr hb)
    (by linarith : 0 ≤ canonicalEnergy k a (s^2) (a^6*v^2+1/(4*s^2))+k/(2*a))
  have ha0 : a ≠ 0 := ne_of_gt ha
  have hi : 2*(k/(2*a))^2=k^2/(2*a^2) := by field_simp <;> ring
  nlinarith

end

end BackreactionCoupled
