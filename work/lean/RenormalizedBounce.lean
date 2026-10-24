import Mathlib

/-! Exact finite-dimensional counterterm-supported FLRW family.
Curvature expressions are explicit reduced formulas, not tensor definitions. -/

namespace RenormalizedBounce

noncomputable section

def familyBeta (u : ℝ) : ℝ := (2*u-1)/(6*u)
def familyLambda (u : ℝ) : ℝ := u-u^2
def familyRadiation (u : ℝ) : ℝ := u*(4*u-3)/2
def familyCounterterm (u z : ℝ) : ℝ := u^2*(6*z-9*z^2)
def familyCurvature (u z : ℝ) : ℝ := 6*u*(2-z)
def familyKretschmann (u z : ℝ) : ℝ := 12*u^2*(1+(1-z)^2)

theorem reduced_counterterm_identity (u z : ℝ) :
    12*u^2*z*(1-z)+(u*(1-z))*(6*u*(2-z))-(6*u*(2-z))^2/12=
      familyCounterterm u z := by
  unfold familyCounterterm
  ring

theorem family_constraint_identity (u z : ℝ) (hu : u ≠ 0) :
    u*(1-z)-(u*(1-z))^2-familyBeta u*familyCounterterm u z-
      familyLambda u-familyRadiation u*z^2=0 := by
  unfold familyBeta familyCounterterm familyLambda familyRadiation
  field_simp
  <;> ring

theorem family_box_curvature_identity (u z : ℝ) :
    -(12*u^2*(3*z^2-2*z))-3*(12*u^2*z*(1-z)) = -12*u^2*z := by ring

theorem family_trace_identity (u z : ℝ) (hu : u ≠ 0) :
    familyCurvature u z-12*familyLambda u-
      12*((u*(1-z))^2+(u*(1-z))*(u*z))+
      3*familyBeta u*(-12*u^2*z)=0 := by
  unfold familyCurvature familyLambda familyBeta
  field_simp
  <;> ring

theorem family_parameter_signs (u : ℝ) (hlo : (3:ℝ)/4 ≤ u) (hhi : u ≤ 1) :
    0 < u ∧ 0 < familyBeta u ∧ 0 ≤ familyLambda u ∧ 0 ≤ familyRadiation u := by
  have hu : 0 < u := by linarith
  constructor
  · exact hu
  constructor
  · unfold familyBeta
    exact div_pos (by linarith) (by positivity)
  constructor
  · unfold familyLambda
    have hp := mul_nonneg (le_of_lt hu) (by linarith : 0 ≤ 1-u)
    nlinarith
  · unfold familyRadiation
    have hp := mul_nonneg (le_of_lt hu) (by linarith : 0 ≤ 4*u-3)
    positivity

theorem family_curvature_bounds (u z : ℝ) (hu : 0 ≤ u)
    (hz : 0 ≤ z) (hz1 : z ≤ 1) :
    6*u ≤ familyCurvature u z ∧ familyCurvature u z ≤ 12*u := by
  unfold familyCurvature
  have h1 := mul_nonneg hu hz
  have h2 := mul_nonneg hu (by linarith : 0 ≤ 1-z)
  constructor <;> nlinarith

theorem family_kretschmann_bounds (u z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1) :
    12*u^2 ≤ familyKretschmann u z ∧ familyKretschmann u z ≤ 24*u^2 := by
  have hs0 : 0 ≤ (1-z)^2 := sq_nonneg _
  have hs1 : (1-z)^2 ≤ 1 := by nlinarith [mul_nonneg hz (sub_nonneg.mpr hz1)]
  have hp0 := mul_nonneg (sq_nonneg u) hs0
  have hp1 := mul_le_mul_of_nonneg_left hs1 (sq_nonneg u)
  unfold familyKretschmann
  constructor <;> nlinarith

theorem unit_parameters : familyBeta 1=1/6 ∧ familyLambda 1=0 ∧ familyRadiation 1=1/2 := by
  norm_num [familyBeta, familyLambda, familyRadiation]

theorem vacuum_parameters :
    familyBeta ((3:ℝ)/4)=1/9 ∧ familyLambda ((3:ℝ)/4)=3/16 ∧
      familyRadiation ((3:ℝ)/4)=0 := by
  norm_num [familyBeta, familyLambda, familyRadiation]

def bounceScale (ab t : ℝ) : ℝ := ab*Real.cosh t
def bounceHubble (t : ℝ) : ℝ := Real.tanh t
def bounceZ (t : ℝ) : ℝ := 1/(Real.cosh t)^2

theorem bounce_scale_minimum (ab t : ℝ) (hab : 0 < ab) :
    ab ≤ bounceScale ab t := by
  have hh := mul_le_mul_of_nonneg_left (Real.one_le_cosh t) (le_of_lt hab)
  simpa [bounceScale] using hh

theorem bounce_scale_positive (ab t : ℝ) (hab : 0 < ab) : 0 < bounceScale ab t := by
  exact lt_of_lt_of_le hab (bounce_scale_minimum ab t hab)

theorem bounce_scale_hasDerivAt (ab t : ℝ) :
    HasDerivAt (bounceScale ab) (bounceHubble t*bounceScale ab t) t := by
  have hc : Real.cosh t ≠ 0 := ne_of_gt (Real.cosh_pos t)
  convert (Real.hasDerivAt_cosh t).const_mul ab using 1 <;>
    dsimp [bounceScale, bounceHubble] <;> rw [Real.tanh_eq_sinh_div_cosh] <;>
    field_simp <;> ring

theorem bounce_hubble_hasDerivAt (t : ℝ) :
    HasDerivAt bounceHubble (bounceZ t) t := by
  have hc : Real.cosh t ≠ 0 := ne_of_gt (Real.cosh_pos t)
  have he : Real.cosh t*Real.cosh t-Real.sinh t*Real.sinh t=1 := by
    nlinarith [Real.cosh_sq_sub_sinh_sq t]
  have heq : bounceHubble=(fun u => Real.sinh u/Real.cosh u) := by
    funext u
    exact Real.tanh_eq_sinh_div_cosh u
  rw [heq]
  simpa only [bounceZ, he] using
    (Real.hasDerivAt_sinh t).div (Real.hasDerivAt_cosh t) hc

theorem bounce_z_hasDerivAt (t : ℝ) :
    HasDerivAt bounceZ (-2*bounceHubble t*bounceZ t) t := by
  have hc : Real.cosh t ≠ 0 := ne_of_gt (Real.cosh_pos t)
  have hc2 : (Real.cosh t)^2 ≠ 0 := pow_ne_zero 2 hc
  convert (hasDerivAt_const t 1).div ((Real.hasDerivAt_cosh t).pow 2) hc2 using 1 <;>
    dsimp [bounceZ, bounceHubble] <;> rw [Real.tanh_eq_sinh_div_cosh] <;>
    field_simp <;> ring

theorem bounce_z_bounds (t : ℝ) : 0 < bounceZ t ∧ bounceZ t ≤ 1 := by
  have hc : 0 < Real.cosh t := Real.cosh_pos t
  have hc1 := Real.one_le_cosh t
  unfold bounceZ
  constructor
  · positivity
  · apply (div_le_iff₀ (by positivity : 0 < (Real.cosh t)^2)).2
    nlinarith

theorem bounce_hubble_square (t : ℝ) : (bounceHubble t)^2=1-bounceZ t := by
  have hc : Real.cosh t ≠ 0 := ne_of_gt (Real.cosh_pos t)
  have hi := Real.cosh_sq_sub_sinh_sq t
  unfold bounceHubble bounceZ
  rw [Real.tanh_eq_sinh_div_cosh]
  field_simp
  nlinarith [hi]

theorem bounce_radiation_scaling (ab t : ℝ) (hab : ab ≠ 0) :
    (ab^4/2)/(bounceScale ab t)^4=(bounceZ t)^2/2 := by
  have hc : Real.cosh t ≠ 0 := ne_of_gt (Real.cosh_pos t)
  unfold bounceScale bounceZ
  field_simp
  <;> ring

theorem exact_unit_bounce_constraint (ab t : ℝ) (hab : ab ≠ 0) :
    (bounceHubble t)^2-(bounceHubble t)^4-
      (1/6)*familyCounterterm 1 (bounceZ t)-(ab^4/2)/(bounceScale ab t)^4=0 := by
  rw [bounce_radiation_scaling ab t hab]
  have hh := bounce_hubble_square t
  have hf := family_constraint_identity 1 (bounceZ t) (by norm_num)
  norm_num [familyBeta, familyLambda, familyRadiation] at hf
  nlinarith [hf]

theorem bounce_initial_data (ab : ℝ) :
    bounceScale ab 0=ab ∧ bounceHubble 0=0 ∧ bounceZ 0=1 := by
  simp [bounceScale, bounceHubble, bounceZ]

def bounceCurvature (t : ℝ) : ℝ := 6*(2-bounceZ t)
def bounceCurvatureVelocity (t : ℝ) : ℝ := 12*bounceHubble t*bounceZ t

theorem bounce_curvature_hasDerivAt (t : ℝ) :
    HasDerivAt bounceCurvature (bounceCurvatureVelocity t) t := by
  convert ((hasDerivAt_const t 2).sub (bounce_z_hasDerivAt t)).const_mul 6 using 1 <;>
    dsimp [bounceCurvature, bounceCurvatureVelocity] <;> ring

theorem bounce_curvature_velocity_hasDerivAt (t : ℝ) :
    HasDerivAt bounceCurvatureVelocity (12*(3*(bounceZ t)^2-2*bounceZ t)) t := by
  have hh := bounce_hubble_square t
  have hall := ((bounce_hubble_hasDerivAt t).mul (bounce_z_hasDerivAt t)).const_mul 12
  convert hall using 1
  · funext u
    dsimp [bounceCurvatureVelocity]
    ring
  · dsimp
    nlinarith [congrArg (fun x : ℝ => x*bounceZ t) hh]

theorem bounce_hubble_curvature_ode (t : ℝ) :
    HasDerivAt bounceHubble (bounceCurvature t/6-2*(bounceHubble t)^2) t := by
  have hh := bounce_hubble_square t
  convert bounce_hubble_hasDerivAt t using 1
  dsimp [bounceCurvature]
  linarith

theorem bounce_curvature_trace_ode (t : ℝ) :
    HasDerivAt bounceCurvatureVelocity
      ((bounceCurvature t+12*(bounceHubble t)^4-
        2*(bounceHubble t)^2*bounceCurvature t)/(3*(1/6))-
        3*bounceHubble t*bounceCurvatureVelocity t) t := by
  have hh := bounce_hubble_square t
  have hh4 : (bounceHubble t)^4=(1-bounceZ t)^2 := by
    calc
      _ = ((bounceHubble t)^2)^2 := by ring
      _ = _ := by rw [hh]
  convert bounce_curvature_velocity_hasDerivAt t using 1
  dsimp [bounceCurvature, bounceCurvatureVelocity]
  rw [hh4]
  ring_nf
  rw [hh]
  ring

def scaledScale (ab u t : ℝ) : ℝ := bounceScale ab (Real.sqrt u*t)
def scaledHubble (u t : ℝ) : ℝ := Real.sqrt u*bounceHubble (Real.sqrt u*t)
def scaledZ (u t : ℝ) : ℝ := bounceZ (Real.sqrt u*t)
def scaledCurvature (u t : ℝ) : ℝ := familyCurvature u (scaledZ u t)
def scaledCurvatureVelocity (u t : ℝ) : ℝ := 12*u*scaledHubble u t*scaledZ u t

theorem scaled_scale_minimum (ab u t : ℝ) (hab : 0 < ab) : ab ≤ scaledScale ab u t := by
  exact bounce_scale_minimum ab (Real.sqrt u*t) hab

theorem scaled_scale_hasDerivAt (ab u t : ℝ) :
    HasDerivAt (scaledScale ab u) (scaledHubble u t*scaledScale ab u t) t := by
  have ht := (hasDerivAt_id t).const_mul (Real.sqrt u)
  convert (bounce_scale_hasDerivAt ab (Real.sqrt u*t)).comp t ht using 1 <;>
    dsimp [scaledScale, scaledHubble] <;> ring

theorem scaled_hubble_hasDerivAt (u t : ℝ) (hu : 0 ≤ u) :
    HasDerivAt (scaledHubble u) (u*scaledZ u t) t := by
  have ht := (hasDerivAt_id t).const_mul (Real.sqrt u)
  have hh := ((bounce_hubble_hasDerivAt (Real.sqrt u*t)).comp t ht).const_mul (Real.sqrt u)
  convert hh using 1 <;> dsimp [scaledHubble, scaledZ]
  nlinarith [congrArg (fun x : ℝ => x*bounceZ (Real.sqrt u*t)) (Real.sq_sqrt hu)]

theorem scaled_z_hasDerivAt (u t : ℝ) :
    HasDerivAt (scaledZ u) (-2*scaledHubble u t*scaledZ u t) t := by
  have ht := (hasDerivAt_id t).const_mul (Real.sqrt u)
  convert (bounce_z_hasDerivAt (Real.sqrt u*t)).comp t ht using 1 <;>
    dsimp [scaledZ, scaledHubble] <;> ring

theorem scaled_hubble_square (u t : ℝ) (hu : 0 ≤ u) :
    (scaledHubble u t)^2=u*(1-scaledZ u t) := by
  unfold scaledHubble scaledZ
  rw [mul_pow, Real.sq_sqrt hu, bounce_hubble_square]

theorem scaled_curvature_hasDerivAt (u t : ℝ) :
    HasDerivAt (scaledCurvature u) (scaledCurvatureVelocity u t) t := by
  convert ((hasDerivAt_const t 2).sub (scaled_z_hasDerivAt u t)).const_mul (6*u) using 1 <;>
    dsimp [scaledCurvature, scaledCurvatureVelocity, familyCurvature] <;> ring

theorem scaled_curvature_velocity_hasDerivAt (u t : ℝ) (hu : 0 ≤ u) :
    HasDerivAt (scaledCurvatureVelocity u) (12*u^2*(3*(scaledZ u t)^2-2*scaledZ u t)) t := by
  have hh := scaled_hubble_square u t hu
  have hall := ((scaled_hubble_hasDerivAt u t hu).mul (scaled_z_hasDerivAt u t)).const_mul (12*u)
  convert hall using 1
  · funext x
    dsimp [scaledCurvatureVelocity]
    ring
  · dsimp
    nlinarith [congrArg (fun x : ℝ => u*x*scaledZ u t) hh]

theorem scaled_counterterm_identity (u t : ℝ) (hu : 0 ≤ u) :
    scaledHubble u t*scaledCurvatureVelocity u t+
      (scaledHubble u t)^2*scaledCurvature u t-(scaledCurvature u t)^2/12=
      familyCounterterm u (scaledZ u t) := by
  have hh := scaled_hubble_square u t hu
  dsimp [scaledCurvatureVelocity, scaledCurvature]
  have hm : scaledHubble u t*(12*u*scaledHubble u t*scaledZ u t)=
      12*u^2*scaledZ u t*(1-scaledZ u t) := by
    nlinarith [congrArg (fun x : ℝ => u*x*scaledZ u t) hh]
  rw [hm, hh]
  exact reduced_counterterm_identity u (scaledZ u t)

theorem scaled_radiation_identity (ab u t : ℝ) (hab : ab ≠ 0) :
    (familyRadiation u*ab^4)/(scaledScale ab u t)^4=
      familyRadiation u*(scaledZ u t)^2 := by
  have hc : Real.cosh (Real.sqrt u*t) ≠ 0 := ne_of_gt (Real.cosh_pos _)
  unfold scaledScale scaledZ bounceScale bounceZ
  field_simp
  <;> ring

theorem exact_scaled_bounce_constraint (ab u t : ℝ) (hab : ab ≠ 0) (hu : 0 < u) :
    (scaledHubble u t)^2-(scaledHubble u t)^4-
      familyBeta u*(scaledHubble u t*scaledCurvatureVelocity u t+
        (scaledHubble u t)^2*scaledCurvature u t-(scaledCurvature u t)^2/12)-
      familyLambda u-(familyRadiation u*ab^4)/(scaledScale ab u t)^4=0 := by
  rw [scaled_counterterm_identity u t (le_of_lt hu), scaled_radiation_identity ab u t hab]
  have hh := scaled_hubble_square u t (le_of_lt hu)
  have hh4 : (scaledHubble u t)^4=(u*(1-scaledZ u t))^2 := by
    calc
      _ = ((scaledHubble u t)^2)^2 := by ring
      _ = _ := by rw [hh]
  rw [hh4, hh]
  exact family_constraint_identity u (scaledZ u t) (ne_of_gt hu)

theorem exact_scaled_bounce_trace (u t : ℝ) (hu : 0 < u) :
    scaledCurvature u t-12*familyLambda u-
      12*((scaledHubble u t)^4+(scaledHubble u t)^2*deriv (scaledHubble u) t)+
      3*familyBeta u*(-deriv (scaledCurvatureVelocity u) t-
        3*scaledHubble u t*scaledCurvatureVelocity u t)=0 := by
  rw [(scaled_hubble_hasDerivAt u t (le_of_lt hu)).deriv,
    (scaled_curvature_velocity_hasDerivAt u t (le_of_lt hu)).deriv]
  have hh := scaled_hubble_square u t (le_of_lt hu)
  have hh4 : (scaledHubble u t)^4=(u*(1-scaledZ u t))^2 := by
    calc
      _ = ((scaledHubble u t)^2)^2 := by ring
      _ = _ := by rw [hh]
  have hm : scaledHubble u t*scaledCurvatureVelocity u t=
      12*u^2*scaledZ u t*(1-scaledZ u t) := by
    dsimp [scaledCurvatureVelocity]
    nlinarith [congrArg (fun x : ℝ => u*x*scaledZ u t) hh]
  have hm3 : 3*scaledHubble u t*scaledCurvatureVelocity u t=
      3*(12*u^2*scaledZ u t*(1-scaledZ u t)) := by nlinarith [hm]
  rw [hh4, hh, hm3, family_box_curvature_identity]
  exact family_trace_identity u (scaledZ u t) (ne_of_gt hu)

end

end RenormalizedBounce
