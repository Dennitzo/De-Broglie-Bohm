import Mathlib

/-! Real angular coupling and the specified kappa=1,m_j=1/2 angular reconstruction. -/
namespace MassiveAngular
noncomputable section

def upper (a b : ℝ) : ℝ := a^2+b^2
def lower (c e : ℝ) : ℝ := c^2+e^2
def density (a b c e : ℝ) : ℝ := upper a b+lower c e
def flux (v a b c e : ℝ) : ℝ := (1-v)*upper a b+(-1-v)*lower c e
def reOverlap (a b c e : ℝ) : ℝ := a*c+b*e
def imOverlap (a b c e : ℝ) : ℝ := a*e-b*c
def transfer (mu w a b c e : ℝ) : ℝ :=
  2*mu*imOverlap a b c e+2*w*reOverlap a b c e

theorem overlap_norm_identity (a b c e : ℝ) :
    (reOverlap a b c e)^2+(imOverlap a b c e)^2=upper a b*lower c e := by
  unfold reOverlap imOverlap upper lower
  ring

theorem square_pair_hasDerivAt (a b : ℝ → ℝ) (x da db : ℝ)
    (ha : HasDerivAt a da x) (hb : HasDerivAt b db x) :
    HasDerivAt (fun t => (a t)^2+(b t)^2) (2*(a x*da+b x*db)) x := by
  convert (ha.pow 2).add (hb.pow 2) using 1 <;> dsimp <;> ring

theorem chiral_balance_hasDerivAt (a b : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r chi vp ar br sa sb : ℝ)
    (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hat : HasDerivAt (fun s => a s r) ((v r-chi)*ar+vp/2*a t r+sa) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-chi)*br+vp/2*b t r+sb) t) :
    deriv (fun s => (a s r)^2+(b s r)^2) t+
      deriv (fun x => (chi-v x)*((a t x)^2+(b t x)^2)) r =
      2*(a t r*sa+b t r*sb) := by
  have ht := square_pair_hasDerivAt (fun s => a s r) (fun s => b s r) t _ _ hat hbt
  have hr := ((hasDerivAt_const r chi).sub hv).mul (square_pair_hasDerivAt (a t) (b t) r ar br har hbr)
  change HasDerivAt (fun x => (chi-v x)*((a t x)^2+(b t x)^2)) _ r at hr
  rw [ht.deriv,hr.deriv]
  dsimp
  ring

theorem angular_upper_mass_exchange (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu w ar br : ℝ) (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hat : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+mu*e t r+w*c t r) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r-mu*c t r+w*e t r) t) :
    deriv (fun s => upper (a s r) (b s r)) t+
      deriv (fun x => (1-v x)*upper (a t x) (b t x)) r =
      transfer mu w (a t r) (b t r) (c t r) (e t r) := by
  have ha : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+(mu*e t r+w*c t r)) t := by
    convert hat using 1 <;> ring
  have hb : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r+(-mu*c t r+w*e t r)) t := by
    convert hbt using 1 <;> ring
  have hh := chiral_balance_hasDerivAt a b v t r 1 vp ar br
    (mu*e t r+w*c t r) (-mu*c t r+w*e t r) hv har hbr ha hb
  convert hh using 1
  unfold transfer reOverlap imOverlap
  ring

theorem angular_lower_mass_exchange (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu w cr er : ℝ) (hv : HasDerivAt v vp r)
    (hcr : HasDerivAt (c t) cr r) (her : HasDerivAt (e t) er r)
    (hct : HasDerivAt (fun s => c s r) ((v r+1)*cr+vp/2*c t r+mu*b t r-w*a t r) t)
    (het : HasDerivAt (fun s => e s r) ((v r+1)*er+vp/2*e t r-mu*a t r-w*b t r) t) :
    deriv (fun s => lower (c s r) (e s r)) t+
      deriv (fun x => (-1-v x)*lower (c t x) (e t x)) r =
      -transfer mu w (a t r) (b t r) (c t r) (e t r) := by
  have hc : HasDerivAt (fun s => c s r) ((v r-(-1))*cr+vp/2*c t r+(mu*b t r-w*a t r)) t := by
    convert hct using 1 <;> ring
  have he : HasDerivAt (fun s => e s r) ((v r-(-1))*er+vp/2*e t r+(-mu*a t r-w*b t r)) t := by
    convert het using 1 <;> ring
  have hh := chiral_balance_hasDerivAt c e v t r (-1) vp cr er
    (mu*b t r-w*a t r) (-mu*a t r-w*b t r) hv hcr her hc he
  convert hh using 1
  unfold transfer reOverlap imOverlap
  ring

theorem angular_dirac_continuity (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu w ar br cr er : ℝ) (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hcr : HasDerivAt (c t) cr r) (her : HasDerivAt (e t) er r)
    (hat : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+mu*e t r+w*c t r) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r-mu*c t r+w*e t r) t)
    (hct : HasDerivAt (fun s => c s r) ((v r+1)*cr+vp/2*c t r+mu*b t r-w*a t r) t)
    (het : HasDerivAt (fun s => e s r) ((v r+1)*er+vp/2*e t r-mu*a t r-w*b t r) t) :
    deriv (fun s => density (a s r) (b s r) (c s r) (e s r)) t+
      deriv (fun x => flux (v x) (a t x) (b t x) (c t x) (e t x)) r=0 := by
  have hu := angular_upper_mass_exchange a b c e v t r vp mu w ar br hv har hbr hat hbt
  have hd := angular_lower_mass_exchange a b c e v t r vp mu w cr er hv hcr her hct het
  have htu := square_pair_hasDerivAt (fun s => a s r) (fun s => b s r) t _ _ hat hbt
  have htd := square_pair_hasDerivAt (fun s => c s r) (fun s => e s r) t _ _ hct het
  have hru := ((hasDerivAt_const r (1 : ℝ)).sub hv).mul (square_pair_hasDerivAt (a t) (b t) r ar br har hbr)
  have hrd := ((hasDerivAt_const r (-1 : ℝ)).sub hv).mul (square_pair_hasDerivAt (c t) (e t) r cr er hcr her)
  have htu' : DifferentiableAt ℝ (fun s => upper (a s r) (b s r)) t := htu.differentiableAt
  have htd' : DifferentiableAt ℝ (fun s => lower (c s r) (e s r)) t := htd.differentiableAt
  have hru' : DifferentiableAt ℝ (fun x => (1-v x)*upper (a t x) (b t x)) r := hru.differentiableAt
  have hrd' : DifferentiableAt ℝ (fun x => (-1-v x)*lower (c t x) (e t x)) r := hrd.differentiableAt
  have ht : deriv (fun s => upper (a s r) (b s r)+lower (c s r) (e s r)) t =
      deriv (fun s => upper (a s r) (b s r)) t+deriv (fun s => lower (c s r) (e s r)) t := by
    simpa using deriv_add htu' htd'
  have hr : deriv (fun x => (1-v x)*upper (a t x) (b t x)+(-1-v x)*lower (c t x) (e t x)) r =
      deriv (fun x => (1-v x)*upper (a t x) (b t x)) r+
      deriv (fun x => (-1-v x)*lower (c t x) (e t x)) r := by
    simpa using deriv_add hru' hrd'
  unfold density flux
  rw [ht,hr]
  linarith

def radialVelocity (v a b c e : ℝ) : ℝ := flux v a b c e/density a b c e
def azimuthVelocity (r a b c e : ℝ) : ℝ := -2*imOverlap a b c e/(r*density a b c e)
def fullNorm (v r theta a b c e : ℝ) : ℝ :=
  -1+(radialVelocity v a b c e+v)^2+
    r^2*(Real.sin theta)^2*(azimuthVelocity r a b c e)^2

theorem radial_norm_identity (v a b c e : ℝ) (hq : density a b c e ≠ 0) :
    -1+(radialVelocity v a b c e+v)^2 =
      -4*upper a b*lower c e/(density a b c e)^2 := by
  unfold radialVelocity
  field_simp
  unfold density flux upper lower
  ring

theorem full_angular_metric_norm (v r theta a b c e : ℝ)
    (hr : r ≠ 0) (hq : density a b c e ≠ 0) :
    fullNorm v r theta a b c e =
      -4*((reOverlap a b c e)^2+(imOverlap a b c e)^2*(Real.cos theta)^2)/(density a b c e)^2 := by
  unfold fullNorm
  rw [radial_norm_identity v a b c e hq]
  unfold azimuthVelocity
  have hh := congrArg (fun z : ℝ => (imOverlap a b c e)^2*z) (Real.sin_sq_add_cos_sq theta)
  have hn := overlap_norm_identity a b c e
  dsimp at hh
  field_simp
  nlinarith [hn,hh]

theorem full_angular_current_causal (v r theta a b c e : ℝ)
    (hr : r ≠ 0) (hq : density a b c e ≠ 0) :
    fullNorm v r theta a b c e ≤ 0 := by
  rw [full_angular_metric_norm v r theta a b c e hr hq]
  apply div_nonpos_of_nonpos_of_nonneg
  · have hs : 0 ≤ (reOverlap a b c e)^2+(imOverlap a b c e)^2*(Real.cos theta)^2 := by positivity
    linarith
  · positivity

theorem real_overlap_strictly_timelike (v r theta a b c e : ℝ)
    (hr : r ≠ 0) (hq : density a b c e ≠ 0) (hz : reOverlap a b c e ≠ 0) :
    fullNorm v r theta a b c e < 0 := by
  rw [full_angular_metric_norm v r theta a b c e hr hq]
  apply div_neg_of_neg_of_pos
  · have hs : 0 < (reOverlap a b c e)^2+(imOverlap a b c e)^2*(Real.cos theta)^2 := by
      have hh := sq_pos_of_ne_zero hz
      positivity
    linarith
  · exact sq_pos_of_ne_zero hq

theorem zero_initial_imaginary_overlap_is_radial (r a b c e : ℝ)
    (hi : imOverlap a b c e=0) : azimuthVelocity r a b c e=0 := by
  simp [azimuthVelocity,hi]

theorem equatorial_imaginary_overlap_null (v r : ℝ) (hr : r ≠ 0) :
    fullNorm v r (Real.pi/2) 1 0 0 1=0 := by
  rw [full_angular_metric_norm v r (Real.pi/2) 1 0 0 1 hr (by norm_num [density,upper,lower])]
  simp [reOverlap,imOverlap]

end
end MassiveAngular
