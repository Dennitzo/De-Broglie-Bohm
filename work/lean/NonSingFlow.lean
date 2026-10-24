import Mathlib

/-! The rescaled PG Dirac current for a positive areal radius. PDE existence is not asserted. -/
namespace NonSingFlow
noncomputable section

def upper (a b : ℝ) : ℝ := a^2+b^2
def lower (c d : ℝ) : ℝ := c^2+d^2
def density (a b c d : ℝ) : ℝ := upper a b+lower c d
def flux (v a b c d : ℝ) : ℝ := (1-v)*upper a b+(-1-v)*lower c d
def realOverlap (a b c d : ℝ) : ℝ := a*c+b*d
def imagOverlap (a b c d : ℝ) : ℝ := a*d-b*c
def fullNorm (v R theta a b c d : ℝ) : ℝ :=
  -1+(flux v a b c d/density a b c d+v)^2+
    R^2*(Real.sin theta)^2*(-2*imagOverlap a b c d/(R*density a b c d))^2

theorem pair_density_hasDerivAt (a b : ℝ → ℝ) (x da db : ℝ)
    (ha : HasDerivAt a da x) (hb : HasDerivAt b db x) :
    HasDerivAt (fun t => (a t)^2+(b t)^2) (2*(a x*da+b x*db)) x := by
  convert (ha.pow 2).add (hb.pow 2) using 1 <;> dsimp <;> ring

theorem hermitian_coupling_preserves_norm (a b c d : ℝ → ℝ) (t mu w : ℝ)
    (ha : HasDerivAt a (mu*d t+w*c t) t)
    (hb : HasDerivAt b (-mu*c t+w*d t) t)
    (hc : HasDerivAt c (mu*b t-w*a t) t)
    (hd : HasDerivAt d (-mu*a t-w*b t) t) :
    HasDerivAt (fun s => density (a s) (b s) (c s) (d s)) 0 t := by
  convert (pair_density_hasDerivAt a b t _ _ ha hb).add
    (pair_density_hasDerivAt c d t _ _ hc hd) using 1 <;>
      dsimp [density,upper,lower] <;> ring

theorem rescaled_dirac_current_continuity (a b c d : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu w ar br cr dr : ℝ) (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hcr : HasDerivAt (c t) cr r) (hdr : HasDerivAt (d t) dr r)
    (hat : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+mu*d t r+w*c t r) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r-mu*c t r+w*d t r) t)
    (hct : HasDerivAt (fun s => c s r) ((v r+1)*cr+vp/2*c t r+mu*b t r-w*a t r) t)
    (hdt : HasDerivAt (fun s => d s r) ((v r+1)*dr+vp/2*d t r-mu*a t r-w*b t r) t) :
    deriv (fun s => density (a s r) (b s r) (c s r) (d s r)) t+
      deriv (fun x => flux (v x) (a t x) (b t x) (c t x) (d t x)) r=0 := by
  have ht := (pair_density_hasDerivAt (fun s => a s r) (fun s => b s r) t _ _ hat hbt).add
    (pair_density_hasDerivAt (fun s => c s r) (fun s => d s r) t _ _ hct hdt)
  have hu := ((hasDerivAt_const r (1 : ℝ)).sub hv).mul
    (pair_density_hasDerivAt (a t) (b t) r ar br har hbr)
  have hd := ((hasDerivAt_const r (-1 : ℝ)).sub hv).mul
    (pair_density_hasDerivAt (c t) (d t) r cr dr hcr hdr)
  have hr := hu.add hd
  change HasDerivAt (fun s => density (a s r) (b s r) (c s r) (d s r)) _ t at ht
  change HasDerivAt (fun x => flux (v x) (a t x) (b t x) (c t x) (d t x)) _ r at hr
  rw [ht.deriv,hr.deriv]
  dsimp
  ring

theorem overlap_identity (a b c d : ℝ) :
    (realOverlap a b c d)^2+(imagOverlap a b c d)^2=upper a b*lower c d := by
  unfold realOverlap imagOverlap upper lower
  ring

theorem radial_norm_identity (v a b c d : ℝ) (hq : density a b c d ≠ 0) :
    -1+(flux v a b c d/density a b c d+v)^2 =
      -4*upper a b*lower c d/(density a b c d)^2 := by
  field_simp
  unfold density flux upper lower
  ring

theorem angular_metric_norm_identity (v R theta a b c d : ℝ)
    (hR : R ≠ 0) (hq : density a b c d ≠ 0) :
    fullNorm v R theta a b c d =
      -4*((realOverlap a b c d)^2+(imagOverlap a b c d)^2*(Real.cos theta)^2)/(density a b c d)^2 := by
  unfold fullNorm
  rw [radial_norm_identity v a b c d hq]
  have hh := congrArg (fun z : ℝ => (imagOverlap a b c d)^2*z) (Real.sin_sq_add_cos_sq theta)
  have hn := overlap_identity a b c d
  dsimp at hh
  field_simp
  nlinarith [hn,hh]

theorem angular_current_causal (v R theta a b c d : ℝ)
    (hR : R ≠ 0) (hq : density a b c d ≠ 0) : fullNorm v R theta a b c d ≤ 0 := by
  rw [angular_metric_norm_identity v R theta a b c d hR hq]
  apply div_nonpos_of_nonpos_of_nonneg
  · have hh : 0 ≤ (realOverlap a b c d)^2+(imagOverlap a b c d)^2*(Real.cos theta)^2 := by positivity
    linarith
  · positivity

end
end NonSingFlow
