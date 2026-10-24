import Mathlib

/-! Reduced causal-current and characteristic-transport identities in
Painleve-Gullstrand coordinates. A Dirac spin geometry is not constructed. -/

namespace BlackHoleCurrent

noncomputable section

theorem causal_velocity_bounds (jT jr v : ℝ) (hT : 0 < jT)
    (hc : (jr+v*jT)^2 ≤ jT^2) :
    -1-v ≤ jr/jT ∧ jr/jT ≤ 1-v := by
  have habs := (sq_le_sq).mp hc
  rw [abs_of_pos hT] at habs
  have hh := abs_le.mp habs
  constructor
  · apply (le_div_iff₀ hT).2
    nlinarith [hh.1]
  · apply (div_le_iff₀ hT).2
    nlinarith [hh.2]

theorem interior_future_current_ingoing (jT jr v : ℝ) (hT : 0 < jT)
    (hc : (jr+v*jT)^2 ≤ jT^2) (hv : 1 < v) : jr/jT < 0 := by
  have hb := (causal_velocity_bounds jT jr v hT hc).2
  linarith

theorem horizon_causal_bounds (jT jr : ℝ) (hT : 0 < jT)
    (hc : (jr+jT)^2 ≤ jT^2) : -2 ≤ jr/jT ∧ jr/jT ≤ 0 := by
  have hh := causal_velocity_bounds jT jr 1 hT (by simpa using hc)
  norm_num at hh
  exact hh

def density (P v g : ℝ → ℝ) (T r : ℝ) : ℝ := P (T+g r)/(1+v r)
def flux (P g : ℝ → ℝ) (T r : ℝ) : ℝ := -P (T+g r)
def ingoingVelocity (v : ℝ) : ℝ := -1-v

theorem guidance_identity (P v g : ℝ → ℝ) (T r : ℝ) (hv : 0 ≤ v r) :
    ingoingVelocity (v r)*density P v g T r=flux P g T r := by
  have hd : 1+v r ≠ 0 := by positivity
  unfold ingoingVelocity density flux
  field_simp
  <;> ring

theorem density_positive (P v g : ℝ → ℝ) (T r : ℝ)
    (hP : 0 < P (T+g r)) (hv : 0 ≤ v r) : 0 < density P v g T r := by
  unfold density
  exact div_pos hP (by positivity)

theorem ingoing_current_null (P v g : ℝ → ℝ) (T r : ℝ) (hv : 0 ≤ v r) :
    (flux P g T r+v r*density P v g T r)^2=(density P v g T r)^2 := by
  have hh := guidance_identity P v g T r hv
  unfold ingoingVelocity at hh
  have he : flux P g T r+v r*density P v g T r= -density P v g T r := by nlinarith [hh]
  rw [he]
  ring

theorem density_time_hasDerivAt (P v g : ℝ → ℝ) (T r p1 : ℝ)
    (hP : HasDerivAt P p1 (T+g r)) :
    HasDerivAt (fun t => density P v g t r) (p1/(1+v r)) T := by
  convert (hP.comp T ((hasDerivAt_id T).add_const (g r))).div_const (1+v r) using 1 <;>
    dsimp [density] <;> ring

theorem flux_radial_hasDerivAt (P v g : ℝ → ℝ) (T r p1 : ℝ)
    (hP : HasDerivAt P p1 (T+g r))
    (hg : HasDerivAt g (1/(1+v r)) r) :
    HasDerivAt (fun x => flux P g T x) (-p1/(1+v r)) r := by
  convert (hP.comp r ((hasDerivAt_const r T).add hg)).neg using 1 <;>
    dsimp [flux] <;> ring

theorem current_continuity (P v g : ℝ → ℝ) (T r p1 : ℝ)
    (hP : HasDerivAt P p1 (T+g r))
    (hg : HasDerivAt g (1/(1+v r)) r) :
    deriv (fun t => density P v g t r) T+deriv (fun x => flux P g T x) r=0 := by
  rw [(density_time_hasDerivAt P v g T r p1 hP).deriv,
    (flux_radial_hasDerivAt P v g T r p1 hP hg).deriv]
  ring

theorem characteristic_coordinate_hasDerivAt (r g : ℝ → ℝ) (T v : ℝ)
    (hv : 1+v ≠ 0)
    (hr : HasDerivAt r (-1-v) T)
    (hg : HasDerivAt g (1/(1+v)) (r T)) :
    HasDerivAt (fun t => t+g (r t)) 0 T := by
  convert (hasDerivAt_id T).add (hg.comp T hr) using 1 <;> dsimp <;> field_simp <;> ring

def pgPrimitive (a r : ℝ) : ℝ :=
  r-2*a*Real.sqrt r+2*a^2*Real.log (Real.sqrt r+a)
def pgShift (a r : ℝ) : ℝ := a/Real.sqrt r

theorem pg_primitive_hasDerivAt (a r : ℝ) (ha : 0 < a) (hr : 0 < r) :
    HasDerivAt (pgPrimitive a) (1/(1+pgShift a r)) r := by
  have hs : 0 < Real.sqrt r := Real.sqrt_pos.2 hr
  have hs0 : Real.sqrt r ≠ 0 := ne_of_gt hs
  have hsum : Real.sqrt r+a ≠ 0 := by positivity
  have hd : 1+a/Real.sqrt r ≠ 0 := by positivity
  have hroot := Real.hasDerivAt_sqrt (ne_of_gt hr)
  have hlog := (Real.hasDerivAt_log hsum).comp r (hroot.add_const a)
  have hall := ((hasDerivAt_id r).sub (hroot.const_mul (2*a))).add (hlog.const_mul (2*a^2))
  convert hall using 1 <;> dsimp [pgPrimitive, pgShift] <;> field_simp <;> ring

theorem pg_shift_positive (a r : ℝ) (ha : 0 < a) (hr : 0 < r) : 0 < pgShift a r := by
  unfold pgShift
  exact div_pos ha (Real.sqrt_pos.2 hr)

theorem pg_current_continuity (P : ℝ → ℝ) (a T r p1 : ℝ)
    (ha : 0 < a) (hr : 0 < r) (hP : HasDerivAt P p1 (T+pgPrimitive a r)) :
    deriv (fun t => density P (pgShift a) (pgPrimitive a) t r) T+
      deriv (fun x => flux P (pgPrimitive a) T x) r=0 := by
  exact current_continuity P (pgShift a) (pgPrimitive a) T r p1 hP
    (pg_primitive_hasDerivAt a r ha hr)

def ingoingAmplitude (F v g : ℝ → ℝ) (T r : ℝ) : ℝ :=
  F (T+g r)/Real.sqrt (1+v r)

theorem amplitude_density_identity (F v g : ℝ → ℝ) (T r : ℝ) (hv : 0 ≤ v r) :
    (ingoingAmplitude F v g T r)^2=density (fun x => (F x)^2) v g T r := by
  unfold ingoingAmplitude density
  rw [div_pow, Real.sq_sqrt (by positivity : 0 ≤ 1+v r)]

theorem amplitude_time_hasDerivAt (F v g : ℝ → ℝ) (T r F1 : ℝ)
    (hF : HasDerivAt F F1 (T+g r)) :
    HasDerivAt (fun t => ingoingAmplitude F v g t r)
      (F1/Real.sqrt (1+v r)) T := by
  convert (hF.comp T ((hasDerivAt_id T).add_const (g r))).div_const
    (Real.sqrt (1+v r)) using 1 <;> dsimp [ingoingAmplitude] <;> ring

theorem amplitude_radial_hasDerivAt (F v g : ℝ → ℝ) (T r F1 v1 : ℝ)
    (hv0 : 0 ≤ v r) (hF : HasDerivAt F F1 (T+g r))
    (hv : HasDerivAt v v1 r) (hg : HasDerivAt g (1/(1+v r)) r) :
    HasDerivAt (fun x => ingoingAmplitude F v g T x)
      (F1/((1+v r)*Real.sqrt (1+v r))-
        F (T+g r)*v1/(2*(1+v r)*Real.sqrt (1+v r))) r := by
  have hd : 0 < 1+v r := by positivity
  have hs : 0 < Real.sqrt (1+v r) := Real.sqrt_pos.2 hd
  have hs0 : Real.sqrt (1+v r) ≠ 0 := ne_of_gt hs
  have hs2 := Real.sq_sqrt (le_of_lt hd)
  have hn := hF.comp r ((hasDerivAt_const r T).add hg)
  have hden := (Real.hasDerivAt_sqrt (ne_of_gt hd)).comp r ((hasDerivAt_const r 1).add hv)
  convert hn.div hden hs0 using 1
  dsimp
  generalize he : Real.sqrt (1+v r)=s at *
  rw [← hs2]
  field_simp
  <;> ring

theorem ingoing_chiral_equation (F v g : ℝ → ℝ) (T r F1 v1 : ℝ)
    (hv0 : 0 ≤ v r) (hF : HasDerivAt F F1 (T+g r))
    (hv : HasDerivAt v v1 r) (hg : HasDerivAt g (1/(1+v r)) r) :
    deriv (fun t => ingoingAmplitude F v g t r) T+
      (-1-v r)*deriv (fun x => ingoingAmplitude F v g T x) r-
      (v1/2)*ingoingAmplitude F v g T r=0 := by
  rw [(amplitude_time_hasDerivAt F v g T r F1 hF).deriv,
    (amplitude_radial_hasDerivAt F v g T r F1 v1 hv0 hF hv hg).deriv]
  have hd : 1+v r ≠ 0 := by positivity
  have hs : Real.sqrt (1+v r) ≠ 0 := ne_of_gt (Real.sqrt_pos.2 (by positivity))
  unfold ingoingAmplitude
  field_simp
  <;> ring

end

end BlackHoleCurrent
