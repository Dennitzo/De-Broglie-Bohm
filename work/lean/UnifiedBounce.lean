import Mathlib

/-! Rational positive-frequency state: regular coefficients and local bounce. -/
namespace UnifiedBounce
noncomputable section
def D (s a p : ℝ) := (s^2+(p+a)^2)*(s^2+(p-a)^2)
def density (s a p : ℝ) := 4*(s^2+p^2)/D s a p
def pa (s a p : ℝ) := -4*s*a*p/D s a p
def pp (s a p : ℝ) := s*((s^2+p^2)^2+4*a^2*p^2-a^4)/((s^2+p^2)*D s a p)
def turnSpeed (s a : ℝ) := (s^2-a^2)/(s*(s^2+a^2))
def streamF (G : ℝ → ℝ) (s a p : ℝ) := G (p+a)-G (p-a)+2*s*a/D s a p
def primitive (s x : ℝ) := x/(2*s*(s^2+x^2))+Real.arctan (x/s)/(2*s^2)

theorem primitive_derivative (s x : ℝ) (hs : 0 < s) :
    HasDerivAt (primitive s) (s/(s^2+x^2)^2) x := by
  have hs0 : s ≠ 0 := ne_of_gt hs
  have hb : s^2+x^2 ≠ 0 := ne_of_gt (show 0 < s^2+x^2 by positivity)
  have hd : 2*s*(s^2+x^2) ≠ 0 := by positivity
  have hden := ((hasDerivAt_const x (s^2)).add ((hasDerivAt_id x).pow 2)).const_mul (2*s)
  have hfirst := (hasDerivAt_id x).div hden hd
  have hsecond := ((Real.hasDerivAt_arctan (x/s)).comp x ((hasDerivAt_id x).div_const s)).div_const (2*s^2)
  convert hfirst.add hsecond using 1
  have hother : 1+(x/s)^2 ≠ 0 := ne_of_gt (show 0 < 1+(x/s)^2 by positivity)
  simp only [Pi.add_apply,Pi.pow_apply,id_eq]
  field_simp
  <;> ring

theorem denominator_forms (s a p : ℝ) :
    D s a p=(s^2+a^2-p^2)^2+4*s^2*p^2 := by
  dsimp [D]; ring

theorem denominator_positive (s a p : ℝ) (hs : 0 < s) : 0 < D s a p := by
  dsimp [D]
  positivity

theorem density_positive (s a p : ℝ) (hs : 0 < s) : 0 < density s a p := by
  dsimp [density]
  exact div_pos (by positivity) (denominator_positive s a p hs)

theorem alpha_momentum_sign (s a p : ℝ) (hs : 0 < s) (ha : 0 < a) (hp : 0 < p) :
    pa s a p < 0 := by
  dsimp [pa]
  have hn : -4*s*a*p < 0 := by nlinarith [mul_pos (mul_pos hs ha) hp]
  exact div_neg_of_neg_of_pos hn (denominator_positive s a p hs)

theorem alpha_momentum_turn (s a : ℝ) : pa s a 0=0 := by
  simp [pa]

theorem phi_momentum_turn (s a : ℝ) (hs : 0 < s) : pp s a 0=turnSpeed s a := by
  have hs0 : s ≠ 0 := ne_of_gt hs
  have hd : s^2+a^2 ≠ 0 := ne_of_gt (show 0 < s^2+a^2 by positivity)
  dsimp [pp,turnSpeed,D]
  field_simp
  <;> ring

theorem turn_speed_positive (s a : ℝ) (ha : 0 < a) (has : a < s) :
    0 < turnSpeed s a := by
  have hs : 0 < s := lt_trans ha has
  dsimp [turnSpeed]
  apply div_pos
  · nlinarith
  · positivity

theorem density_even_phi (s a p : ℝ) : density s a (-p)=density s a p := by
  have hd : D s a (-p)=D s a p := by dsimp [D]; ring
  dsimp [density]
  rw [hd]
  ring

theorem alpha_momentum_odd_phi (s a p : ℝ) : pa s a (-p)= -pa s a p := by
  have hd : D s a (-p)=D s a p := by dsimp [D]; ring
  dsimp [pa]
  rw [hd]
  ring

theorem phi_momentum_even_phi (s a p : ℝ) : pp s a (-p)=pp s a p := by
  have hd : D s a (-p)=D s a p := by dsimp [D]; ring
  dsimp [pp]
  rw [hd]
  ring

theorem quantum_potential_turn_negative (s a : ℝ) (ha : 0 < a) (has : a < s) :
    ((pa s a 0)^2-(pp s a 0)^2)/2 < 0 := by
  have hs : 0 < s := lt_trans ha has
  rw [alpha_momentum_turn,phi_momentum_turn s a hs]
  have hp := turn_speed_positive s a ha has
  nlinarith [sq_pos_of_pos hp]

theorem alpha_momentum_turn_derivative (a p : ℝ → ℝ) (t s a0 v : ℝ)
    (hs : 0 < s) (ha0 : a t=a0) (hp0 : p t=0)
    (ha : HasDerivAt a 0 t) (hp : HasDerivAt p v t) :
    HasDerivAt (fun x => pa s (a x) (p x))
      (-4*s*a0*v/(s^2+a0^2)^2) t := by
  have hd : D s (a t) (p t) ≠ 0 := ne_of_gt (denominator_positive s (a t) (p t) hs)
  have hden := ((hasDerivAt_const t (s^2)).add ((hp.add ha).pow 2)).mul
    ((hasDerivAt_const t (s^2)).add ((hp.sub ha).pow 2))
  have hnum := (ha.const_mul (-4*s)).mul hp
  have hquot := hnum.div hden hd
  convert hquot using 1
  simp only [Pi.pow_apply,Pi.mul_apply,Pi.add_apply,Pi.sub_apply,ha0,hp0]
  have hne : s^2+a0^2 ≠ 0 := ne_of_gt (show 0 < s^2+a0^2 by positivity)
  dsimp [D]
  field_simp
  <;> ring

theorem proper_bounce_derivative (a p u : ℝ → ℝ) (t s a0 v b : ℝ)
    (hs : 0 < s) (ha0 : a t=a0) (hp0 : p t=0)
    (ha : HasDerivAt a 0 t) (hp : HasDerivAt p v t)
    (hu : HasDerivAt u 0 t) :
    HasDerivAt (fun x => -b*u x*pa s (a x) (p x))
      (4*b*u t*s*a0*v/(s^2+a0^2)^2) t := by
  have hd := alpha_momentum_turn_derivative a p t s a0 v hs ha0 hp0 ha hp
  convert (hu.const_mul (-b)).mul hd using 1
  ring

theorem proper_bounce_acceleration_formula (s a b u : ℝ) (hs : 0 < s) :
    4*b*u*s*a*(b*u*turnSpeed s a)/(s^2+a^2)^2=
      4*b^2*u^2*a*(s^2-a^2)/(s^2+a^2)^3 := by
  have hs0 : s ≠ 0 := ne_of_gt hs
  have hd : s^2+a^2 ≠ 0 := ne_of_gt (show 0 < s^2+a^2 by positivity)
  dsimp [turnSpeed]
  field_simp
  <;> ring

theorem proper_bounce_acceleration_positive (s a b u : ℝ)
    (ha : 0 < a) (has : a < s) (hb : 0 < b) (hu : 0 < u) :
    0 < 4*b^2*u^2*a*(s^2-a^2)/(s^2+a^2)^3 := by
  have hs : 0 < s := lt_trans ha has
  have hsa : 0 < s^2-a^2 := by nlinarith
  positivity

theorem stream_alpha_derivative (G : ℝ → ℝ) (s a p : ℝ) (hs : 0 < s)
    (hplus : HasDerivAt G (s/(s^2+(p+a)^2)^2) (p+a))
    (hminus : HasDerivAt G (s/(s^2+(p-a)^2)^2) (p-a)) :
    HasDerivAt (fun x => streamF G s x p) (density s a p*pp s a p) a := by
  have hpa := (hasDerivAt_const a p).add (hasDerivAt_id a)
  have hma := (hasDerivAt_const a p).sub (hasDerivAt_id a)
  have hden := ((hasDerivAt_const a (s^2)).add (hpa.pow 2)).mul
    ((hasDerivAt_const a (s^2)).add (hma.pow 2))
  have hd : D s a p ≠ 0 := ne_of_gt (denominator_positive s a p hs)
  have hquot := ((hasDerivAt_id a).const_mul (2*s)).div hden hd
  have hx := ((hplus.comp a hpa).sub (hminus.comp a hma)).add hquot
  convert hx using 1
  have h1 : s^2+(p+a)^2 ≠ 0 := ne_of_gt (show 0 < s^2+(p+a)^2 by positivity)
  have h2 : s^2+(p-a)^2 ≠ 0 := ne_of_gt (show 0 < s^2+(p-a)^2 by positivity)
  have h3 : s^2+p^2 ≠ 0 := ne_of_gt (show 0 < s^2+p^2 by positivity)
  dsimp [density,pp,D]
  field_simp
  <;> ring

theorem stream_phi_derivative (G : ℝ → ℝ) (s a p : ℝ) (hs : 0 < s)
    (hplus : HasDerivAt G (s/(s^2+(p+a)^2)^2) (p+a))
    (hminus : HasDerivAt G (s/(s^2+(p-a)^2)^2) (p-a)) :
    HasDerivAt (fun x => streamF G s a x) (density s a p*pa s a p) p := by
  have hpa := (hasDerivAt_id p).add_const a
  have hma := (hasDerivAt_id p).sub_const a
  have hden := ((hasDerivAt_const p (s^2)).add (hpa.pow 2)).mul
    ((hasDerivAt_const p (s^2)).add (hma.pow 2))
  have hd : D s a p ≠ 0 := ne_of_gt (denominator_positive s a p hs)
  have hquot := (hasDerivAt_const p (2*s*a)).div hden hd
  have hx := ((hplus.comp p hpa).sub (hminus.comp p hma)).add hquot
  convert hx using 1
  have h1 : s^2+(p+a)^2 ≠ 0 := ne_of_gt (show 0 < s^2+(p+a)^2 by positivity)
  have h2 : s^2+(p-a)^2 ≠ 0 := ne_of_gt (show 0 < s^2+(p-a)^2 by positivity)
  dsimp [density,pa,D]
  field_simp
  <;> ring

theorem invariant_along_guidance (F : ℝ → ℝ) (t r sa sp b u : ℝ)
    (hchain : HasDerivAt F ((r*sp)*(-b*u*sa)+(r*sa)*(b*u*sp)) t) :
    HasDerivAt F 0 t := by
  convert hchain using 1
  ring

end
end UnifiedBounce
