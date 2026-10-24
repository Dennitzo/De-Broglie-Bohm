import Mathlib

/-! Radiation/Schutz minisuperspace bounce with lapse N=a. The scale factor
is a quantum configuration variable in this truncation. Derivatives below
are with respect to conformal fluid time T; H=(da/dT)/a^2 is cosmic H. -/

namespace GravityBounce

noncomputable section

def scale (ab tau T : ℝ) : ℝ := ab*Real.sqrt (1+T^2/tau^2)

theorem radicand_positive (tau T : ℝ) (htau : 0 < tau) :
    0 < 1+T^2/tau^2 := by positivity

theorem scale_minimum (ab tau T : ℝ) (hab : 0 < ab) (htau : 0 < tau) :
    ab ≤ scale ab tau T := by
  have hrad := radicand_positive tau T htau
  have hs := Real.sq_sqrt (le_of_lt hrad)
  have hn := Real.sqrt_nonneg (1+T^2/tau^2)
  have hq : 0 ≤ T^2/tau^2 := by positivity
  have hroot : 1 ≤ Real.sqrt (1+T^2/tau^2) := by nlinarith
  unfold scale
  nlinarith

theorem scale_positive (ab tau T : ℝ) (hab : 0 < ab) (htau : 0 < tau) :
    0 < scale ab tau T := lt_of_lt_of_le hab (scale_minimum ab tau T hab htau)

theorem scale_squared (ab tau T : ℝ) (htau : 0 < tau) :
    (scale ab tau T)^2=ab^2*(1+T^2/tau^2) := by
  unfold scale
  rw [mul_pow, Real.sq_sqrt (le_of_lt (radicand_positive tau T htau))]

theorem scale_hasDerivAt (ab tau T : ℝ) (htau : 0 < tau) :
    HasDerivAt (scale ab tau) (ab*T/(tau^2*Real.sqrt (1+T^2/tau^2))) T := by
  have hi := (((hasDerivAt_id T).pow 2).div_const (tau^2)).const_add 1
  have h := (hi.sqrt (ne_of_gt (radicand_positive tau T htau))).const_mul ab
  convert h using 1 <;> dsimp [scale] <;> field_simp <;> ring

theorem scale_momentum_identity (ab tau T : ℝ) (htau : 0 < tau) :
    deriv (scale ab tau) T*scale ab tau T=ab^2*T/tau^2 := by
  have hn : Real.sqrt (1+T^2/tau^2) ≠ 0 :=
    ne_of_gt (Real.sqrt_pos.mpr (radicand_positive tau T htau))
  rw [(scale_hasDerivAt ab tau T htau).deriv]
  unfold scale
  field_simp
  <;> ring

theorem scale_bohm_guidance (ab tau T : ℝ) (htau : 0 < tau) :
    deriv (scale ab tau) T=scale ab tau T*T/(tau^2+T^2) := by
  have ht : tau ≠ 0 := ne_of_gt htau
  have hr : Real.sqrt (1+T^2/tau^2) ≠ 0 :=
    ne_of_gt (Real.sqrt_pos.mpr (radicand_positive tau T htau))
  have hd : tau^2+T^2 ≠ 0 := by positivity
  have hs : tau^2*(Real.sqrt (1+T^2/tau^2))^2=tau^2+T^2 := by
    rw [Real.sq_sqrt (le_of_lt (radicand_positive tau T htau))]
    field_simp
  rw [(scale_hasDerivAt ab tau T htau).deriv]
  unfold scale
  generalize hroot : Real.sqrt (1+T^2/tau^2)=s at *
  field_simp
  nlinarith [congrArg (fun z : ℝ => ab*T*z) hs]

theorem modified_friedmann (ab tau T : ℝ) (hab : 0 < ab) (htau : 0 < tau) :
    (deriv (scale ab tau) T/(scale ab tau T)^2)^2=
      (ab^2/tau^2)/(scale ab tau T)^4-(ab^4/tau^2)/(scale ab tau T)^6 := by
  have ha : scale ab tau T ≠ 0 := ne_of_gt (scale_positive ab tau T hab htau)
  have ht : tau ≠ 0 := ne_of_gt htau
  have hm := scale_momentum_identity ab tau T htau
  have hs := scale_squared ab tau T htau
  have hcore : (deriv (scale ab tau) T*scale ab tau T)^2=
      (ab^2/tau^2)*(scale ab tau T)^2-ab^4/tau^2 := by
    rw [hm, hs]
    field_simp
    <;> ring
  field_simp at hcore ⊢
  nlinarith [hcore]

theorem bounce_nec_negative (ab tau : ℝ) (hab : 0 < ab) (htau : 0 < tau) :
    (4/3:ℝ)*(ab^2/tau^2)/ab^4-2*(ab^4/tau^2)/ab^6 < 0 := by
  have ha : ab ≠ 0 := ne_of_gt hab
  have ht : tau ≠ 0 := ne_of_gt htau
  have he : (4/3:ℝ)*(ab^2/tau^2)/ab^4-2*(ab^4/tau^2)/ab^6=
      -2/(3*ab^2*tau^2) := by field_simp; ring
  rw [he]
  have hd : 0 < 3*ab^2*tau^2 := by positivity
  exact div_neg_of_neg_of_pos (by norm_num) hd

theorem negative_kinetic_hamilton_jacobi (tau T a : ℝ) (htau : 0 < tau) :
    (tau/(2*(tau^2+T^2))-a^2*(tau^2-T^2)/(tau^2+T^2)^2)-
      (-2*a*T/(tau^2+T^2))^2/4+
      (-tau/(2*(tau^2+T^2))+a^2*tau^2/(tau^2+T^2)^2)=0 := by
  have hd : tau^2+T^2 ≠ 0 := by positivity
  field_simp
  <;> ring

end

end GravityBounce
