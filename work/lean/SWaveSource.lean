import Mathlib

/-! Specified dilaton-coupled static source in conformal radial coordinate.
The common positive normalization is omitted. The angular component is
fixed by the four-dimensional radial conservation equation. -/
namespace SWaveSource
noncomputable section

def factor (rho r : ℝ) : ℝ := Real.exp (-2*rho)/r^2
def energy (rho a b c d : ℝ) : ℝ := 2*d-b^2+6*b*a+(6*rho+3)*a^2-3*c
def pressure (rho a b c : ℝ) : ℝ := -b^2+6*b*a+(6*rho-3)*a^2+3*c
def angular (rho a b c d e : ℝ) : ℝ :=
  -(3/2 : ℝ)*(e/a-2*c+2*d+4*b*a+4*rho*c)
def pressureRate (rho a b c d e : ℝ) : ℝ :=
  -2*b*d+6*d*a+6*b*c+6*b*a^2+(12*rho-6)*a*c+3*e

theorem factor_hasDerivAt (rho r : ℝ → ℝ) (x a b : ℝ)
    (hρ : HasDerivAt rho b x) (hr : HasDerivAt r (-a*r x) x) (hr0 : r x ≠ 0) :
    HasDerivAt (fun y => factor (rho y) (r y))
      ((2*a-2*b)*factor (rho x) (r x)) x := by
  convert ((hρ.const_mul (-2)).exp).div (hr.pow 2) (pow_ne_zero 2 hr0) using 1 <;>
    dsimp [factor] <;> field_simp <;> ring

theorem pressure_hasDerivAt (rho a b c : ℝ → ℝ) (x d e : ℝ)
    (hρ : HasDerivAt rho (b x) x) (ha : HasDerivAt a (c x) x)
    (hb : HasDerivAt b d x) (hc : HasDerivAt c e x) :
    HasDerivAt (fun y => pressure (rho y) (a y) (b y) (c y))
      (pressureRate (rho x) (a x) (b x) (c x) d e) x := by
  have h := (((hb.pow 2).neg.add ((hb.mul ha).const_mul 6)).add
    (((hρ.const_mul 6).sub_const 3).mul (ha.pow 2))).add (hc.const_mul 3)
  convert h using 1 <;> dsimp [pressure, pressureRate] <;> (try funext y) <;> (try dsimp) <;> ring

theorem reduced_ward_identity (rho a b c d e : ℝ) (ha : a ≠ 0) :
    pressureRate rho a b c d e+b*(energy rho a b c d-pressure rho a b c)+
      2*a*angular rho a b c d e=0 := by
  unfold pressureRate energy pressure angular
  field_simp
  <;> ring

theorem full_source_conservation (rho r a b c : ℝ → ℝ) (x d e : ℝ)
    (hρ : HasDerivAt rho (b x) x) (hr : HasDerivAt r (-(a x)*r x) x)
    (ha : HasDerivAt a (c x) x) (hb : HasDerivAt b d x) (hc : HasDerivAt c e x)
    (hr0 : r x ≠ 0) (ha0 : a x ≠ 0) :
    deriv (fun y => factor (rho y) (r y)*pressure (rho y) (a y) (b y) (c y)) x+
      b x*(factor (rho x) (r x)*energy (rho x) (a x) (b x) (c x) d+
        factor (rho x) (r x)*pressure (rho x) (a x) (b x) (c x))-
      2*a x*(factor (rho x) (r x)*pressure (rho x) (a x) (b x) (c x)-
        factor (rho x) (r x)*angular (rho x) (a x) (b x) (c x) d e)=0 := by
  have hF := factor_hasDerivAt rho r x (a x) (b x) hρ hr hr0
  have hp := pressure_hasDerivAt rho a b c x d e hρ ha hb hc
  have hd : deriv (fun y => factor (rho y) (r y)*pressure (rho y) (a y) (b y) (c y)) x =
      ((2*a x-2*b x)*factor (rho x) (r x))*pressure (rho x) (a x) (b x) (c x)+
      factor (rho x) (r x)*pressureRate (rho x) (a x) (b x) (c x) d e := by
    simpa using (hF.mul hp).deriv
  rw [hd]
  have hw := reduced_ward_identity (rho x) (a x) (b x) (c x) d e ha0
  have hz := congrArg (fun y : ℝ => factor (rho x) (r x)*y) hw
  dsimp at hz
  nlinarith [hz]

theorem angular_source_unique (rho a b c d e p : ℝ) (ha : a ≠ 0)
    (hward : pressureRate rho a b c d e+
      b*(energy rho a b c d-pressure rho a b c)+2*a*p=0) :
    p=angular rho a b c d e := by
  have hw := reduced_ward_identity rho a b c d e ha
  have hh : 2*a*(p-angular rho a b c d e)=0 := by linarith
  exact (mul_eq_zero.mp hh).resolve_left (mul_ne_zero (by norm_num) ha) |> sub_eq_zero.mp

end
end SWaveSource
