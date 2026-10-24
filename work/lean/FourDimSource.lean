import Mathlib

/-! Specified zero-temperature massless order-reduced AHS source.
The equations are reduced equations, not a quantum-state construction. -/
namespace FourDimSource
noncomputable section

def criticalCoupling : ℝ := 11/60

def A (xi h : ℝ) : ℝ :=
  (xi-11/60)*h^2+2*(xi-13/60)*h+5*(xi-83/300)
def Q (xi h : ℝ) : ℝ :=
  (xi-11/60)*h^2+6*(xi-11/60)*h+9*(xi-5/36)
def Qh (xi h : ℝ) : ℝ := 2*(xi-11/60)*h+6*(xi-11/60)

def timeSource (xi r h : ℝ) : ℝ := -(h-1)^2*A xi h/(8*h^2*r^4)
def pressure (xi r h : ℝ) : ℝ := (h-1)^2*Q xi h/(24*h^2*r^4)
def hRate (L xi r h : ℝ) : ℝ := -h*(h-1)/r+L*(h-1)^2*A xi h/r^3
def logFRate (L xi r h : ℝ) : ℝ :=
  (h-1)/r+L*(h-1)^2*Q xi h/(3*h*r^3)
def pressureRate (xi r h dh : ℝ) : ℝ :=
  dh*(h-1)*(2*Q xi h+h*(h-1)*Qh xi h)/(24*h^3*r^4)-4*pressure xi r h/r
def angular (xi r h dh g : ℝ) : ℝ :=
  pressure xi r h+r*pressureRate xi r h dh/2+r*g*(pressure xi r h-timeSource xi r h)/4

theorem time_source_original (xi r h : ℝ) :
    timeSource xi r h = (h-1)^2/(h^2*r^4)*
      ((h^2+6*h+33)/480-(xi-1/6)*(h^2+2*h+5)/8) := by
  unfold timeSource A
  ring

theorem pressure_original (xi r h : ℝ) :
    pressure xi r h = (h-1)^2/(h^2*r^4)*
      (-(h^2+6*h-15)/1440+(xi-1/6)*(h+3)^2/24) := by
  unfold pressure Q
  ring

theorem pressure_hasDerivAt (h : ℝ → ℝ) (xi r dh : ℝ)
    (hh : HasDerivAt h dh r) (hr : r ≠ 0) (hh0 : h r ≠ 0) :
    HasDerivAt (fun x => pressure xi x (h x)) (pressureRate xi r (h r) dh) r := by
  have hQ : HasDerivAt (fun x => Q xi (h x)) (Qh xi (h r)*dh) r := by
    have hq := (((hh.pow 2).const_mul (xi-11/60)).add (hh.const_mul (6*(xi-11/60)))).add_const (9*(xi-5/36))
    convert hq using 1 <;> dsimp [Q,Qh] <;> (try funext x) <;> (try dsimp) <;> ring
  have hn := ((hh.sub_const 1).pow 2).mul hQ
  have hd := ((hh.pow 2).mul ((hasDerivAt_id r).pow 4)).const_mul 24
  convert hn.div hd (by dsimp; positivity) using 1 <;>
    dsimp [pressure,pressureRate,Q,Qh] <;> (try funext x) <;> (try dsimp) <;>
    field_simp <;> ring

theorem full_source_conservation (h : ℝ → ℝ) (xi r dh g : ℝ)
    (hh : HasDerivAt h dh r) (hr : r ≠ 0) (hh0 : h r ≠ 0) :
    deriv (fun x => pressure xi x (h x)) r+
      g*(pressure xi r (h r)-timeSource xi r (h r))/2+
      2/r*(pressure xi r (h r)-angular xi r (h r) dh g)=0 := by
  rw [(pressure_hasDerivAt h xi r dh hh hr hh0).deriv]
  unfold angular
  field_simp
  <;> ring

theorem angular_source_unique (xi r h dh g p : ℝ) (hr : r ≠ 0)
    (hp : pressureRate xi r h dh+g*(pressure xi r h-timeSource xi r h)/2+
      2/r*(pressure xi r h-p)=0) : p=angular xi r h dh g := by
  unfold angular
  field_simp at hp ⊢
  nlinarith

theorem time_einstein_equation (L xi r h : ℝ) (hr : r ≠ 0) (hh : h ≠ 0) :
    (h*(1-h)-r*hRate L xi r h)/(h^2*r^2)=8*L*timeSource xi r h := by
  unfold hRate timeSource
  field_simp
  <;> ring

theorem radial_einstein_equation (L xi r h : ℝ) (hr : r ≠ 0) (hh : h ≠ 0) :
    (r*logFRate L xi r h+1-h)/(h*r^2)=8*L*pressure xi r h := by
  unfold logFRate pressure
  field_simp
  <;> ring

theorem critical_polynomial_decomposition (xi h : ℝ) :
    A xi h=(xi-criticalCoupling)*(h^2+2*h+5)-(h+7)/15 ∧
    Q xi h=(xi-criticalCoupling)*(h+3)^2+2/5 := by
  unfold A Q criticalCoupling
  constructor <;> ring

theorem interior_h_rate_negative (L xi r h : ℝ)
    (hL : 0 ≤ L) (hxi : xi ≤ criticalCoupling) (hr : 0 < r) (hh : 1 < h) :
    hRate L xi r h < 0 := by
  have ha := (critical_polynomial_decomposition xi h).1
  have hpoly : 0 < h^2+2*h+5 := by positivity
  have hq : (xi-criticalCoupling)*(h^2+2*h+5) ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg (sub_nonpos.mpr hxi) (le_of_lt hpoly)
  have hA : A xi h < 0 := by rw [ha]; linarith
  have hclass : -h*(h-1)/r < 0 := div_neg_of_neg_of_pos (by nlinarith) hr
  have hquant : L*(h-1)^2*A xi h/r^3 ≤ 0 :=
    div_nonpos_of_nonpos_of_nonneg (mul_nonpos_of_nonneg_of_nonpos (by positivity) (le_of_lt hA)) (by positivity)
  unfold hRate
  linarith

theorem minimal_conformal_subcritical : (0 : ℝ) < criticalCoupling ∧ (1/6 : ℝ) < criticalCoupling := by
  norm_num [criticalCoupling]

end
end FourDimSource
