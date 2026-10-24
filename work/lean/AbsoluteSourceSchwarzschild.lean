import Mathlib

/-! Leading large-mass Schwinger--DeWitt source and first-order metric
   response on Schwarzschild. All equalities below concern this truncation. -/
namespace AbsoluteSourceSchwarzschild
noncomputable section

def prefactor (M mu r : ℝ) : ℝ := M^2/(20160*Real.pi^2*mu^2*r^9)
def tt (M mu r : ℝ) : ℝ := prefactor M mu r*(2384*M-1080*r)
def rr (M mu r : ℝ) : ℝ := prefactor M mu r*(504*r-784*M)
def angular (M mu r : ℝ) : ℝ := prefactor M mu r*(3536*M-1512*r)
def metric (M r : ℝ) : ℝ := 1-2*M/r

def massShift (M mu G r : ℝ) : ℝ :=
  G*M^2*(149*M-81*r)/(1890*Real.pi*mu^2*r^6)
def lapseShift (M mu G r : ℝ) : ℝ :=
  -11*G*M^2/(210*Real.pi*mu^2*r^6)

theorem prefactor_hasDerivAt (M mu r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0) :
    HasDerivAt (prefactor M mu) (-9*M^2/(20160*Real.pi^2*mu^2*r^10)) r := by
  have hpi := Real.pi_ne_zero
  convert (hasDerivAt_const r (M^2)).div
    (((hasDerivAt_id r).pow 9).const_mul (20160*Real.pi^2*mu^2))
    (by positivity : 20160*Real.pi^2*mu^2*r^9 ≠ 0) using 1 <;>
    dsimp [prefactor] <;> field_simp <;> ring

theorem rr_hasDerivAt (M mu r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0) :
    HasDerivAt (rr M mu)
      (M^2*(7056*M-4032*r)/(20160*Real.pi^2*mu^2*r^10)) r := by
  have hpi := Real.pi_ne_zero
  convert (prefactor_hasDerivAt M mu r hmu hr).mul
    (((hasDerivAt_id r).const_mul 504).sub_const (784*M)) using 1 <;>
    dsimp [rr,prefactor] <;> field_simp <;> ring

theorem metric_hasDerivAt (M r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (metric M) (2*M/r^2) r := by
  convert (hasDerivAt_const r 1).sub
    ((hasDerivAt_const r (2*M)).div (hasDerivAt_id r) hr) using 1 <;>
    dsimp [metric] <;> field_simp <;> ring

theorem source_ward_identity (M mu r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0)
    (hhor : r ≠ 2*M) :
    deriv (rr M mu) r + deriv (metric M) r/(2*metric M r)*(rr M mu r-tt M mu r)
      +2/r*(rr M mu r-angular M mu r)=0 := by
  rw [(rr_hasDerivAt M mu r hmu hr).deriv,(metric_hasDerivAt M r hr).deriv]
  have hpi := Real.pi_ne_zero
  have hh : r-2*M ≠ 0 := sub_ne_zero.mpr hhor
  have he : metric M r=(r-2*M)/r := by unfold metric; field_simp <;> ring
  rw [he]
  unfold rr tt angular prefactor
  field_simp
  ring_nf
  have hnorm : -(M*2)+r ≠ 0 := by
    intro hz
    apply hh
    linarith
  field_simp [hnorm]
  <;> ring

theorem mass_shift_hasDerivAt (M mu G r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0) :
    HasDerivAt (massShift M mu G) (-4*Real.pi*G*r^2*tt M mu r) r := by
  have hpi := Real.pi_ne_zero
  have hn : HasDerivAt (fun x : ℝ => G*M^2*(149*M-81*x)) (-81*G*M^2) r := by
    convert (((hasDerivAt_const r (149*M)).sub ((hasDerivAt_id r).const_mul 81)).const_mul
      (G*M^2)) using 1 <;> dsimp <;> ring
  have hd := ((hasDerivAt_id r).pow 6).const_mul (1890*Real.pi*mu^2)
  convert hn.div hd (by positivity : 1890*Real.pi*mu^2*r^6 ≠ 0) using 1 <;>
    dsimp [massShift,tt,prefactor] <;> field_simp <;> ring

theorem radial_stress_difference (M mu r : ℝ) :
    rr M mu r-tt M mu r=1584*prefactor M mu r*(r-2*M) := by
  unfold rr tt
  ring

theorem quotient_extends_to_horizon (M mu r : ℝ) (hr : r ≠ 0) (hhor : r ≠ 2*M) :
    (rr M mu r-tt M mu r)/metric M r=1584*r*prefactor M mu r := by
  rw [radial_stress_difference]
  unfold metric
  field_simp
  <;> ring

theorem lapse_shift_hasDerivAt (M mu G r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0) :
    HasDerivAt (lapseShift M mu G)
      (4*Real.pi*G*r*(1584*r*prefactor M mu r)) r := by
  have hpi := Real.pi_ne_zero
  convert (hasDerivAt_const r (-11*G*M^2)).div
    (((hasDerivAt_id r).pow 6).const_mul (210*Real.pi*mu^2))
    (by positivity : 210*Real.pi*mu^2*r^6 ≠ 0) using 1 <;>
    dsimp [lapseShift,prefactor] <;> field_simp <;> ring

theorem lapse_einstein_equation (M mu G r : ℝ) (hmu : mu ≠ 0) (hr : r ≠ 0)
    (hhor : r ≠ 2*M) :
    deriv (lapseShift M mu G) r=4*Real.pi*G*r*(rr M mu r-tt M mu r)/metric M r := by
  rw [(lapse_shift_hasDerivAt M mu G r hmu hr).deriv]
  rw [mul_div_assoc,quotient_extends_to_horizon M mu r hr hhor]

theorem horizon_shift_value (M mu G : ℝ) (hM : M ≠ 0) (hmu : mu ≠ 0) :
    2*massShift M mu G (2*M)= -13*G/(60480*Real.pi*mu^2*M^3) := by
  have hpi := Real.pi_ne_zero
  unfold massShift
  field_simp
  ring

theorem horizon_shift_negative (M mu G : ℝ) (hM : 0 < M) (hmu : 0 < mu) (hG : 0 < G) :
    2*massShift M mu G (2*M)<0 := by
  rw [horizon_shift_value M mu G (ne_of_gt hM) (ne_of_gt hmu)]
  exact div_neg_of_neg_of_pos (by nlinarith) (by positivity)

theorem horizon_mass_derivative (M mu G : ℝ) (hM : M ≠ 0) (hmu : mu ≠ 0) :
    deriv (massShift M mu G) (2*M)= -G/(2880*Real.pi*mu^2*M^4) := by
  rw [(mass_shift_hasDerivAt M mu G (2*M) hmu (mul_ne_zero (by norm_num) hM)).deriv]
  have hpi := Real.pi_ne_zero
  unfold tt prefactor
  field_simp <;> ring

theorem horizon_lapse (M mu G : ℝ) (hM : M ≠ 0) (hmu : mu ≠ 0) :
    lapseShift M mu G (2*M)= -11*G/(13440*Real.pi*mu^2*M^4) := by
  have hpi := Real.pi_ne_zero
  unfold lapseShift
  field_simp <;> ring

theorem surface_gravity_shift (M mu G : ℝ) (hM : M ≠ 0) (hmu : mu ≠ 0) :
    lapseShift M mu G (2*M)/(4*M)-deriv (massShift M mu G) (2*M)/(2*M)
      -massShift M mu G (2*M)/(4*M^2)= -G/(241920*Real.pi*mu^2*M^5) := by
  rw [horizon_mass_derivative M mu G hM hmu,horizon_lapse M mu G hM hmu]
  have hpi := Real.pi_ne_zero
  unfold massShift
  field_simp <;> ring

end
end AbsoluteSourceSchwarzschild
