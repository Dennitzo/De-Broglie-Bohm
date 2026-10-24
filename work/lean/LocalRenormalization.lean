import Mathlib

/-! Local-lapse ratios and thermodynamic consequences of positive Casimir sums.
The Bessel sum derivative and positivity are explicit hypotheses; this file
does not formalize convergence of the infinite mode sum. -/
namespace LocalRenormalization
noncomputable section

theorem local_lapse_ratio (m k w : ℝ) (hm : 0 < m) (hk : 0 < k)
    (hw : w^2=m^2+k^2) (hwpos : 0 < w) :
    m*(w-m)/k^2=m/(w+m) := by
  have hk0 := ne_of_gt hk
  have hden : w+m ≠ 0 := ne_of_gt (add_pos hwpos hm)
  field_simp
  nlinarith [hw]

theorem local_lapse_ratio_bounds (m k w : ℝ) (hm : 0 < m) (hk : 0 < k)
    (hw : w^2=m^2+k^2) (hwpos : 0 < w) :
    0 < m/(w+m) ∧ m/(w+m)<1/2 := by
  have hwm : m < w := by nlinarith [sq_pos_of_pos hk]
  have hden := add_pos hwpos hm
  constructor
  · exact div_pos hm hden
  · apply (div_lt_iff₀ hden).mpr
    linarith

theorem local_lapse_ratios_distinct (m w z : ℝ) (hm : 0 < m)
    (hw : 0 < w) (hz : 0 < z) (hwz : w ≠ z) :
    m/(w+m) ≠ m/(z+m) := by
  intro heq
  have hd1 := ne_of_gt (add_pos hw hm)
  have hd2 := ne_of_gt (add_pos hz hm)
  have hmul := (div_eq_div_iff hd1 hd2).mp heq
  have he : m*(z-w)=0 := by nlinarith [hmul]
  have hzme : z-w=0 := (mul_eq_zero.mp he).resolve_left (ne_of_gt hm)
  exact hwz (by linarith)

theorem no_common_gradient_repair (m w z delta : ℝ) (hm : 0 < m)
    (hw : 0 < w) (hz : 0 < z) (hwz : w ≠ z) :
    ¬(m/(w+m)+delta=1 ∧ m/(z+m)+delta=1) := by
  rintro ⟨h1,h2⟩
  exact local_lapse_ratios_distinct m w z hm hw hz hwz (by linarith)

def density (c s l : ℝ) : ℝ := -c*s/l^2
def pressure (c m s t l : ℝ) : ℝ := -c*(s+m*l*t)/(3*l^2)

theorem casimir_density_derivative (s : ℝ → ℝ) (c m l t : ℝ)
    (hl : l ≠ 0) (hs : HasDerivAt s (-m*t-2*s l/l) l) :
    HasDerivAt (fun x => density c (s x) x)
      (c*(4*s l+m*l*t)/l^3) l := by
  convert (hs.const_mul (-c)).div ((hasDerivAt_id l).pow 2) (pow_ne_zero 2 hl)
    using 1 <;> dsimp [density] <;> field_simp <;> ring

theorem casimir_thermodynamic_pressure (s : ℝ → ℝ) (c m l t : ℝ)
    (hl : l ≠ 0) (hs : HasDerivAt s (-m*t-2*s l/l) l) :
    -density c (s l) l-l/3*deriv (fun x => density c (s x) x) l=
      pressure c m (s l) t l := by
  rw [(casimir_density_derivative s c m l t hl hs).deriv]
  dsimp [density,pressure]
  field_simp
  <;> ring

theorem casimir_enthalpy (c m s t l : ℝ) :
    density c s l+pressure c m s t l= -c*(4*s+m*l*t)/(3*l^2) := by
  dsimp [density,pressure]
  ring

theorem casimir_enthalpy_negative (c m s t l : ℝ)
    (hc : 0 < c) (hm : 0 ≤ m) (hs : 0 < s) (ht : 0 ≤ t) (hl : 0 < l) :
    density c s l+pressure c m s t l < 0 := by
  rw [casimir_enthalpy]
  have hsum : 0 < 4*s+m*l*t := by positivity
  apply div_neg_of_neg_of_pos
  · exact mul_neg_of_neg_of_pos (neg_neg_of_pos hc) hsum
  · positivity

theorem cosmological_shift_enthalpy (rho p lambda : ℝ) :
    (rho+lambda)+(p-lambda)=rho+p := by ring

theorem flat_static_obstruction (rho p lambda : ℝ) (h : rho+p < 0) :
    ¬(rho+lambda=0 ∧ p-lambda=0) := by
  rintro ⟨hr,hp⟩
  linarith

theorem casimir_flat_static_obstruction (c m s t l lambda : ℝ)
    (hc : 0 < c) (hm : 0 ≤ m) (hs : 0 < s) (ht : 0 ≤ t) (hl : 0 < l) :
    ¬(density c s l+lambda=0 ∧ pressure c m s t l-lambda=0) :=
  flat_static_obstruction _ _ lambda (casimir_enthalpy_negative c m s t l hc hm hs ht hl)

end
end LocalRenormalization
