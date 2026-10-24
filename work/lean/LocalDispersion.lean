import Mathlib

/-! Finite Fourier-symbol and independent-sector identities.
The variables below are real Fourier components, not a continuum regulator. -/
namespace LocalDispersion
noncomputable section

def dotMass (m k1 k2 k3 p1 p2 p3 : ℝ) : ℝ :=
  m^2+k1*(k1+p1)+k2*(k2+p2)+k3*(k3+p3)

def crossSquare (k1 k2 k3 p1 p2 p3 : ℝ) : ℝ :=
  (k2*p3-k3*p2)^2+(k3*p1-k1*p3)^2+(k1*p2-k2*p1)^2

theorem independent_sector_coefficients (a b c d z : ℝ) :
    (∀ x y : ℝ, a*b*x+c*d*y=z*(x+y)) ↔ a*b=z ∧ c*d=z := by
  constructor
  · intro h
    have h1 := h 1 0
    have h2 := h 0 1
    constructor <;> nlinarith
  · rintro ⟨h1,h2⟩ x y
    rw [h1,h2]
    ring

theorem classical_sector_coefficients (a b c d : ℝ) :
    (∀ x y : ℝ, a*b*x+c*d*y=x+y) ↔ a*b=1 ∧ c*d=1 := by
  simpa using independent_sector_coefficients a b c d 1

theorem dispersion_gap_product
    (m k1 k2 k3 p1 p2 p3 w z : ℝ)
    (hw : w^2=m^2+k1^2+k2^2+k3^2)
    (hz : z^2=m^2+(k1+p1)^2+(k2+p2)^2+(k3+p3)^2) :
    (w*z-dotMass m k1 k2 k3 p1 p2 p3)*
      (w*z+dotMass m k1 k2 k3 p1 p2 p3)=
      m^2*(p1^2+p2^2+p3^2)+crossSquare k1 k2 k3 p1 p2 p3 := by
  calc
    _ = w^2*z^2-(dotMass m k1 k2 k3 p1 p2 p3)^2 := by ring
    _ = _ := by rw [hw,hz]; dsimp [dotMass,crossSquare]; ring

theorem momentum_square_positive (p1 p2 p3 : ℝ)
    (hp : p1 ≠ 0 ∨ p2 ≠ 0 ∨ p3 ≠ 0) : 0 < p1^2+p2^2+p3^2 := by
  rcases hp with h | h | h
  · nlinarith [sq_pos_of_ne_zero h, sq_nonneg p2, sq_nonneg p3]
  · nlinarith [sq_pos_of_ne_zero h, sq_nonneg p1, sq_nonneg p3]
  · nlinarith [sq_pos_of_ne_zero h, sq_nonneg p1, sq_nonneg p2]

theorem dispersion_gap_positive
    (m k1 k2 k3 p1 p2 p3 w z : ℝ)
    (hm : 0 < m) (hwpos : 0 < w) (hzpos : 0 < z)
    (hp : p1 ≠ 0 ∨ p2 ≠ 0 ∨ p3 ≠ 0)
    (hw : w^2=m^2+k1^2+k2^2+k3^2)
    (hz : z^2=m^2+(k1+p1)^2+(k2+p2)^2+(k3+p3)^2) :
    0 < w*z-dotMass m k1 k2 k3 p1 p2 p3 := by
  have hidentity := dispersion_gap_product m k1 k2 k3 p1 p2 p3 w z hw hz
  have hmass := mul_pos (sq_pos_of_pos hm) (momentum_square_positive p1 p2 p3 hp)
  have hcross : 0 ≤ crossSquare k1 k2 k3 p1 p2 p3 := by
    dsimp [crossSquare]
    positivity
  have hprod : 0 < (w*z-dotMass m k1 k2 k3 p1 p2 p3)*
      (w*z+dotMass m k1 k2 k3 p1 p2 p3) := by rw [hidentity]; linarith
  have hwz := mul_pos hwpos hzpos
  rcases (mul_pos_iff.mp hprod) with h | h
  · exact h.1
  · linarith [h.1,h.2]

theorem dispersion_denominator_positive
    (m k1 k2 k3 p1 p2 p3 w z : ℝ)
    (hm : 0 < m) (hwpos : 0 < w) (hzpos : 0 < z)
    (hp : p1 ≠ 0 ∨ p2 ≠ 0 ∨ p3 ≠ 0)
    (hw : w^2=m^2+k1^2+k2^2+k3^2)
    (hz : z^2=m^2+(k1+p1)^2+(k2+p2)^2+(k3+p3)^2) :
    0 < w*z+dotMass m k1 k2 k3 p1 p2 p3 := by
  have hgap := dispersion_gap_positive m k1 k2 k3 p1 p2 p3 w z hm hwpos hzpos hp hw hz
  have hidentity := dispersion_gap_product m k1 k2 k3 p1 p2 p3 w z hw hz
  have hmass := mul_pos (sq_pos_of_pos hm) (momentum_square_positive p1 p2 p3 hp)
  have hcross : 0 ≤ crossSquare k1 k2 k3 p1 p2 p3 := by dsimp [crossSquare]; positivity
  have hprod : 0 < (w*z-dotMass m k1 k2 k3 p1 p2 p3)*
      (w*z+dotMass m k1 k2 k3 p1 p2 p3) := by rw [hidentity]; linarith
  exact (mul_pos_iff_of_pos_left hgap).mp hprod

theorem dispersion_gap_rationalized
    (m k1 k2 k3 p1 p2 p3 w z : ℝ)
    (hm : 0 < m) (hwpos : 0 < w) (hzpos : 0 < z)
    (hp : p1 ≠ 0 ∨ p2 ≠ 0 ∨ p3 ≠ 0)
    (hw : w^2=m^2+k1^2+k2^2+k3^2)
    (hz : z^2=m^2+(k1+p1)^2+(k2+p2)^2+(k3+p3)^2) :
    w*z-dotMass m k1 k2 k3 p1 p2 p3=
      (m^2*(p1^2+p2^2+p3^2)+crossSquare k1 k2 k3 p1 p2 p3)/
        (w*z+dotMass m k1 k2 k3 p1 p2 p3) := by
  apply (eq_div_iff (ne_of_gt (dispersion_denominator_positive
    m k1 k2 k3 p1 p2 p3 w z hm hwpos hzpos hp hw hz))).mpr
  exact dispersion_gap_product m k1 k2 k3 p1 p2 p3 w z hw hz

theorem homogeneous_lapse_gap (m k1 k2 k3 w : ℝ)
    (hw : w^2=m^2+k1^2+k2^2+k3^2) :
    w*w-dotMass m k1 k2 k3 0 0 0=0 := by
  dsimp [dotMass]
  nlinarith

theorem massive_nonzero_mode_witness :
    0 < (1:ℝ)*Real.sqrt 2-dotMass 1 0 0 0 1 0 0 := by
  apply dispersion_gap_positive 1 0 0 0 1 0 0 1 (Real.sqrt 2)
  · norm_num
  · norm_num
  · positivity
  · exact Or.inl (by norm_num)
  · norm_num
  · norm_num

end
end LocalDispersion
