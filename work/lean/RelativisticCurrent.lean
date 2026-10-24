import Mathlib

/-!
Finite-dimensional two-particle Dirac current in 1+1 dimensions.
The basis diagonalizes alpha = diag(1,-1); a,b,c,d are its four squared
spinor amplitudes, ordered ++,+-,-+,--. Normals are future timelike
COVECTORS; hence n0 > |n1|, independently of the sign of n1.
These results prove pointwise algebra, not existence of PDE solutions.
-/

namespace RelativisticCurrent

def j00 (a b c d : ℝ) : ℝ := a + b + c + d
def j10 (a b c d : ℝ) : ℝ := a + b - c - d
def j01 (a b c d : ℝ) : ℝ := a - b + c - d
def j11 (a b c d : ℝ) : ℝ := a - b - c + d

def density (a b c d n0 n1 m0 m1 : ℝ) : ℝ :=
  j00 a b c d * n0 * m0 + j10 a b c d * n1 * m0 +
    j01 a b c d * n0 * m1 + j11 a b c d * n1 * m1

def firstTime (a b c d m0 m1 : ℝ) : ℝ :=
  j00 a b c d * m0 + j01 a b c d * m1

def firstSpace (a b c d m0 m1 : ℝ) : ℝ :=
  j10 a b c d * m0 + j11 a b c d * m1

theorem normal_factors_positive (n0 n1 : ℝ) (hn : |n1| < n0) :
    0 < n0 + n1 ∧ 0 < n0 - n1 := by
  have h := abs_lt.mp hn
  constructor <;> linarith

theorem density_factorization (a b c d n0 n1 m0 m1 : ℝ) :
    density a b c d n0 n1 m0 m1 =
      a * (n0 + n1) * (m0 + m1) + b * (n0 + n1) * (m0 - m1) +
      c * (n0 - n1) * (m0 + m1) + d * (n0 - n1) * (m0 - m1) := by
  unfold density j00 j10 j01 j11
  ring

theorem four_weighted_positive (a b c d A B C D : ℝ)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d)
    (hs : 0 < a + b + c + d)
    (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D) :
    0 < a * A + b * B + c * C + d * D := by
  have hpos : 0 < a ∨ 0 < b ∨ 0 < c ∨ 0 < d := by
    by_contra h
    push Not at h
    linarith
  have haA := mul_nonneg ha (le_of_lt hA)
  have hbB := mul_nonneg hb (le_of_lt hB)
  have hcC := mul_nonneg hc (le_of_lt hC)
  have hdD := mul_nonneg hd (le_of_lt hD)
  rcases hpos with h | h | h | h
  · linarith [mul_pos h hA]
  · linarith [mul_pos h hB]
  · linarith [mul_pos h hC]
  · linarith [mul_pos h hD]

theorem density_positive_weights (a b c d n0 n1 m0 m1 : ℝ)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d)
    (hs : 0 < a + b + c + d) (hn : |n1| < n0) (hm : |m1| < m0) :
    0 < density a b c d n0 n1 m0 m1 := by
  obtain ⟨hnp, hnm⟩ := normal_factors_positive n0 n1 hn
  obtain ⟨hmp, hmm⟩ := normal_factors_positive m0 m1 hm
  rw [density_factorization]
  simpa only [mul_assoc] using
    four_weighted_positive a b c d ((n0+n1)*(m0+m1)) ((n0+n1)*(m0-m1))
      ((n0-n1)*(m0+m1)) ((n0-n1)*(m0-m1)) ha hb hc hd hs
      (mul_pos hnp hmp) (mul_pos hnp hmm) (mul_pos hnm hmp) (mul_pos hnm hmm)

theorem spinor_total_positive (a b c d : ℂ)
    (hne : a ≠ 0 ∨ b ≠ 0 ∨ c ≠ 0 ∨ d ≠ 0) :
    0 < Complex.normSq a + Complex.normSq b + Complex.normSq c + Complex.normSq d := by
  have ha := Complex.normSq_nonneg a
  have hb := Complex.normSq_nonneg b
  have hc := Complex.normSq_nonneg c
  have hd := Complex.normSq_nonneg d
  rcases hne with h | h | h | h
  · linarith [Complex.normSq_pos.mpr h]
  · linarith [Complex.normSq_pos.mpr h]
  · linarith [Complex.normSq_pos.mpr h]
  · linarith [Complex.normSq_pos.mpr h]

theorem spinor_hypersurface_density_positive (a b c d : ℂ) (n0 n1 m0 m1 : ℝ)
    (hne : a ≠ 0 ∨ b ≠ 0 ∨ c ≠ 0 ∨ d ≠ 0)
    (hn : |n1| < n0) (hm : |m1| < m0) :
    0 < density (Complex.normSq a) (Complex.normSq b) (Complex.normSq c)
      (Complex.normSq d) n0 n1 m0 m1 := by
  exact density_positive_weights _ _ _ _ _ _ _ _
    (Complex.normSq_nonneg a) (Complex.normSq_nonneg b)
    (Complex.normSq_nonneg c) (Complex.normSq_nonneg d)
    (spinor_total_positive a b c d hne) hn hm

theorem first_current_future_causal (a b c d m0 m1 : ℝ)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d)
    (hm : |m1| ≤ m0) :
    |firstSpace a b c d m0 m1| ≤ firstTime a b c d m0 m1 := by
  obtain ⟨hlo, hhi⟩ := abs_le.mp hm
  have hp : 0 ≤ m0 + m1 := by linarith
  have hn : 0 ≤ m0 - m1 := by linarith
  have ha' := mul_nonneg ha hp
  have hb' := mul_nonneg hb hn
  have hc' := mul_nonneg hc hp
  have hd' := mul_nonneg hd hn
  apply abs_le.mpr
  unfold firstSpace firstTime j00 j10 j01 j11
  constructor <;> nlinarith

theorem first_current_time_positive (a b c d m0 m1 : ℝ)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d)
    (hs : 0 < a + b + c + d) (hm : |m1| < m0) :
    0 < firstTime a b c d m0 m1 := by
  have h := density_positive_weights a b c d 1 0 m0 m1
    ha hb hc hd hs (by norm_num) hm
  simpa [density, firstTime] using h

theorem spinor_guidance_causal (a b c d : ℂ) (m0 m1 : ℝ)
    (hne : a ≠ 0 ∨ b ≠ 0 ∨ c ≠ 0 ∨ d ≠ 0) (hm : |m1| < m0) :
    0 < firstTime (Complex.normSq a) (Complex.normSq b) (Complex.normSq c)
          (Complex.normSq d) m0 m1 ∧
    |firstSpace (Complex.normSq a) (Complex.normSq b) (Complex.normSq c)
          (Complex.normSq d) m0 m1 /
      firstTime (Complex.normSq a) (Complex.normSq b) (Complex.normSq c)
          (Complex.normSq d) m0 m1| ≤ 1 := by
  have ht := first_current_time_positive _ _ _ _ m0 m1
    (Complex.normSq_nonneg a) (Complex.normSq_nonneg b)
    (Complex.normSq_nonneg c) (Complex.normSq_nonneg d)
    (spinor_total_positive a b c d hne) hm
  have hc := first_current_future_causal _ _ _ _ m0 m1
    (Complex.normSq_nonneg a) (Complex.normSq_nonneg b)
    (Complex.normSq_nonneg c) (Complex.normSq_nonneg d) (le_of_lt hm)
  constructor
  · exact ht
  · rw [abs_div, abs_of_pos ht]
    exact (div_le_one ht).mpr hc

theorem second_current_is_first_transposed (a b c d n0 n1 : ℝ) :
    j00 a b c d * n0 + j10 a b c d * n1 = firstTime a c b d n0 n1 ∧
    j01 a b c d * n0 + j11 a b c d * n1 = firstSpace a c b d n0 n1 := by
  simp only [firstTime, firstSpace, j00, j10, j01, j11]
  constructor <;> ring

theorem current_minkowski_square (a b c d m0 m1 : ℝ) :
    (firstTime a b c d m0 m1)^2 - (firstSpace a b c d m0 m1)^2 =
      4 * (a*(m0+m1) + b*(m0-m1)) * (c*(m0+m1) + d*(m0-m1)) := by
  unfold firstTime firstSpace j00 j10 j01 j11
  ring

end RelativisticCurrent
