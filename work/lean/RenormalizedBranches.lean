import Mathlib

/-! Exact branches for the specified zero-counterterm FLRW source.
All differentiations below state their time-evolution hypotheses.
The module does not construct a renormalized quantum field state. -/

namespace RenormalizedBranches

noncomputable section

def friedmannPolynomial (alpha lam radiation y : ℝ) : ℝ :=
  alpha*y^2-y+lam+radiation
def discriminant (alpha lam radiation : ℝ) : ℝ := 1-4*alpha*(lam+radiation)
def upperBranch (alpha lam radiation : ℝ) : ℝ :=
  (1+Real.sqrt (discriminant alpha lam radiation))/(2*alpha)
def lowerBranch (alpha lam radiation : ℝ) : ℝ :=
  (1-Real.sqrt (discriminant alpha lam radiation))/(2*alpha)

theorem constraint_polynomial_equivalence (c A C a H Lambda : ℝ) :
    (H^2=c*(A*H^4+C/a^4)+Lambda/3) ↔
      friedmannPolynomial (c*A) (Lambda/3) (c*C/a^4) (H^2)=0 := by
  unfold friedmannPolynomial
  have he : c*A*(H^2)^2-H^2+Lambda/3+c*C/a^4=
      c*(A*H^4+C/a^4)+Lambda/3-H^2 := by ring
  rw [he]
  constructor <;> intro hh <;> linarith

theorem root_discriminant_square (alpha lam radiation y : ℝ)
    (hroot : friedmannPolynomial alpha lam radiation y=0) :
    (2*alpha*y-1)^2=discriminant alpha lam radiation := by
  unfold friedmannPolynomial at hroot
  unfold discriminant
  have hm := congrArg (fun x : ℝ => 4*alpha*x) hroot
  dsimp at hm
  nlinarith [hm]

theorem real_root_discriminant_nonnegative (alpha lam radiation y : ℝ)
    (hroot : friedmannPolynomial alpha lam radiation y=0) :
    0 ≤ discriminant alpha lam radiation := by
  rw [← root_discriminant_square alpha lam radiation y hroot]
  exact sq_nonneg _

theorem upper_branch_factor (alpha lam radiation : ℝ) (ha : alpha ≠ 0) :
    2*alpha*upperBranch alpha lam radiation-1=
      Real.sqrt (discriminant alpha lam radiation) := by
  unfold upperBranch
  field_simp
  <;> ring

theorem lower_branch_factor (alpha lam radiation : ℝ) (ha : alpha ≠ 0) :
    2*alpha*lowerBranch alpha lam radiation-1=
      -Real.sqrt (discriminant alpha lam radiation) := by
  unfold lowerBranch
  field_simp
  <;> ring

theorem upper_branch_is_root (alpha lam radiation : ℝ) (ha : alpha ≠ 0)
    (hd : 0 ≤ discriminant alpha lam radiation) :
    friedmannPolynomial alpha lam radiation (upperBranch alpha lam radiation)=0 := by
  have hs := Real.sq_sqrt hd
  have hf := upper_branch_factor alpha lam radiation ha
  have hfsq := congrArg (fun x : ℝ => x^2) hf
  dsimp at hfsq
  rw [hs] at hfsq
  unfold discriminant at hfsq
  unfold friedmannPolynomial
  have hmul : alpha*(alpha*(upperBranch alpha lam radiation)^2-
      upperBranch alpha lam radiation+lam+radiation)=0 := by nlinarith [hfsq]
  exact (mul_eq_zero.mp hmul).resolve_left ha

theorem lower_branch_is_root (alpha lam radiation : ℝ) (ha : alpha ≠ 0)
    (hd : 0 ≤ discriminant alpha lam radiation) :
    friedmannPolynomial alpha lam radiation (lowerBranch alpha lam radiation)=0 := by
  have hs := Real.sq_sqrt hd
  have hf := lower_branch_factor alpha lam radiation ha
  have hfsq := congrArg (fun x : ℝ => x^2) hf
  dsimp at hfsq
  rw [neg_sq, hs] at hfsq
  unfold discriminant at hfsq
  unfold friedmannPolynomial
  have hmul : alpha*(alpha*(lowerBranch alpha lam radiation)^2-
      lowerBranch alpha lam radiation+lam+radiation)=0 := by nlinarith [hfsq]
  exact (mul_eq_zero.mp hmul).resolve_left ha

theorem real_root_classification (alpha lam radiation y : ℝ) (ha : alpha ≠ 0)
    (hroot : friedmannPolynomial alpha lam radiation y=0) :
    y=lowerBranch alpha lam radiation ∨ y=upperBranch alpha lam radiation := by
  have hd := real_root_discriminant_nonnegative alpha lam radiation y hroot
  have hs := Real.sq_sqrt hd
  have hsq := root_discriminant_square alpha lam radiation y hroot
  have hprod : (2*alpha*y-1+Real.sqrt (discriminant alpha lam radiation))*
      (2*alpha*y-1-Real.sqrt (discriminant alpha lam radiation))=0 := by nlinarith
  have hden : 2*alpha ≠ 0 := mul_ne_zero (by norm_num) ha
  rcases mul_eq_zero.mp hprod with hh | hh
  · left
    unfold lowerBranch
    apply (eq_div_iff hden).2
    linarith
  · right
    unfold upperBranch
    apply (eq_div_iff hden).2
    linarith

theorem branch_bounds (alpha lam radiation : ℝ) (ha : 0 < alpha)
    (hsum : 0 ≤ lam+radiation) (hd : 0 ≤ discriminant alpha lam radiation) :
    0 ≤ lowerBranch alpha lam radiation ∧
    lowerBranch alpha lam radiation ≤ 1/(2*alpha) ∧
    1/(2*alpha) ≤ upperBranch alpha lam radiation ∧
    upperBranch alpha lam radiation ≤ 1/alpha := by
  have hs0 := Real.sqrt_nonneg (discriminant alpha lam radiation)
  have hsq := Real.sq_sqrt hd
  have hs1 : Real.sqrt (discriminant alpha lam radiation) ≤ 1 := by
    have hp := mul_nonneg (le_of_lt ha) hsum
    have hd1 : discriminant alpha lam radiation ≤ 1 := by unfold discriminant; linarith
    nlinarith [hd1]
  have hden : 0 < 2*alpha := by positivity
  unfold lowerBranch upperBranch
  constructor
  · exact div_nonneg (by linarith) (le_of_lt hden)
  constructor
  · exact (div_le_div_iff_of_pos_right hden).2 (by linarith)
  constructor
  · exact (div_le_div_iff_of_pos_right hden).2 (by linarith)
  · apply (div_le_iff₀ hden).2
    have he : 1/alpha*(2*alpha)=2 := by field_simp <;> ring
    rw [he]
    linarith

theorem lower_branch_positive (alpha lam radiation : ℝ) (ha : 0 < alpha)
    (hsum : 0 < lam+radiation) (hd : 0 ≤ discriminant alpha lam radiation) :
    0 < lowerBranch alpha lam radiation := by
  have hs0 := Real.sqrt_nonneg (discriminant alpha lam radiation)
  have hsq := Real.sq_sqrt hd
  have hs1 : Real.sqrt (discriminant alpha lam radiation) < 1 := by
    have hp := mul_pos ha hsum
    have hd1 : discriminant alpha lam radiation < 1 := by unfold discriminant; linarith
    nlinarith [hd1]
  unfold lowerBranch
  exact div_pos (by linarith) (by positivity)

theorem no_hubble_zero (alpha lam radiation H : ℝ)
    (hl : 0 ≤ lam) (hr : 0 < radiation)
    (hroot : friedmannPolynomial alpha lam radiation (H^2)=0) : H ≠ 0 := by
  intro hh
  subst H
  unfold friedmannPolynomial at hroot
  nlinarith

theorem positive_radiation_no_bounce (c A C a H Lambda : ℝ)
    (hc : 0 < c) (hC : 0 < C) (ha : 0 < a) (hl : 0 ≤ Lambda)
    (hconstraint : H^2=c*(A*H^4+C/a^4)+Lambda/3) : H ≠ 0 := by
  apply no_hubble_zero (c*A) (Lambda/3) (c*C/a^4) H (by positivity) (by positivity)
  exact (constraint_polynomial_equivalence c A C a H Lambda).mp hconstraint

theorem vacuum_de_sitter_root (alpha : ℝ) (ha : alpha ≠ 0) :
    friedmannPolynomial alpha 0 0 (1/alpha)=0 := by
  unfold friedmannPolynomial
  field_simp
  <;> ring

theorem vacuum_roots (alpha y : ℝ) (ha : alpha ≠ 0)
    (hroot : friedmannPolynomial alpha 0 0 y=0) : y=0 ∨ y=1/alpha := by
  have hp : y*(alpha*y-1)=0 := by
    unfold friedmannPolynomial at hroot
    nlinarith [hroot]
  rcases mul_eq_zero.mp hp with hh | hh
  · exact Or.inl hh
  · right
    apply (eq_div_iff ha).2
    linarith

theorem radiation_hasDerivAt (a : ℝ → ℝ) (t H c C : ℝ)
    (ha : HasDerivAt a (H*a t) t) (ha0 : a t ≠ 0) :
    HasDerivAt (fun u => c*C/(a u)^4) (-4*H*(c*C/(a t)^4)) t := by
  have hn : (a t)^4 ≠ 0 := pow_ne_zero 4 ha0
  convert (hasDerivAt_const t (c*C)).div (ha.pow 4) hn using 1 <;>
    dsimp <;> field_simp <;> ring

theorem polynomial_hasDerivAt (y radiation : ℝ → ℝ) (t alpha lam dy H : ℝ)
    (hy : HasDerivAt y dy t)
    (hr : HasDerivAt radiation (-4*H*radiation t) t) :
    HasDerivAt (fun u => friedmannPolynomial alpha lam (radiation u) (y u))
      ((2*alpha*y t-1)*dy-4*H*radiation t) t := by
  convert ((((hy.pow 2).const_mul alpha).sub hy).add_const lam).add hr using 1 <;>
    dsimp [friedmannPolynomial] <;> ring

theorem differentiated_constraint (y radiation : ℝ → ℝ) (t alpha lam dy H : ℝ)
    (hy : HasDerivAt y dy t)
    (hr : HasDerivAt radiation (-4*H*radiation t) t)
    (hzero : HasDerivAt (fun u => friedmannPolynomial alpha lam (radiation u) (y u)) 0 t) :
    (2*alpha*y t-1)*dy=4*H*radiation t := by
  have hh := (polynomial_hasDerivAt y radiation t alpha lam dy H hy hr).unique hzero
  linarith

theorem differentiated_local_roots (y radiation : ℝ → ℝ) (t alpha lam dy H : ℝ)
    (hy : HasDerivAt y dy t)
    (hr : HasDerivAt radiation (-4*H*radiation t) t)
    (hroot : Filter.Eventually
      (fun u => friedmannPolynomial alpha lam (radiation u) (y u)=0) (nhds t)) :
    (2*alpha*y t-1)*dy=4*H*radiation t := by
  exact differentiated_constraint y radiation t alpha lam dy H hy hr
    ((hasDerivAt_const t 0).congr_of_eventuallyEq hroot)

theorem branch_rate_formula (alpha y H radiation dy : ℝ)
    (hf : (2*alpha*y-1)*dy=4*H*radiation) (hne : 2*alpha*y-1 ≠ 0) :
    dy=4*H*radiation/(2*alpha*y-1) := by
  apply (eq_div_iff hne).2
  nlinarith [hf]

theorem merger_excludes_finite_rate (alpha y H radiation dy : ℝ)
    (hH : H ≠ 0) (hr : radiation ≠ 0) (hmerge : 2*alpha*y-1=0) :
    (2*alpha*y-1)*dy ≠ 4*H*radiation := by
  rw [hmerge, zero_mul]
  exact Ne.symm (mul_ne_zero (mul_ne_zero (by norm_num) hH) hr)

theorem zero_discriminant_merger (alpha lam radiation y : ℝ)
    (hroot : friedmannPolynomial alpha lam radiation y=0)
    (hd : discriminant alpha lam radiation=0) : 2*alpha*y-1=0 := by
  have hs := root_discriminant_square alpha lam radiation y hroot
  rw [hd] at hs
  nlinarith [sq_nonneg (2*alpha*y-1)]

theorem merger_not_differentiable (y radiation : ℝ → ℝ) (t alpha lam H : ℝ)
    (hH : H ≠ 0) (hr0 : radiation t ≠ 0)
    (hr : HasDerivAt radiation (-4*H*radiation t) t)
    (hroot : Filter.Eventually
      (fun u => friedmannPolynomial alpha lam (radiation u) (y u)=0) (nhds t))
    (hd : discriminant alpha lam (radiation t)=0) :
    ¬ DifferentiableAt ℝ y t := by
  intro hy
  have hrt : friedmannPolynomial alpha lam (radiation t) (y t)=0 := hroot.self_of_nhds
  have hm := zero_discriminant_merger alpha lam (radiation t) (y t) hrt hd
  have hf := differentiated_local_roots y radiation t alpha lam (deriv y t) H
    hy.hasDerivAt hr hroot
  exact merger_excludes_finite_rate alpha (y t) H (radiation t) (deriv y t) hH hr0 hm hf

end

end RenormalizedBranches
