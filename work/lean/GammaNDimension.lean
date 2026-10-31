import Mathlib

open scoped Real

/--
For every finite field number `n`, the logarithmic derivative of the coefficient
`-β^((1-n)/2) exp(3(n-1)α/2)` is `3(n-1)/2`.
For `β > 0` this coefficient equals `μ_n M^{αα}` by the elementary calculation in the manuscript.
-/
theorem gamma_coeff_deriv
    (β : ℝ) (hβ : 0 < β) (n : ℕ) (α : ℝ) :
    HasDerivAt
      (fun x : ℝ =>
        -(β ^ ((1 - (n : ℝ)) / 2)) *
        Real.exp ((3 * ((n : ℝ) - 1)) / 2 * x))
      ((3 * ((n : ℝ) - 1)) / 2 *
        (-(β ^ ((1 - (n : ℝ)) / 2)) *
        Real.exp ((3 * ((n : ℝ) - 1)) / 2 * α)))
      α := by
  have hexp :
      HasDerivAt
        (fun x : ℝ => Real.exp ((3 * ((n : ℝ) - 1)) / 2 * x))
        (Real.exp ((3 * ((n : ℝ) - 1)) / 2 * α) *
          ((3 * ((n : ℝ) - 1)) / 2))
        α := by
    simpa using
      (Real.hasDerivAt_exp ((3 * ((n : ℝ) - 1)) / 2 * α)).comp α
        ((hasDerivAt_id α).const_mul ((3 * ((n : ℝ) - 1)) / 2))
  convert hexp.const_mul (-(β ^ ((1 - (n : ℝ)) / 2))) using 1
  ring

#print axioms gamma_coeff_deriv
