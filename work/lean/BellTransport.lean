import Mathlib

/-!
# Conservation of a density ratio along a common differentiable flow

Here rho and q are densities evaluated along the same trajectory. Their
material derivatives obey the same compression law, with theta the local
velocity divergence. This proves a derivative identity under these explicit
hypotheses; it does not assert existence or regularity of a Bohmian flow.
-/

namespace BohmTransport

theorem continuity_ratio_numerator (rho q drho dq dxrho dxq v theta : ℝ)
    (hrho : drho + v * dxrho + rho * theta = 0)
    (hq : dq + v * dxq + q * theta = 0) :
    (drho * q - rho * dq) + v * (dxrho * q - rho * dxq) = 0 := by
  linear_combination q * hrho - rho * hq

theorem continuity_ratio_transport (rho q drho dq dxrho dxq v theta : ℝ)
    (hrho : drho + v * dxrho + rho * theta = 0)
    (hq : dq + v * dxq + q * theta = 0) :
    (drho * q - rho * dq) / q ^ 2 +
      v * ((dxrho * q - rho * dxq) / q ^ 2) = 0 := by
  calc
    (drho * q - rho * dq) / q ^ 2 +
        v * ((dxrho * q - rho * dxq) / q ^ 2)
        = ((drho * q - rho * dq) + v * (dxrho * q - rho * dxq)) / q ^ 2 := by
          ring
    _ = 0 := by
      rw [continuity_ratio_numerator rho q drho dq dxrho dxq v theta hrho hq]
      simp

theorem ratio_hasDerivAt_zero (rho q : ℝ → ℝ) (t theta : ℝ)
    (hrho : HasDerivAt rho (-rho t * theta) t)
    (hq : HasDerivAt q (-q t * theta) t)
    (hq_nonzero : q t ≠ 0) :
    HasDerivAt (fun u => rho u / q u) 0 t := by
  convert hrho.div hq hq_nonzero using 1
  ring

/- Algebraic product-rule expression for Pi = C * ((R')^2 - R * R'').
   R, Rp, Rpp and Rppp represent R and its first three spatial derivatives;
   this first identity does not itself establish those derivative relations. -/
theorem quantum_stress_product_rule_algebra (C R Rp Rpp Rppp : ℝ) :
    C * (2 * Rp * Rpp - (Rp * Rpp + R * Rppp)) =
      C * (Rp * Rpp - R * Rppp) := by
  ring

/- With rho = R^2 and Q = -C * R''/R, the quotient-rule expression for
   rho * Q' equals the product-rule expression for Pi'. Nodes are excluded. -/
theorem quantum_force_stress_balance (C R Rp Rpp Rppp : ℝ) (hR : R ≠ 0) :
    R ^ 2 * (-C * ((Rppp * R - Rpp * Rp) / R ^ 2)) =
      C * (Rp * Rpp - R * Rppp) := by
  field_simp [hR]
  ring

end BohmTransport

#print axioms BohmTransport.continuity_ratio_numerator
#print axioms BohmTransport.continuity_ratio_transport
#print axioms BohmTransport.ratio_hasDerivAt_zero
#print axioms BohmTransport.quantum_stress_product_rule_algebra
#print axioms BohmTransport.quantum_force_stress_balance
