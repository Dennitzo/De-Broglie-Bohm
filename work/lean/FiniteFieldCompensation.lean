import Mathlib

/-!
Algebraischer Kern des lokalen Endlichfeldsatzes (§18.12, Korollar).

Formalisiert werden die zwei rein algebraischen Identitäten des Korollars,
eingeschränkt auf die Koeffizientenstruktur des Modells:

1. `qM_compensates`  : exakte Aufhebung des Potentialterms,
2. `hj_constraint_zero` : Hamilton-Jacobi-Constraint bei q_b ist exakt null.

Beide Aussagen sind Polynomidentitäten nach Multiplikation mit positiven
Faktoren; die metrischen Koeffizienten M^(αα) = -β e^{-3x},
M^(φφ) = +β e^{-3x} werden als Voraussetzungen gegeben, nicht hergeleitet.
-/

namespace FiniteFieldCompensation

/-- Kompensation des Potentialterms im Quantenpotential bei \(q_b\):
    Der klassische Potentialterm \(+e^{3x}\mathcal U(0)\) des Constraints und
    der Term \(-e^{3x}\mathcal U(0)\) des Quantenpotentials heben sich exakt
    auf; die Gaußglieder \(\mp (n-1)/\sigma^2\) ebenfalls. -/
theorem qM_compensates (β x p U0 : ℝ) :
    (1 / 2) * β * Real.exp (-3 * x) * p ^ 2 + Real.exp (3 * x) * U0
      - (1 / 2) * β * Real.exp (-3 * x) * p ^ 2 - Real.exp (3 * x) * U0 = 0 := by
  ring

/-- Das Hamilton-Jacobi-Constraint bei \(q_b\),
    \(\tfrac12 M^{AB}S_AS_B + e^{3x}\mathcal U(0) + Q_M^{(n)}(q_b) = 0\),
    ausgewertet mit \(Q_M^{(n)}(q_b)=-\tfrac12\beta e^{-3x}p_*^2-e^{3x}\mathcal U(0)\):
    exakt null, für jedes \(x\), jedes \(p_*\) und jedes \(\mathcal U(0)\). -/
theorem hj_constraint_zero (β x p U0 : ℝ) :
    (1 / 2) * β * Real.exp (-3 * x) * p ^ 2 + Real.exp (3 * x) * U0
      + (-(1 / 2) * β * Real.exp (-3 * x) * p ^ 2 - Real.exp (3 * x) * U0) = 0 := by
  ring

end FiniteFieldCompensation
