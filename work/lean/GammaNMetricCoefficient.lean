import Mathlib

open scoped Real

noncomputable section

/--
`muMetricCoeffExtra β α m` ist $\mu_{m+1}M^{\alpha\alpha}$, wobei $m$ die Anzahl
der zusätzlichen Skalarfeldkoordinaten neben der Uhrkoordinate bezeichnet.
Für $m=n-1$ ist das genau der Faktor des $\alpha$-Teils im
Laplace-Beltrami-Operator bei $n$ Skalarfeldern.
-/
def muMetricCoeffExtra (β α : ℝ) (m : ℕ) : ℝ :=
  (β⁻¹ * Real.exp (3 * α)) ^ (((m : ℝ) + 2) / 2) *
    (-β * Real.exp (-3 * α))

/--
Die geschlossene Form von `muMetricCoeffExtra`:
\[
-\beta^{-m/2}\exp(3m\alpha/2).
\]
-/
def simplifiedCoeffExtra (β α : ℝ) (m : ℕ) : ℝ :=
  -(β ^ (-((m : ℝ) / 2)) *
    Real.exp (((3 : ℝ) / 2) * (m : ℝ) * α))

/--
Für jedes endliche $m$ folgt aus $\mu_{m+1}=(\beta^{-1}e^{3\alpha})^{(m+2)/2}$
und $M^{\alpha\alpha}=-\beta e^{-3\alpha}$ die geschlossene Form.
-/
theorem muMetricCoeffExtra_eq_simplified
    (hβ : 0 < β) (α : ℝ) (m : ℕ) :
    muMetricCoeffExtra β α m = simplifiedCoeffExtra β α m := by
  set a : ℝ := ((m : ℝ) + 2) / 2 with ha
  set c : ℝ := ((3 : ℝ) / 2) * (m : ℝ) * α with hc
  have hbase : 0 < β⁻¹ * Real.exp (3 * α) := by positivity
  have hβexp : 0 < β * Real.exp (-3 * α) := by positivity
  have hrpow : 0 < (β⁻¹ * Real.exp (3 * α)) ^ a :=
    Real.rpow_pos_of_pos hbase a
  have hlhs_pos : 0 < β * Real.exp (-3 * α) *
      (β⁻¹ * Real.exp (3 * α)) ^ a :=
    mul_pos hβexp hrpow
  have hrhs_rpow : 0 < β ^ (-((m : ℝ) / 2)) :=
    Real.rpow_pos_of_pos hβ (-((m : ℝ) / 2))
  have hrhs_pos : 0 < β ^ (-((m : ℝ) / 2)) * Real.exp c :=
    mul_pos hrhs_rpow (Real.exp_pos c)

  have hlog_lhs :
      Real.log (β * Real.exp (-3 * α) *
        (β⁻¹ * Real.exp (3 * α)) ^ a) =
        Real.log β - 3 * α + a * (-Real.log β + 3 * α) := by
    rw [Real.log_mul (ne_of_gt hβexp) (ne_of_gt hrpow),
        Real.log_mul (ne_of_gt hβ) (ne_of_gt (Real.exp_pos _)),
        Real.log_rpow hbase,
        Real.log_exp,
        Real.log_mul (ne_of_gt (inv_pos.mpr hβ))
          (ne_of_gt (Real.exp_pos _)),
        Real.log_inv,
        Real.log_exp]
    ring

  have hlog_rhs :
      Real.log (β ^ (-((m : ℝ) / 2)) * Real.exp c) =
        -((m : ℝ) / 2) * Real.log β + c := by
    rw [Real.log_mul (ne_of_gt hrhs_rpow)
          (ne_of_gt (Real.exp_pos c)),
        Real.log_rpow hβ,
        Real.log_exp]

  have hpos_eq :
      β * Real.exp (-3 * α) *
        (β⁻¹ * Real.exp (3 * α)) ^ a =
        β ^ (-((m : ℝ) / 2)) * Real.exp c := by
    rw [← Real.exp_log hlhs_pos,
        ← Real.exp_log hrhs_pos,
        hlog_lhs, hlog_rhs, ha, hc]
    ring

  unfold muMetricCoeffExtra simplifiedCoeffExtra
  rw [← ha, ← hc]
  have hpull :
      (β⁻¹ * Real.exp (3 * α)) ^ a * (-β * Real.exp (-3 * α)) =
        -(β * Real.exp (-3 * α) *
          (β⁻¹ * Real.exp (3 * α)) ^ a) := by
    simp only [neg_mul, mul_neg, mul_comm, mul_left_comm]
  rw [hpull, hpos_eq]

/--
Für jedes endliche $m$ hat die geschlossene Form die logarithmische
$\alpha$-Ableitung $3m/2$.  Für $n$ Skalarfelder ist $m=n-1$, also
$\gamma_n=\frac32(n-1)$.
-/
theorem gamma_coeff_deriv_extra
    (β : ℝ) (_hβ : 0 < β) (m : ℕ) (α : ℝ) :
    HasDerivAt
      (fun x : ℝ => simplifiedCoeffExtra β x m)
      (((3 : ℝ) / 2) * (m : ℝ) * simplifiedCoeffExtra β α m)
      α := by
  unfold simplifiedCoeffExtra
  have hexp :
      HasDerivAt
        (fun x : ℝ => Real.exp (((3 : ℝ) / 2) * (m : ℝ) * x))
        (Real.exp (((3 : ℝ) / 2) * (m : ℝ) * α) *
          (((3 : ℝ) / 2) * (m : ℝ)))
        α := by
    simpa using
      (Real.hasDerivAt_exp (((3 : ℝ) / 2) * (m : ℝ) * α)).comp α
        ((hasDerivAt_id α).const_mul (((3 : ℝ) / 2) * (m : ℝ)))
  convert hexp.const_mul (-(β ^ (-((m : ℝ) / 2)))) using 1
  · funext x
    rw [neg_mul]
  · rw [mul_neg, neg_mul,
        mul_comm (((3 : ℝ) / 2) * (m : ℝ))
          (β ^ (-((m : ℝ) / 2)) *
            Real.exp (((3 : ℝ) / 2) * (m : ℝ) * α)),
        mul_assoc]

/--
Endliche-Indexmengen-Fassung: Ist $J$ die endliche Menge der zusätzlichen
Skalarfeldkoordinaten, so ist der Operator-Koeffizient
\[
\gamma(J)=\frac{3}{2}\,|J|.
\]
Bei $n$ Skalarfeldern ist $J$ die Menge der $n-1$ Koordinaten neben der Uhr,
also $|J|=n-1$.
-/
theorem gamma_coeff_deriv_for_fintype {J : Type} [Fintype J]
    (β : ℝ) (_hβ : 0 < β) (α : ℝ) :
    HasDerivAt
      (fun x : ℝ =>
        simplifiedCoeffExtra β x (Fintype.card J))
      (((3 : ℝ) / 2) * (Fintype.card J : ℝ) *
        simplifiedCoeffExtra β α (Fintype.card J))
      α := by
  simpa using gamma_coeff_deriv_extra β _hβ (Fintype.card J) α

/--
Übersetzung auf die im Manuskript benutzte Zählung: Bei $n\ge1$
Skalarfeldern ist $m=n-1$ und damit $\gamma_n=\frac32(n-1)$.
-/
theorem gamma_coeff_deriv_original
    (β : ℝ) (_hβ : 0 < β) (n : ℕ) (hn : 1 ≤ n) (α : ℝ) :
    HasDerivAt
      (fun x : ℝ => simplifiedCoeffExtra β x (n - 1))
      (((3 : ℝ) / 2) * ((n : ℝ) - 1) *
        simplifiedCoeffExtra β α (n - 1))
      α := by
  have hcast : ((n - 1 : ℕ) : ℝ) = (n : ℝ) - 1 := by
    simpa using Nat.cast_sub hn
  simpa [hcast] using gamma_coeff_deriv_extra β _hβ (n - 1) α

#print axioms muMetricCoeffExtra_eq_simplified
#print axioms gamma_coeff_deriv_extra
#print axioms gamma_coeff_deriv_for_fintype
#print axioms gamma_coeff_deriv_original
