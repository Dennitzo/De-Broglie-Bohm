import Mathlib

/-!
# ADM-Spuridentitäten und kinetische Legendre-Transformation in drei Dimensionen

Formalisierung von Abschnitt 29.6.1 der Publikation: Beziehung zwischen extrinsischer
Krümmung und kanonischem Impuls, Umkehrung dieser Beziehung, Kontraktion der
DeWitt-Supermetrik, Polarisation der kinetischen Quadratform (das ist die Hamiltonsche
Ableitung `∂H_kin/∂π^{ij}` in Richtung eines Testtensors) sowie die **korrigierte**
Identität `π^{ij} K_{ij} = L_kin / N` statt der früheren, um den Faktor `N` falschen
Aussage `π^{ij} K_{ij} = 2 L_kin`.

## Was formalisiert ist und was nicht

Alle Sätze gelten in einem **orthonormalen Rahmen** der Raummetrik (`h_{ij} = δ_{ij}`);
dort ist Indeshochziehen die Identität, und die Tensor­kontraktion `A_{ij} B_{ij}` fällt
mit `Matrix.trace (A * B)` zusammen. Die allgemein kovariante Fassung steht im Manuskript.

Der Verschiebungsanteil `£_{\vec N} h` wird **nicht** formalisiert: die Sätze über
`extCurvature`, `lagrangianKinetic` und `euler_shiftFree` erfassen den reinen
Lapse-Anteil `\dot h = 2 N K` (Verschiebung null). Die Behauptung
`π^{ij}(\dot h - £_{\vec N}h)_{ij} = 2 L_kin` des Manuskripts ist damit nur in diesem
Spezialfall bewiesen; das allgemeine Lie-Ableitungs­glied bleibt mathematischer Text.

Es werden ausschließlich endliche Matrixidentitäten über `Matrix (Fin 3) (Fin 3) ℝ`,
Algebra ohne Analysis bis auf die Polarisations­identität, und keine numerischen
Berechnungen verwendet; keine `sorry`.
-/

namespace ADMConstraintLegendre
open Matrix
noncomputable section

/-! ### Kontraktion als Matrixspur -/

/-- Die Kontraktion `A_{ij} B_{ij}`. Im orthonormalen Rahmen stimmt sie für symmetrische
Tensoren mit der Spur des Matrixprodukts überein; wir verwenden diese Spur als Definition,
dadurch sind alle folgenden Identitäten reine Algebra in den Einträgen. -/
def con (A B : Matrix (Fin 3) (Fin 3) ℝ) : ℝ := Matrix.trace (A * B)

theorem con_add_left (A B C : Matrix (Fin 3) (Fin 3) ℝ) :
    con A (B + C) = con A B + con A C := by
  unfold con
  rw [Matrix.mul_add, Matrix.trace_add]

theorem con_add_right (A B C : Matrix (Fin 3) (Fin 3) ℝ) :
    con (A + B) C = con A C + con B C := by
  unfold con
  rw [Matrix.add_mul, Matrix.trace_add]

theorem con_sub_left (A B C : Matrix (Fin 3) (Fin 3) ℝ) :
    con A (B - C) = con A B - con A C := by
  unfold con
  rw [Matrix.mul_sub, Matrix.trace_sub]

theorem con_sub_right (A B C : Matrix (Fin 3) (Fin 3) ℝ) :
    con (A - B) C = con A C - con B C := by
  unfold con
  rw [Matrix.sub_mul, Matrix.trace_sub]

/-- Die Kontraktion ist symmetrisch (Zyklizität der Spur). -/
theorem con_comm (A B : Matrix (Fin 3) (Fin 3) ℝ) : con A B = con B A := by
  unfold con
  rw [Matrix.trace_mul_comm]

theorem con_smul_left (r : ℝ) (A B : Matrix (Fin 3) (Fin 3) ℝ) :
    con (r • A) B = r * con A B := by
  unfold con
  rw [smul_mul_assoc, Matrix.trace_smul, smul_eq_mul]

theorem con_smul_right (r : ℝ) (A B : Matrix (Fin 3) (Fin 3) ℝ) :
    con A (r • B) = r * con A B := by
  unfold con
  rw [Matrix.trace_mul_comm, smul_mul_assoc, Matrix.trace_smul, smul_eq_mul,
    Matrix.trace_mul_comm]

/-- `con A 1` ist die Spur. -/
theorem con_one_right (A : Matrix (Fin 3) (Fin 3) ℝ) : con A 1 = Matrix.trace A := by
  unfold con
  rw [Matrix.mul_one]

theorem con_one_left (A : Matrix (Fin 3) (Fin 3) ℝ) : con 1 A = Matrix.trace A := by
  unfold con
  rw [Matrix.one_mul]

/-- Die Raumdimension als Spur der Einheitsmatrix; getrennt bewiesen, weil `trace_one`
über `Fintype.card (Fin 3)` geht und diese Zahl sonst in algebraischen Zielen stehen bleibt. -/
theorem trace_one_three : Matrix.trace (1 : Matrix (Fin 3) (Fin 3) ℝ) = 3 := by
  norm_num

/-- Die Selbstkontraktion der Einheitsmatrix ist die Raumdimension `3`. -/
theorem con_one_one : con (1 : Matrix (Fin 3) (Fin 3) ℝ) 1 = 3 := by
  rw [con_one_right, trace_one_three]

/-! ### Spurumkehr, kanonischer Impuls und seine Umkehrung -/

/-- `T_c(A) = A - c · tr(A) · 1`. Die ADM-Spurumkehr ist der Fall `c = 1`
(`A^{ij} - h^{ij} A`), die Hamiltonsche Ableitung der kinetischen Quadratform führt
dagegen auf `c = 1/2`. -/
def traceReverse (c : ℝ) (A : Matrix (Fin 3) (Fin 3) ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  A - c • Matrix.trace A • (1 : Matrix (Fin 3) (Fin 3) ℝ)

theorem trace_traceReverse (c : ℝ) (A : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix.trace (traceReverse c A) = (1 - 3 * c) * Matrix.trace A := by
  unfold traceReverse
  rw [Matrix.trace_sub, Matrix.trace_smul, Matrix.trace_smul, trace_one_three]
  ring

/-- Die Kontraktion einer spurumgekehrten Matrix mit einem Testtensor. -/
theorem con_traceReverse_left (c : ℝ) (A B : Matrix (Fin 3) (Fin 3) ℝ) :
    con (traceReverse c A) B = con A B - c * Matrix.trace A * Matrix.trace B := by
  unfold traceReverse
  rw [con_sub_right, con_smul_left, con_smul_left, con_one_left]
  ring

theorem con_traceReverse_right (c : ℝ) (A B : Matrix (Fin 3) (Fin 3) ℝ) :
    con A (traceReverse c B) = con A B - c * Matrix.trace A * Matrix.trace B := by
  unfold traceReverse
  rw [con_sub_left, con_smul_right, con_smul_right, con_one_right]
  ring

/-- Selbstkontraktion der ADM-Spurumkehr (`c = 1`). -/
theorem con_traceReverse_self (K : Matrix (Fin 3) (Fin 3) ℝ) :
    con (traceReverse 1 K) (traceReverse 1 K) = con K K + (Matrix.trace K) ^ 2 := by
  rw [con_traceReverse_left, con_traceReverse_right, trace_traceReverse]
  ring

/-- Kontraktion der ADM-Spurumkehr mit `K` selbst; das ist `π_{ij}K^{ij}/c`. -/
theorem con_traceReverse_K (K : Matrix (Fin 3) (Fin 3) ℝ) :
    con (traceReverse 1 K) K = con K K - (Matrix.trace K) ^ 2 := by
  rw [con_traceReverse_left]
  ring

/-- Der kanonische Impuls `π^{ij} = c (K^{ij} - h^{ij} K)` mit `c = √h/(2κ)`; das ist die
Spurumkehr mit Koeffizient `1`. -/
def piOf (c : ℝ) (K : Matrix (Fin 3) (Fin 3) ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  c • traceReverse 1 K

/-- Die Impulsspur `π = h_{ij}π^{ij} = -2c · tr K`, mit `c = √h/(2κ)` also
`π = -√h K/κ` (Manuskript, Rechnung nach Gl. (29.6.1d)). Bemerkenswert ist die Asymmetrie:
vorwärts steht der Spurumkehr-Koeffizient `1`, rückwärts `½`; genau das ist die
DeWitt-Kontraktion unten. -/
theorem trace_piOf (c : ℝ) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix.trace (piOf c K) = -2 * c * Matrix.trace K := by
  unfold piOf
  rw [Matrix.trace_smul, trace_traceReverse]
  ring

/-- **Umkehrung der Spurumkehr** mit Impulsspur (Manuskript, Rechnung vor Gl. (29.6.1f)):
`K_{ij} = c⁻¹ (π_{ij} - ½ h_{ij} π)`. -/
theorem piOf_roundtrip (c : ℝ) (hc : c ≠ 0) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    c⁻¹ • (piOf c K - (Matrix.trace (piOf c K) / 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ)) = K := by
  rw [trace_piOf c K]
  unfold piOf traceReverse
  ext i j
  simp only [Matrix.smul_apply, Matrix.sub_apply, Matrix.one_apply, smul_eq_mul]
  field_simp [hc]
  ring

/-! ### Die kinetischen Quadratformen -/

/-- `K_{ij}K^{ij} - K^{2}` (im Manuskript die Größe in Gl. (29.6.1b)). -/
def kineticForm (K : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  con K K - (Matrix.trace K) ^ 2

/-- `π^{ij}π_{ij} - ½π²` (die Größe in Gl. (29.6.1f), ohne den Vorfaktor). -/
def quadForm (P : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  con P P - (Matrix.trace P) ^ 2 / 2

/-- **Kernidentität** (Manuskript, Rechnung vor Gl. (29.6.1f)):
`π^{ij}π_{ij} - ½π² = c² (K_{ij}K^{ij} - K²)`. Mit `c = √h/(2κ)` ist das genau
`(h/4κ²)(K_{ij}K^{ij} - K²)`, womit die Impulsdarstellung des kinetischen Hamiltonians
Gl. (29.6.1k) begründet ist. -/
theorem quadForm_piOf (c : ℝ) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    quadForm (piOf c K) = c ^ 2 * kineticForm K := by
  unfold quadForm kineticForm
  rw [trace_piOf c K]
  unfold piOf
  rw [con_smul_left, con_smul_right, con_traceReverse_self]
  ring

/-- Polarisation von `quadForm` mit ganzzahligen Koeffizienten (`2·quadForm`, damit `ring`
ohne `field_simp` auskommt). -/
def quadFormTwice (P : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  2 * con P P - (Matrix.trace P) ^ 2

/-- **Polarisationsidentität.** Lineiner Term in Richtung `S`:
`2Q(π + tS) = 2Q(π) + 4t · (π_{ij} - ½h_{ij}π) S^{ij} + t² · 2Q(S)`.
Da der lineine Term die Richtungsableitung ist, folgt daraus für den kinetischen
Hamiltonian `H_kin = (N/(2c))·2Q(π)` von Gl. (29.6.1k):
`∂H_kin/∂π^{ij} = (2/c)(π_{ij} - ½h_{ij}π) = (4κ/√h)(π_{ij} - ½h_{ij}π) = 2K_{ij}`,
also Gl. (29.6.1m) des Manuskripts; zusammen mit dem DeWitt-Satz ist diese Ableitung
genau die Kontraktion der Supermetrik. -/
theorem quadFormTwice_polar (P S : Matrix (Fin 3) (Fin 3) ℝ) (t : ℝ) :
    quadFormTwice (P + t • S) =
      quadFormTwice P + t * (4 * con (traceReverse (1 / 2) P) S)
        + t ^ 2 * quadFormTwice S := by
  have hc : con (P + t • S) (P + t • S)
      = con P P + t * (con P S + con S P) + t * t * con S S := by
    simp only [con_add_right, con_add_left, con_smul_left, con_smul_right]
    ring
  have hcs : con S P = con P S := con_comm S P
  have ht : Matrix.trace (P + t • S) = Matrix.trace P + t * Matrix.trace S := by
    simp only [Matrix.trace_add, Matrix.trace_smul, smul_eq_mul]
  unfold quadFormTwice
  rw [hc, hcs, ht, con_traceReverse_left]
  ring

/-! ### DeWitt-Supermetrik und ihre Kontraktion -/

/-- Die DeWitt-Supermetrik im orthonormalen Rahmen (Dichtegewicht `1`):
`G_{ijkl} = ½(δ_{ik}δ_{jl} + δ_{il}δ_{jk} - δ_{ij}δ_{kl})`. -/
def dewitt : Fin 3 → Fin 3 → Fin 3 → Fin 3 → ℝ := fun i j k l =>
  (1 / 2 : ℝ) * ((if i = k then if j = l then 1 else 0 else 0)
    + (if i = l then if j = k then 1 else 0 else 0)
    - (if i = j then if k = l then 1 else 0 else 0))

/-- Kontraktion `G_{ijkl} π^{kl}` über beide Indizes. -/
def dewittContract (P : Matrix (Fin 3) (Fin 3) ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => ∑ k : Fin 3, ∑ l : Fin 3, dewitt i j k l * P k l

/-- **Auswertung der Delta-Kontraktion** Eintrag für Eintrag. Ohne jede Symmetrieannahme
ergibt die Kontraktion von `G_{ijkl} = ½(δ_{ik}δ_{jl}+δ_{il}δ_{jk}-δ_{ij}δ_{kl})` über
`kl` gerade `½(π_{ij} + π_{ji}) - ½δ_{ij}π`. -/
theorem dewittContract_apply (P : Matrix (Fin 3) (Fin 3) ℝ) (i j : Fin 3) :
    dewittContract P i j = (1 / 2 : ℝ) * P i j + (1 / 2 : ℝ) * P j i
      - (if i = j then (1 / 2 : ℝ) * Matrix.trace P else 0) := by
  classical
  unfold dewittContract dewitt
  fin_cases i <;> fin_cases j <;>
    simp only [Fin.sum_univ_three, Matrix.trace, mul_add, mul_sub, sub_mul, ite_mul] <;>
    (try decide) <;>
    simp <;>
    ring

/-- **Kontraktion der DeWitt-Supermetrik** (Gl. (29.6.1o)):
`G_{ijkl}π^{kl} = π_{ij} - ½ h_{ij} π`, im orthonormalen Rahmen also
`G·π = traceReverse (1/2) π`.

Die Symmetrieannahme ist **wesentlich**, kein Formalisierungstrick: nach
`dewittContract_apply` liefert die Kontraktion für `i ≠ j` den Wert
`½(π_{ij} + π_{ji})`, und das ist genau dann `π_{ij}`, wenn `π_{ij} = π_{ji}`. Der
kanonische Impuls ist symmetrisch, weil `K_{ij}` es ist; für einen nicht symmetrischen
Tensor ist die Behauptung falsch. -/
theorem dewittContract_eq_traceReverse (P : Matrix (Fin 3) (Fin 3) ℝ) (hP : P.IsSymm) :
    dewittContract P = traceReverse (1 / 2) P := by
  classical
  ext i j
  rw [dewittContract_apply, traceReverse]
  by_cases h : i = j
  · subst h
    simp only [if_true, Matrix.smul_apply, Matrix.sub_apply, Matrix.one_apply, smul_eq_mul]
    ring
  · rw [show P j i = P i j from hP.apply i j]
    simp only [if_neg h, Matrix.smul_apply, Matrix.sub_apply, Matrix.one_apply, smul_eq_mul]
    ring

/-! ### Extrinsische Krümmung, kinetische Lagrange-Dichte und die korrigierte Identität -/

/-- `K_{ij} = (\dot h_{ij} - £_{\vec N}h_{ij})/(2N)` im Spezialfall ohne Verschiebung:
`K = \dot h/(2N)`. -/
def extCurvature (N : ℝ) (hdot : Matrix (Fin 3) (Fin 3) ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  (1 / (2 * N)) • hdot

/-- Die kinetische Lagrange-Dichte `L_kin = N c (K_{ij}K^{ij} - K²)` von Gl. (29.6.1b)
mit `c = √h/(2κ)`; hier ohne `\sqrt h`-Gewicht, da orthonormaler Rahmen (`\sqrt h = 1`). -/
def lagrangianKinetic (N c : ℝ) (hdot : Matrix (Fin 3) (Fin 3) ℝ) : ℝ :=
  N * c * kineticForm (extCurvature N hdot)

/-- `\dot h = 2NK` ist konsistent mit der Definition von `K`: einsetzbar in allen
folgenden Sätzen. -/
theorem extCurvature_two_N (N : ℝ) (hN : N ≠ 0) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    extCurvature N ((2 * N) • K) = K := by
  unfold extCurvature
  rw [smul_smul]
  field_simp [hN]
  rw [one_smul]

/-- **Euler-Identität ohne Verschiebung** (Manuskript, Gl. (29.6.1e)):
`π^{ij}(\dot h)_{ij} = 2 L_kin` für `\dot h = 2 N K`. Ohne den Faktor `N` in der
Lapse-Kontraktion wäre das falsch; das ist gerade die Stelle des früheren Fehlers. -/
theorem euler_shiftFree (N c : ℝ) (hN : N ≠ 0) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    con (piOf c K) ((2 * N) • K) = 2 * lagrangianKinetic N c ((2 * N) • K) := by
  unfold lagrangianKinetic kineticForm
  rw [extCurvature_two_N N hN K]
  unfold piOf
  rw [con_smul_left, con_smul_right, con_traceReverse_K]
  ring

/-- **Die korrigierte Identität** (Manuskript Gl. (29.6.1d)):
`π^{ij}K_{ij} = L_kin/N`, also **nicht** `2L_kin`. -/
theorem piK_eq_lagrangian_div_N (N c : ℝ) (hN : N ≠ 0) (K : Matrix (Fin 3) (Fin 3) ℝ) :
    con (piOf c K) K = lagrangianKinetic N c ((2 * N) • K) / N := by
  unfold lagrangianKinetic kineticForm
  rw [extCurvature_two_N N hN K]
  unfold piOf
  rw [con_smul_left, con_traceReverse_K]
  field_simp [hN]

/-- **Die frühere Behauptung ist falsch.** Für `N = 2`, `c = 1` und `K = δ_{ij}`
(konsistent mit `\dot h = 2N K = 4 δ_{ij}`) gilt `π^{ij}K_{ij} = -6`, aber
`2 L_kin = -24`. Das Gegenbeispiel sichert, dass die Korrektur von Gl. (29.6.1d) keine
Konvention, sondern eine echte Fehlerbehebung ist. -/
theorem old_claim_is_false :
    ¬ (con (piOf 1 (extCurvature 2 (((2 : ℝ) * 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ))))
        (extCurvature 2 (((2 : ℝ) * 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ)))) = 2 *
      lagrangianKinetic 2 1 (((2 : ℝ) * 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ)) := by
  intro h
  have hN₂ : (2 : ℝ) ≠ 0 := by norm_num
  have h₁ : extCurvature 2 (((2 : ℝ) * 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ))
      = (1 : Matrix (Fin 3) (Fin 3) ℝ) := extCurvature_two_N 2 hN₂ _
  rw [h₁] at h
  have h₂ : con (piOf 1 (1 : Matrix (Fin 3) (Fin 3) ℝ)) (1 : Matrix (Fin 3) (Fin 3) ℝ) = -6 := by
    unfold piOf
    rw [con_smul_left, con_traceReverse_left, trace_one_three, con_one_one]
    norm_num
  have h₃ : lagrangianKinetic 2 1 (((2 : ℝ) * 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ)) = -12 := by
    unfold lagrangianKinetic
    rw [h₁]
    have h₄ : kineticForm (1 : Matrix (Fin 3) (Fin 3) ℝ) = -6 := by
      unfold kineticForm
      rw [con_one_one, trace_one_three]
      norm_num
    rw [h₄]
    norm_num
  rw [h₂, h₃] at h
  norm_num at h

end
end ADMConstraintLegendre

#print axioms ADMConstraintLegendre.piOf_roundtrip
#print axioms ADMConstraintLegendre.quadForm_piOf
#print axioms ADMConstraintLegendre.dewittContract_eq_traceReverse
#print axioms ADMConstraintLegendre.piK_eq_lagrangian_div_N
#print axioms ADMConstraintLegendre.old_claim_is_false
