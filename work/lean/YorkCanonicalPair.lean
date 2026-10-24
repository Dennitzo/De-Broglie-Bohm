import Mathlib

/-!
# Das kanonische York-Paar `{T, P_T} = 1` und die reduzierte York-Dynamik

Formalisierung der Abschnitte 28.9 und 29.15 der Publikation. Ausgehend vom kanonischen Paar
`(a, p_a)` wird die York-Zeit

```
T = p_a / (3a²) ,        P_T = −a³                                     (28.9)
```

betrachtet. **Bewiesen** wird, dass `(T, P_T)` erneut kanonisch ist, und zwar durch
tatsächliches Differenzieren beider Definitionen nach beiden Variablen:

```
∂T/∂a = −2p_a/(3a³) ,  ∂T/∂p_a = 1/(3a²) ,  ∂P_T/∂a = −3a² ,  ∂P_T/∂p_a = 0
{T, P_T} = (∂T/∂a)(∂P_T/∂p_a) − (∂T/∂p_a)(∂P_T/∂a) = 0 − (1/(3a²))·(−3a²) = +1 .
```

Der Vorfaktor `3` in `T = p_a/(3a²)` ist genau die Normierung, die den Wert `+1` liefert; er
ist also nicht frei wählbar, sondern wird durch die Kanonizität erzwungen.

Dazu kommen die Algebra des Hamilton-Zwangs
`C = −p_φ²/(2P_T) + 6πG T² P_T − P_T V(φ)` (Multiplikation mit `2P_T` gibt `P_T²A − p_φ²`
mit `A = 12πGT² − 2V`), die Verifikation, dass der Zweig `P_T = −p_φ/√A` den Zwang erfüllt,
und die beiden Hamilton-Gleichungen auf dem Zweig `σ = sgn(p_φ)` durch unabhängige
Differentiation:

```
∂H/∂p_φ = σ/√A ,      ∂H/∂φ = σ p_φ V'(φ)/A^{3/2} ,     p_φ' = −∂H/∂φ .
```

Das ist die Rechnung, die im Manuskript zur Korrektur des gedruckten Vorzeichens in
[RoserValentini2014] Gl. (3.11) führt. Für das freie masselose Feld wird schließlich
`d/dT [ln(−T)/√(12πG) + c] = 1/(√(12πG)·T)` für `T < 0` nachgerechnet und mit
`−1/√(12πGT²)` identifiziert.

## Abgrenzung

Formalisiert ist die **kanonische Struktur des Minisuperspace-Modells** und die Algebra der
reduzierten Dynamik. Nicht formalisiert: die Herleitung von `T = p_a/(3a²)` aus der ADM-Zerlegung,
die Existenz einer CMC-Foliation im allgemeinen Fall, die Identität `dT/dt = Φ̇²` (sie hängt an
den Friedmann-Gleichungen und Konventionsfragen zu `p_a`), und alles Quantenmechanische.
Kein numerischer Test, kein `sorry`.
-/

namespace YorkCanonicalPair

open Real

noncomputable section

/-! ### Die Transformation und ihre partiellen Ableitungen -/

/-- York-Zeit `T = p_a/(3a²)` als Funktion des kanonischen Paars `(a, p_a)`. -/
noncomputable def yorkTime (a pa : ℝ) : ℝ := pa / (3 * a ^ 2)

/-- Der zu `T` gehörende Impuls `P_T = −a³`, also das negative räumliche Volumen. -/
def volMom (a : ℝ) : ℝ := -(a ^ 3)

/-- Poisson-Klammer zweier Phasenraumfunktionen, ausgedrückt in ihren vier partiellen
Ableitungen an einem Punkt. -/
def bracket (d₁f d₂f d₁g d₂g : ℝ) : ℝ := d₁f * d₂g - d₂f * d₁g

/-- `∂T/∂a = −2p_a/(3a³)`. Quotientenregel mit konstantem Zähler; benötigt `a ≠ 0`. -/
theorem hasDerivAt_yorkTime_fst (pa a : ℝ) (ha : a ≠ 0) :
    HasDerivAt (fun z => yorkTime z pa) (-2 * pa / (3 * a ^ 3)) a := by
  unfold yorkTime
  have hne : (3 : ℝ) * a ^ 2 ≠ 0 :=
    mul_ne_zero (show (3 : ℝ) ≠ 0 by norm_num) (pow_ne_zero 2 ha)
  have hden : HasDerivAt (fun z : ℝ => 3 * z ^ 2) (6 * a) a := by
    convert ((hasDerivAt_id a).pow 2).const_mul 3 using 1
    simp only [id_eq]
    ring
  convert (hasDerivAt_const a pa).div hden hne using 1
  try simp only [id_eq]
  field_simp [ha]
  ring

/-- `∂T/∂p_a = 1/(3a²)`: die York-Zeit ist linear im Impuls. -/
theorem hasDerivAt_yorkTime_snd (pa a : ℝ) (ha : a ≠ 0) :
    HasDerivAt (fun z => yorkTime a z) (1 / (3 * a ^ 2)) pa := by
  unfold yorkTime
  have hne : (3 : ℝ) * a ^ 2 ≠ 0 :=
    mul_ne_zero (show (3 : ℝ) ≠ 0 by norm_num) (pow_ne_zero 2 ha)
  convert ((hasDerivAt_id pa).div (hasDerivAt_const pa (3 * a ^ 2)) hne) using 1
  try simp only [id_eq]
  field_simp [hne]
  ring

/-- `∂P_T/∂a = −3a²`. -/
theorem hasDerivAt_volMom_fst (a : ℝ) :
    HasDerivAt (fun z => volMom z) (-3 * a ^ 2) a := by
  unfold volMom
  convert ((hasDerivAt_pow 3 a).const_mul (-1)) using 1
  · funext z
    ring
  · norm_num

/-- `∂P_T/∂p_a = 0`: das Volumen hängt nicht vom Impuls ab. -/
theorem hasDerivAt_volMom_snd (pa a : ℝ) :
    HasDerivAt (fun _z => volMom a) 0 pa := hasDerivAt_const pa (volMom a)

/-- **Das kanonische York-Paar.** Für `a ≠ 0` ist
`{T, P_T} = (∂T/∂a)(∂P_T/∂p_a) − (∂T/∂p_a)(∂P_T/∂a) = +1`.
Der Beweis setzt die vier oben nachgerechneten partiellen Ableitungen ein; es wird keine
Kanonizität vorausgesetzt, sondern aus `T = p_a/(3a²)` und `P_T = −a³` abgeleitet. -/
theorem york_pair_canonical (pa a : ℝ) (ha : a ≠ 0) :
    bracket (-2 * pa / (3 * a ^ 3)) (1 / (3 * a ^ 2)) (-3 * a ^ 2) 0 = 1 := by
  unfold bracket
  field_simp [ha]
  ring

/-! ### Die Algebra des Hamilton-Zwangs -/

/-- Der Radikand `A(T, φ) = 12πG T² − 2V(φ)`. -/
def radikand (G : ℝ) (T : ℝ) (V : ℝ → ℝ) (phi : ℝ) : ℝ :=
  12 * Real.pi * G * T ^ 2 - 2 * V phi

/-- Der Hamilton-Zwang in York-Variablen, Gl. (28.9):
`C = −p_φ²/(2P_T) + 6πG T² P_T − P_T V(φ)`. -/
noncomputable def hamConstraint (G : ℝ) (T P pphi : ℝ) (V : ℝ → ℝ) (phi : ℝ) : ℝ :=
  -pphi ^ 2 / (2 * P) + 6 * Real.pi * G * T ^ 2 * P - P * V phi

/-- Multiplikation des Zwangs mit `2P_T` gibt `P_T²·A − p_φ²`; das ist die im Manuskript
angegebene Umformung, nach der der Zwang quadratisch in `P_T` ist. -/
theorem two_mul_P_hamConstraint (G : ℝ) (T P pphi : ℝ) (V : ℝ → ℝ) (phi : ℝ) (hP : P ≠ 0) :
    2 * P * hamConstraint G T P pphi V phi
      = P ^ 2 * radikand G T V phi - pphi ^ 2 := by
  unfold hamConstraint radikand
  have h2 : (2 : ℝ) * P ≠ 0 := mul_ne_zero two_ne_zero hP
  field_simp [h2]
  ring

/-- Äquivalenz des Zwangs mit der quadratischen Form: für `P_T ≠ 0` ist `C = 0` genau dann,
wenn `P_T²·A = p_φ²`. -/
theorem hamConstraint_eq_zero_iff (G : ℝ) (T P pphi : ℝ) (V : ℝ → ℝ) (phi : ℝ) (hP : P ≠ 0) :
    hamConstraint G T P pphi V phi = 0 ↔ P ^ 2 * radikand G T V phi = pphi ^ 2 := by
  have h := two_mul_P_hamConstraint G T P pphi V phi hP
  constructor
  · intro hC
    rw [hC, mul_zero] at h
    linarith
  · intro hEq
    apply mul_left_cancel₀ (mul_ne_zero two_ne_zero hP)
    rw [h, hEq, sub_self, mul_zero]

/-- **Algebraischer Kern der Zweigverifikation.** Für `s ≠ 0`, `p ≠ 0` und `s² = 2(c − V)`
erfüllt `P = −p/s` die Gleichung `−p²/(2P) + cP − PV = 0`. Nach Beseitigung der Nenner bleibt
`p·(s² − 2(c − V)) = 0` übrig - die Rechnung ist also exakt äquivalent zur Quadratwurzelwahl. -/
theorem branch_core (c V p s : ℝ) (hs : s ≠ 0) (hp : p ≠ 0) (hs2 : s * s = 2 * (c - V)) :
    -p ^ 2 / (2 * (-p / s)) + c * (-p / s) - (-p / s) * V = 0 := by
  have key : -p ^ 2 / (2 * (-p / s)) = p * s / 2 := by field_simp [hs, hp]
  rw [key]
  apply mul_left_cancel₀ (mul_ne_zero two_ne_zero hs)
  have h : 2 * s * (p * s / 2 + c * (-p / s) - (-p / s) * V) = p * (s * s - 2 * (c - V)) := by
    field_simp [hs, hp]
    ring
  rw [h, hs2, sub_self, mul_zero, mul_zero]

/-- **Zweigverifikation.** Ist `s > 0` mit `s² = A(T,φ)` und `P_T = −p_φ/s`, so erfüllt
dieser `P_T` den Zwang: `C(T, P_T, p_φ) = 0`. Die Rechnung benutzt nur `s² = A`; die
Wurzelfunktion kommt erst im Korollar ins Spiel. -/
theorem hamConstraint_eq_zero_branch (G : ℝ) (T P phi pphi s : ℝ) (V : ℝ → ℝ)
    (hs : 0 < s) (hs2 : s * s = radikand G T V phi) (hP : P = -pphi / s)
    (hp : pphi ≠ 0) :
    hamConstraint G T P pphi V phi = 0 := by
  unfold hamConstraint
  rw [hP]
  refine branch_core (6 * Real.pi * G * T ^ 2) (V phi) pphi s hs.ne' hp ?_
  have h : radikand G T V phi = 2 * (6 * Real.pi * G * T ^ 2 - V phi) := by
    unfold radikand
    ring
  rw [hs2, h]

/-- Korollar mit der positiven Quadratwurzel: auf dem Zweig `P_T = −p_φ/√A` mit `A > 0`
verschwindet der Zwang. Das ist die im Manuskript gewählte Wurzel (`P_T < 0` für `p_φ > 0`,
und umgekehrt). -/
theorem hamConstraint_eq_zero_sqrt (G : ℝ) (T phi pphi : ℝ) (V : ℝ → ℝ)
    (hA : 0 < radikand G T V phi) (hp : pphi ≠ 0) :
    hamConstraint G T (-pphi / Real.sqrt (radikand G T V phi)) pphi V phi = 0 :=
  hamConstraint_eq_zero_branch G T (-pphi / Real.sqrt (radikand G T V phi)) phi pphi
    (Real.sqrt (radikand G T V phi)) V
    (Real.sqrt_pos.mpr hA) (Real.mul_self_sqrt (le_of_lt hA)) rfl hp

/-! ### Die reduzierte Hamilton-Funktion und ihre unabhängige Differentiation -/

/-- Der reduzierte Hamiltonian auf dem Impulszweig `σ = sgn(p_φ)`:
`H_σ = σ·p_φ/√A`. Für `σ = −1` ist das die im Manuskript benutzte Form `H_T = −p_φ/√A`. -/
noncomputable def hamBranch (G σ T : ℝ) (V : ℝ → ℝ) (phi pphi : ℝ) : ℝ :=
  σ * pphi / Real.sqrt (radikand G T V phi)

/-- `∂H_σ/∂p_φ = σ/√A`: der Hamiltonian ist linear im Impuls, die Ableitung ist eine
Konstante - hier tatsächlich als `HasDerivAt` nachgewiesen. -/
theorem hasDerivAt_hamBranch_pphi (G σ T phi pphi : ℝ) (V : ℝ → ℝ)
    (hA : 0 < radikand G T V phi) :
    HasDerivAt (fun z => hamBranch G σ T V phi z) (σ / Real.sqrt (radikand G T V phi)) pphi := by
  unfold hamBranch
  have hs : Real.sqrt (radikand G T V phi) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hA)
  convert ((hasDerivAt_id pphi).const_mul σ).div (hasDerivAt_const pphi _) hs using 1
  try simp only [id_eq, mul_zero, sub_zero, mul_one]
  field_simp [hs]

/-- `∂A/∂φ = −2V'(φ)`: der Geometrieteil von `A` ist φ-unabhängig. -/
theorem hasDerivAt_radikand_phi (G T phi : ℝ) (V V' : ℝ → ℝ) (hV : ∀ x, HasDerivAt V (V' x) x) :
    HasDerivAt (fun z => radikand G T V z) (-2 * V' phi) phi := by
  unfold radikand
  convert (hasDerivAt_const phi (12 * Real.pi * G * T ^ 2)).sub ((hV phi).const_mul 2) using 1
  ring

/-- **Ableitung des Hamiltonians nach dem Feld.** Mit `A^{3/2} := A·√A` gilt
`∂H_σ/∂φ = σ·p_φ·V'(φ)/(A·√A)`. Bewiesen mit Ketten- und Quotientenregel für
`√A` und `1/√A`, also durch unabhängige Differentiation - das ist die Stelle, an der das
Vorzeichen in [RoserValentini2014] Gl. (3.11) überprüft wird. -/
theorem hasDerivAt_hamBranch_phi (G σ T phi pphi : ℝ) (V V' : ℝ → ℝ)
    (hV : ∀ x, HasDerivAt V (V' x) x) (hA : 0 < radikand G T V phi) :
    HasDerivAt (fun z => hamBranch G σ T V z pphi)
      (σ * pphi * V' phi / (radikand G T V phi * Real.sqrt (radikand G T V phi))) phi := by
  have hA' : HasDerivAt (fun z => radikand G T V z) (-2 * V' phi) phi :=
    hasDerivAt_radikand_phi G T phi V V' hV
  have hs : Real.sqrt (radikand G T V phi) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hA)
  have hAne : radikand G T V phi ≠ 0 := ne_of_gt hA
  -- `√A` ableiten, dann `1/√A`, dann mit der Konstanten `σ·p_φ` multiplizieren.
  have hsqrt : HasDerivAt (fun z => Real.sqrt (radikand G T V z))
      (-2 * V' phi / (2 * Real.sqrt (radikand G T V phi))) phi := (hA').sqrt hAne
  have hinv : HasDerivAt (fun z => (Real.sqrt (radikand G T V z))⁻¹)
      (-(-2 * V' phi / (2 * Real.sqrt (radikand G T V phi)))
        / (Real.sqrt (radikand G T V phi) ^ 2)) phi := hsqrt.inv hs
  have heq : (fun z => hamBranch G σ T V z pphi)
      = fun z => (σ * pphi) * (Real.sqrt (radikand G T V z))⁻¹ := by
    funext z
    unfold hamBranch
    field_simp [ne_of_gt (Real.sqrt_pos.mpr hA)]
  rw [heq]
  convert (hinv.const_mul (σ * pphi)) using 1
  field_simp [hs, hA.ne']
  simp only [Real.sq_sqrt (le_of_lt hA)]

/-- **Hamilton-Gleichungen auf dem Zweig `σ`**, Gl. (28.9): mit den beiden eben bewiesenen
partiellen Ableitungen gilt `φ' = ∂H/∂p_φ = σ/√A` und
`p_φ' = −∂H/∂φ = −σ·p_φ·V'(φ)/(A·√A)`. Auf dem Manuskript-Zweig `σ = −1` wird daraus
`φ' = −1/√A` und `p_φ' = +p_φ V'/(A·√A)` - das Plus ist Ergebnis der Differentiation. -/
theorem york_hamilton_equations (G σ T phi pphi : ℝ) (V V' : ℝ → ℝ)
    (hV : ∀ x, HasDerivAt V (V' x) x) (hA : 0 < radikand G T V phi) :
    HasDerivAt (fun z => hamBranch G σ T V phi z) (σ / Real.sqrt (radikand G T V phi)) pphi ∧
      HasDerivAt (fun z => hamBranch G σ T V z pphi)
        (σ * pphi * V' phi / (radikand G T V phi * Real.sqrt (radikand G T V phi))) phi :=
  ⟨hasDerivAt_hamBranch_pphi G σ T phi pphi V hA,
    hasDerivAt_hamBranch_phi G σ T phi pphi V V' hV hA⟩

/-- **Vorzeichen auf dem Manuskriptzweig `σ = −1`.** Die rechte Seite der
`p_φ`-Bewegungsgleichung ist `−∂H/∂φ`, also `+p_φ V'(φ)/(A·√A)`. Das Plus ist damit
nachgerechnet; es folgt aus `∂A/∂φ = −2V'(φ)` und der Kettenregel, nicht aus einer Konvention. -/
theorem branch_minus_momentum_force (G T phi pphi : ℝ) (V V' : ℝ → ℝ) :
    -((-1 : ℝ) * pphi * V' phi / (radikand G T V phi * Real.sqrt (radikand G T V phi)))
      = pphi * V' phi / (radikand G T V phi * Real.sqrt (radikand G T V phi)) := by
  ring

/-! ### Das freie masselose Feld: Lösung und Geschwindigkeit -/

/-- `√(c·T²) = √c·|T|` für `c ≥ 0`. Wird im freien Fall mit `c = 12πG` gebraucht. -/
theorem sqrt_mul_sq (c T : ℝ) (hc : 0 ≤ c) :
    Real.sqrt (c * T ^ 2) = Real.sqrt c * |T| := by
  rw [Real.sqrt_mul hc, Real.sqrt_sq_eq_abs]

/-- Die im Manuskript angegebene Lösung der freien York-Dynamik,
`φ(T) = ln(−T)/√(12πG) + c`, hat für `T < 0` die Ableitung `1/(√(12πG)·T)`, also ein
negatives Vorzeichen auf dem Expansionszweig (`T < 0`). -/
theorem hasDerivAt_free_solution (G c T : ℝ) (hG : 0 < G) (hT : T < 0) :
    HasDerivAt (fun z => Real.log (-z) / Real.sqrt (12 * Real.pi * G) + c)
      (1 / (Real.sqrt (12 * Real.pi * G) * T)) T := by
  have hT' : -T ≠ 0 := ne_of_gt (neg_pos.mpr hT)
  have hc : 0 < Real.sqrt (12 * Real.pi * G) := Real.sqrt_pos.mpr (by positivity)
  -- Als `Konstante * log(-z)` geschrieben, vermeidet die Rechnung die Quotientenregel und
  -- damit Quadrate der Wurzel.
  have hlog : HasDerivAt (fun z : ℝ => Real.log (-z)) (-1 / -T) T :=
    ((hasDerivAt_id T).neg).log hT'
  have h1 : HasDerivAt (fun z : ℝ => (1 / Real.sqrt (12 * Real.pi * G)) * Real.log (-z))
      ((1 / Real.sqrt (12 * Real.pi * G)) * (-1 / -T)) T := hlog.const_mul _
  have heq : (fun z => Real.log (-z) / Real.sqrt (12 * Real.pi * G) + c)
      = fun z => (1 / Real.sqrt (12 * Real.pi * G)) * Real.log (-z) + c := by
    funext z
    field_simp [hc.ne']
  rw [heq]
  convert h1.add (hasDerivAt_const T c) using 1
  field_simp [hT', hc.ne']
  ring

/-- Die zweite Form der freien Geschwindigkeit, wie sie im Manuskript für den Zweig
`p_φ < 0` dasteht: `∂H/∂p_φ = −1/√(12πGT²)`, ausgewertet für `V = 0`. Zusammen mit dem
vorigen Satz und `free_branch_velocity` ist die Lösung damit verifiziert. -/
theorem free_branch_dH_dphi (G T phi pphi : ℝ) :
    hamBranch G (-1) T (fun _ => 0) phi pphi = -pphi / Real.sqrt (12 * Real.pi * G * T ^ 2) := by
  unfold hamBranch radikand
  rw [show (-1 : ℝ) * pphi = -pphi by ring]
  simp only [mul_zero, sub_zero]


/-- Identifikation der beiden im Manuskript nebeneinander stehenden Formen der freien
Geschwindigkeit: `-1/√(12πGT²) = 1/(√(12πG)·T)` für `T < 0`. Zusammen mit dem vorigen Satz
hat die angegebene Lösung damit exakt die aus `∂H/∂p_φ` folgende Geschwindigkeit. -/
theorem free_branch_velocity (G T : ℝ) (hG : 0 < G) (hT : T < 0) :
    -1 / Real.sqrt (12 * Real.pi * G * T ^ 2) = 1 / (Real.sqrt (12 * Real.pi * G) * T) := by
  have hc : 0 ≤ (12 : ℝ) * Real.pi * G := le_of_lt (by positivity)
  rw [sqrt_mul_sq ((12 : ℝ) * Real.pi * G) T hc, abs_of_neg hT]
  have hT' : T ≠ 0 := ne_of_lt hT
  have hs : Real.sqrt (12 * Real.pi * G) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr (by positivity))
  field_simp [hT', hs]

end

end YorkCanonicalPair

#print axioms YorkCanonicalPair.york_pair_canonical
#print axioms YorkCanonicalPair.two_mul_P_hamConstraint
#print axioms YorkCanonicalPair.hamConstraint_eq_zero_iff
#print axioms YorkCanonicalPair.branch_core
#print axioms YorkCanonicalPair.hamConstraint_eq_zero_branch
#print axioms YorkCanonicalPair.hamConstraint_eq_zero_sqrt
#print axioms YorkCanonicalPair.hasDerivAt_hamBranch_pphi
#print axioms YorkCanonicalPair.hasDerivAt_hamBranch_phi
#print axioms YorkCanonicalPair.york_hamilton_equations
#print axioms YorkCanonicalPair.hasDerivAt_free_solution
#print axioms YorkCanonicalPair.free_branch_velocity
#print axioms YorkCanonicalPair.free_branch_dH_dphi
