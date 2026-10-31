import Mathlib

/-!
Phase 5: Brueckensatze.

Fuenf Verbindungen, die im Manuskript analytisch begruendet, aber nicht als eigene
formale Saetze ausgewiesen waren:

0. Endliche Aggregation: Fuer eine beliebige endliche Indexmenge ist die Ableitung
   der Summe die Summe der Ableitungen; zusammen mit der Potenzregel ergibt das die
   Ableitung einer aggregierten Dichte.

1. Integralgleichung und differentielle Fuehrung: Fuer stetige Besetzungsfelder ist
   die gefuehrte Differentialgleichung gamma'(t) = V(gamma(t)) gleichwertig zur
   Picard-Integralgleichung; beide Richtungen sind als Saetze formuliert.

2. Constraintfortpflanzung auf zusammenhaengendem Intervall: Aus A' = H und
   C' = -3 H C folgt auf Icc 0 T die gewichtete Konstanz (exp(3A) C)' = 0, und
   verschwindet C am Anfangspunkt, so verschwindet C auf dem ganzen Intervall.

3. Instanziierung V = M^{-1} dS: Fuer die invertierbare Konfigurationsmetrik G ist
   V = G^{-1} mulVec dS die eindeutige Loesung von G mulVec V = dS; das ist die
   musikalische Abbildung des Gradienten im endlichen Modell.

4. Stromidentitaet unter der Zustandsgleichung: Erfuellt der Strom die
   Kontinuitaetsformel mit einer rechten Seite proportional zum Imaginaerteil von
   Psi * (C Psi), so folgt aus C Psi = 0 punktweise verschwindende Divergenz.

Abgrenzung: Es sind algebraische und ableitungsanalytische Kernsaetze fuer das
endliche Modell; keine Existenz- oder Eindeutigkeitsresultate jenseits dessen, was
aus den angegebenen Voraussetzungen folgt, keine Randbedingungen und keine
unendlichdimensionale Theorie.
-/

open Set
open scoped BigOperators

namespace PhaseFiveBridges

/-! 0. Endliche Aggregation: Summenregel ueber beliebige endliche Mengen. -/

theorem hasDerivAt_finset_sum {ι : Type} [DecidableEq ι] (s : Finset ι) {q : ι → ℝ → ℝ} {v : ι → ℝ}
    {t : ℝ} (hq : ∀ i ∈ s, HasDerivAt (q i) (v i) t) :
    HasDerivAt (fun x => Finset.sum s fun i => q i x) (Finset.sum s v) t := by
  induction s using Finset.induction_on with
  | empty => simpa using hasDerivAt_const t (0 : ℝ)
  | @insert a s hasa ih =>
      have h1 := hq a (Finset.mem_insert_self a s)
      have h2 := ih (fun i hi => hq i (Finset.mem_insert_of_mem hi))
      simp only [Finset.sum_insert hasa]
      exact h1.add h2

/-- Aggregierte Dichte: Summe der Quadrate endlicher Felder. -/
theorem hasDerivAt_finset_sum_sq {ι : Type} [DecidableEq ι] (s : Finset ι)
    {q : ι → ℝ → ℝ} {v : ι → ℝ} {t : ℝ} (hq : ∀ i ∈ s, HasDerivAt (q i) (v i) t) :
    HasDerivAt (fun x => Finset.sum s fun i => (q i x) ^ 2)
      (Finset.sum s fun i => 2 * q i t * v i) t :=
  hasDerivAt_finset_sum s fun i hi => by
    convert (hq i hi).pow 2 using 1
    simp [pow_one]

/-! 1. Integralgleichung und differentielle Fuehrung (endliches Modell). -/

theorem guided_hasDerivAt_of_integral
    {n : ℕ} {V : (Fin n → ℝ) → Fin n → ℝ} {γ : ℝ → Fin n → ℝ} {t s : ℝ}
    (hVγ : Continuous (fun r => V (γ r)))
    (hInt : ∀ u, γ u = γ t + ∫ r in t..u, V (γ r)) :
    HasDerivAt γ (V (γ s)) s := by
  have h := intervalIntegral.integral_hasDerivAt_right
    (hVγ.intervalIntegrable t s)
    (hVγ.stronglyMeasurableAtFilter MeasureTheory.volume (nhds s))
    hVγ.continuousAt
  have hγfun : γ = fun u => γ t + ∫ r in t..u, V (γ r) :=
    (funext hInt : γ = fun u => γ t + ∫ r in t..u, V (γ r))
  rw [hγfun]
  convert ((hasDerivAt_const s (γ t)).add h) using 1
  rw [zero_add]
  exact congrArg V (congrFun hγfun s).symm

theorem guided_deriv_of_integral
    {n : ℕ} {V : (Fin n → ℝ) → Fin n → ℝ} {γ : ℝ → Fin n → ℝ} {t s : ℝ}
    (hVγ : Continuous (fun r => V (γ r)))
    (hInt : ∀ u, γ u = γ t + ∫ r in t..u, V (γ r)) :
    deriv γ s = V (γ s) :=
  (guided_hasDerivAt_of_integral hVγ hInt).deriv

theorem integral_of_guided_hasDerivAt
    {n : ℕ} {V : (Fin n → ℝ) → Fin n → ℝ} {γ : ℝ → Fin n → ℝ} {a b : ℝ}
    (hVγ : Continuous (fun r => V (γ r)))
    (hdiff : ∀ t, HasDerivAt γ (V (γ t)) t)
    (hab : a ≤ b) :
    γ b = γ a + ∫ t in a..b, V (γ t) := by
  have hcont : Continuous γ :=
    continuous_iff_continuousAt.mpr fun t => (hdiff t).continuousAt
  have h := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (f' := fun t => V (γ t)) hab
    hcont.continuousOn (fun t _ => hdiff t) (hVγ.intervalIntegrable a b)
  rw [h]
  abel

/-! 2. Constraintfortpflanzung auf zusammenhaengendem Intervall. -/

theorem weighted_constraint_const_of_diff
    {A C : ℝ → ℝ} {H : ℝ → ℝ} {T : ℝ}
    (hA : ∀ t ∈ Icc 0 T, HasDerivAt A (H t) t)
    (hC : ∀ t ∈ Icc 0 T, HasDerivAt C (-3 * H t * C t) t)
    (hcont : ContinuousOn (fun t => Real.exp (3 * A t) * C t) (Icc 0 T)) :
    ∀ t ∈ Icc 0 T, Real.exp (3 * A t) * C t = Real.exp (3 * A 0) * C 0 := by
  intro t ht
  have hzero : ∀ s ∈ Ioo 0 t, HasDerivAt (fun u => Real.exp (3 * A u) * C u) 0 s := by
    intro s hs
    have hst : s ∈ Icc 0 T := ⟨le_of_lt hs.1, le_of_lt (lt_of_lt_of_le hs.2 ht.2)⟩
    have h1 : HasDerivAt (fun u => Real.exp (3 * A u)) (Real.exp (3 * A s) * (3 * H s)) s :=
      ((hA s hst).const_mul 3).exp
    have hz : HasDerivAt (fun u => Real.exp (3 * A u) * C u)
        (Real.exp (3 * A s) * (3 * H s) * C s + Real.exp (3 * A s) * (-3 * H s * C s)) s :=
      h1.mul (hC s hst)
    convert hz using 1
    ring
  have hcont' : ContinuousOn (fun u => Real.exp (3 * A u) * C u) (Icc 0 t) :=
    hcont.mono fun _ hx => ⟨hx.1, le_trans hx.2 ht.2⟩
  have hkey := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (f' := fun _ => (0 : ℝ)) ht.1
    hcont' hzero
    (ContinuousOn.intervalIntegrable_of_Icc (μ := MeasureTheory.volume) ht.1
      (continuousOn_const (s := Icc (0 : ℝ) t)) :
      IntervalIntegrable (fun _ : ℝ => (0 : ℝ)) MeasureTheory.volume 0 t)
  simp only [intervalIntegral.integral_zero] at hkey
  linarith

theorem constraint_eq_zero_of_zero_at_start
    {A C : ℝ → ℝ} {H : ℝ → ℝ} {T : ℝ}
    (hA : ∀ t ∈ Icc 0 T, HasDerivAt A (H t) t)
    (hC : ∀ t ∈ Icc 0 T, HasDerivAt C (-3 * H t * C t) t)
    (hcont : ContinuousOn (fun t => Real.exp (3 * A t) * C t) (Icc 0 T))
    (h0 : C 0 = 0) :
    ∀ t ∈ Icc 0 T, C t = 0 := by
  intro t ht
  have h := weighted_constraint_const_of_diff hA hC hcont t ht
  rw [h0, mul_zero] at h
  rcases mul_eq_zero.mp h with hexp | hCt
  · exact absurd hexp (Real.exp_ne_zero _)
  · exact hCt

/-! 3. Instanziierung V = G^{-1} mulVec dS (musikalische Abbildung, endliche Dimension). -/

theorem inverse_metric_velocity
    {n : ℕ} (G : Matrix (Fin n) (Fin n) ℝ) (hG : IsUnit G.det) (dS : Fin n → ℝ) :
    G.mulVec (G⁻¹.mulVec dS) = dS := by
  calc G.mulVec (G⁻¹.mulVec dS)
      = (G * G⁻¹).mulVec dS := by rw [← Matrix.mulVec_mulVec]
    _ = (1 : Matrix (Fin n) (Fin n) ℝ).mulVec dS := by
        rw [show G * G⁻¹ = (1 : Matrix (Fin n) (Fin n) ℝ) from Matrix.mul_nonsing_inv G hG]
    _ = dS := by simp

theorem inverse_metric_velocity_unique
    {n : ℕ} (G : Matrix (Fin n) (Fin n) ℝ) (hG : IsUnit G.det) (dS : Fin n → ℝ)
    (V W : Fin n → ℝ) (hV : G.mulVec V = dS) (hW : G.mulVec W = dS) :
    V = W := by
  have hleft : ∀ (U : Fin n → ℝ), G⁻¹.mulVec (G.mulVec U) = U := by
    intro U
    rw [Matrix.mulVec_mulVec,
      show G⁻¹ * G = (1 : Matrix (Fin n) (Fin n) ℝ) from Matrix.nonsing_inv_mul G hG]
    simp
  calc V = G⁻¹.mulVec (G.mulVec V) := (hleft V).symm
    _ = G⁻¹.mulVec (G.mulVec W) := by rw [hV, hW]
    _ = W := hleft W

/-! 4. Stromidentitaet unter der Zustandsgleichung C Psi = 0. -/

theorem current_div_iff_imag
    {μ ħ : ℝ} (hμ : μ ≠ 0) (hħ : ħ ≠ 0)
    {J : ℝ → ℝ} {Ψ : ℝ → ℂ} {CΨ : ℝ → ℂ}
    (hcont : ∀ x, HasDerivAt J (-(2 * μ / ħ) * (Ψ x * CΨ x).im) x) (x : ℝ) :
    deriv J x = 0 ↔ (Ψ x * CΨ x).im = 0 := by
  have h2 : (2 : ℝ) ≠ 0 := by norm_num
  have hc : -(2 * μ / ħ) ≠ 0 :=
    neg_ne_zero.mpr (div_ne_zero (mul_ne_zero h2 hμ) hħ)
  have hx := (hcont x).deriv
  constructor
  · intro hd
    rw [hx] at hd
    rw [mul_eq_zero] at hd
    exact hd.resolve_left hc
  · intro hi
    rw [hx, hi, mul_zero]

theorem current_deriv_vanishes_of_constraint
    {μ ħ : ℝ}
    {J : ℝ → ℝ} {Ψ : ℝ → ℂ} {CΨ : ℝ → ℂ}
    (hcont : ∀ x, HasDerivAt J (-(2 * μ / ħ) * (Ψ x * CΨ x).im) x)
    (hconstraint : ∀ x, CΨ x = 0) (x : ℝ) :
    deriv J x = 0 := by
  have hx := hcont x
  rw [hconstraint x] at hx
  simp at hx
  exact hx.deriv

theorem imag_term_vanishes_of_constraint
    {Ψ : ℝ → ℂ} {CΨ : ℝ → ℂ} (hconstraint : ∀ x, CΨ x = 0) (x : ℝ) :
    (Ψ x * CΨ x).im = 0 := by
  rw [hconstraint x]
  simp

end PhaseFiveBridges
