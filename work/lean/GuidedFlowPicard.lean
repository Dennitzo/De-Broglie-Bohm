import Mathlib

/-!
Allgemeine lokale Existenz und Eindeutigkeit der gefuehrten Bahn

Formalisiert die in 18.6 als analytisches Picard-Lindelöf-Lemma ausgewiesene
Aussage fuer ein *beliebiges* glattes gefuehrtes Feld, ohne Bezug auf das
planare Referenzfeld von 18.7.

Sei V : E -> E die gefuehrte Geschwindigkeit V^A = M^{AB} S_B des endlichen
homogenen Geometrie-Feld-Modells (E = R^(n+1), Koordinaten (alpha, phi^I)) und
sei V am Startpunkt x_0 von Klasse C^1. Dann existiert ein eps > 0 und eine
Funktion gamma : R -> E mit gamma 0 = x_0, die auf Icc (-eps) eps die
Picard-Integralgleichung

    gamma t = x_0 + integral_{s in 0..t} V (gamma s) ds

erfuellt; das ist die lokale Loesung der gefuehrten Gleichung qdot = V(q).
Zwei Loesungen, die in derselben Picard-Kapsel Fixpunkte des Picard-Iterators
sind, stimmen ueberein, weil eine Iterierte dieses Iterators kontrahierend ist.

Abgrenzung. Bewiesen wird lokale Existenz (in Integralgleichungsform) und
lokale Eindeutigkeit innerhalb der Picard-Kapsel. Nicht bewiesen und nicht
behauptet: globale Existenz oder zeitliche Vollstaendigkeit, ein Verlassen-
Verbot vorgegebener Gebiete, die Existenz einer Bahn an einem Knoten (dort ist
V^A = M^{AB} S_B nicht definiert), und eine separate Lean-Fassung fuer jede
konkrete Feldzahl n (die Existenzaussage gilt fuer jeden Banach-Raum E; die
endlichen Modelle sind Instanzen). Die Aequivalenz "Integralgleichung <->
Differentialgleichung qdot = V(q)" ist analytisch
(intervalIntegral.integral_hasDerivAt_right plus Stetigkeit von
s |-> V (gamma s)); sie wird hier nicht als eigener Lean-Satz ausgewiesen.
-/

open Set Metric

namespace GuidedFlowPicard

/-- Fixpunkte einer Abbildung sind Fixpunkte jeder Iterierten. -/
theorem fixed_of_iterate {X : Type} {f : X → X} {p : X} (hp : f p = p) :
    ∀ n : ℕ, (f^[n]) p = p
  | 0 => rfl
  | (k + 1) => by rw [Function.iterate_succ_apply', fixed_of_iterate hp k, hp]

/-- Eindeutigkeitskern: Ist eine Iterierte kontrahierend mit Faktor < 1, so hat
die Abbildung hoechstens einen Fixpunkt. -/
theorem unique_of_contracting_iterate_dist {X : Type} [MetricSpace X] {f : X → X}
    {n : ℕ} {c : NNReal} (hc1 : c < 1)
    (hc : ∀ x y, dist ((f^[n]) x) ((f^[n]) y) ≤ c * dist x y)
    {p q : X} (hp : f p = p) (hq : f q = q) : p = q := by
  have hp' : (f^[n]) p = p := fixed_of_iterate hp n
  have hq' : (f^[n]) q = q := fixed_of_iterate hq n
  have h := hc p q
  rw [hp', hq'] at h
  by_contra hne
  have hpos : (0 : ℝ) < dist p q := dist_pos.mpr hne
  have h2 : (dist p q : ℝ) ≤ (c : ℝ) * dist p q := by exact_mod_cast h
  have hcr : (c : ℝ) < 1 := by exact_mod_cast hc1
  nlinarith

/-- Zwei Fixpunkte des Picard-Iterators derselben Kapsel sind gleich. -/
theorem unique_picard_fixed_pt {E : Type} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [CompleteSpace E] {f : ℝ → E → E} {tmin tmax : ℝ} {t₀ : Set.Icc tmin tmax}
    {x₀ x : E} {a r L K : NNReal}
    (hf : IsPicardLindelof f t₀ x₀ a r L K) (hx : x ∈ closedBall x₀ ↑r)
    (α β : ODE.FunSpace t₀ x₀ r L)
    (ha : Function.IsFixedPt (ODE.FunSpace.next hf hx) α)
    (hb : Function.IsFixedPt (ODE.FunSpace.next hf hx) β) : α = β := by
  obtain ⟨n, c, hc⟩ := ODE.FunSpace.exists_contractingWith_iterate_next hf
  have hctr := hc x hx
  exact unique_of_contracting_iterate_dist hctr.1
    (fun y z => LipschitzWith.dist_le_mul hctr.2 y z) ha hb

theorem mem_closedBall_self {E : Type} [PseudoMetricSpace E] (x₀ : E) (r : NNReal)
    (h : (0 : ℝ) ≤ (r : ℝ)) : x₀ ∈ closedBall x₀ r := by
  rw [mem_closedBall, dist_self]; exact h

/-- **Lokale Existenz der gefuehrten Bahn.** Ist die gefuehrte Geschwindigkeit
am Startpunkt C^1, so gibt es eps > 0 und eine Funktion gamma mit gamma 0 = x_0,
die auf Icc (-eps) eps die Picard-Integralgleichung erfuellt. -/
theorem exists_guided_integral {E : Type} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [CompleteSpace E] (V : E → E) (x₀ : E) (hV : ContDiffAt ℝ 1 V x₀) :
    ∃ (ε : ℝ) (hε : 0 < ε) (γ : ℝ → E),
      γ 0 = x₀ ∧
      ∀ t ∈ Icc (-ε) ε, γ t = x₀ + ∫ s in (0 : ℝ)..t, V (γ s) := by
  obtain ⟨μ, _hμ, a, r, L, K, hr, h⟩ := IsPicardLindelof.of_contDiffAt_one hV
  obtain ⟨γ, hγ0, hγ⟩ :=
    IsPicardLindelof.exists_eq_forall_mem_Icc_eq_picard (h 0)
      (mem_closedBall_self x₀ r (NNReal.coe_nonneg _))
  refine ⟨μ, _hμ, γ, ?_, ?_⟩
  · simpa using hγ0
  · intro t ht
    have hco : t ∈ Set.Icc ((0 : ℝ) - μ) (0 + μ) := by simpa using ht
    have := hγ t hco
    rw [ODE.picard_apply] at this
    simpa using this

end GuidedFlowPicard

#print axioms GuidedFlowPicard.fixed_of_iterate
#print axioms GuidedFlowPicard.unique_of_contracting_iterate_dist
#print axioms GuidedFlowPicard.unique_picard_fixed_pt
#print axioms GuidedFlowPicard.exists_guided_integral
