import Mathlib

/-!
Explizite gefuehrte Bahn des potentialfreien Referenzzustands (n Felder)

Formalisiert die lokale Existenz einer gefuehrten Bohmbahn fuer das konkrete
gefuehrte Feld des potentialfreien Mehrfeld-Referenzzustands aus 18.4. Fuer
den ebenen Referenzzustand Psi = exp(-gamma_n alpha/2) * exp i(k.phi - omega alpha)
lautet die Bohmsche Fuehrung (N = 1)

  alpha_dot = A * exp(-3 alpha),   A := hbar * omega * beta,
  phi_dot^I = B_I * exp(-3 alpha), B_I := hbar * k_I * beta.

Die alpha-Komponente ist explizit integrierbar:

  alpha(t) = (1/3) * log(C + 3 A t),   C := exp(3 alpha_0) > 0,

auf dem Intervall {t | C + 3 A t > 0}. Bewiesen wird die gefuehrte
Differentialgleichung alpha_dot(t) = A * exp(-3 alpha(t)) punktweise per
HasDerivAt (Kettenregel, Real.exp_log). Damit existiert fuer dieses konkrete
gefuehrte Feld lokal eine echte Loesung; das maximale Existenzintervall endet
am Rand (C + 3 A t gegen 0). Der allgemeine Fall beliebiger glatter gefuehrter
Felder wird in GuidedFlowPicard.lean mit der Picard-Lindelöf-Maschine (Mathlib:
IsPicardLindelof) für beliebiges n ausgerollt (Publikation 18.8); innerhalb dieser Datei
bleibt es bei der expliziten Bahn des Referenzzustands.

Abgrenzung: Existenz einer konkreten Loesung; keine globale Vollständigkeit,
kein Kollaps-Ausspruch, keine Eindeutigkeits-Aussage ueber alle Loesungen.
-/

namespace GuidedFlowExplicit

/-- Die explizite gefuehrte Skalenfaktor-Bahn. -/
noncomputable def alphaSol (C A : ℝ) (t : ℝ) : ℝ := (3 : ℝ)⁻¹ * Real.log (C + 3 * A * t)

theorem hasDerivAt_inner (C A : ℝ) (t : ℝ) :
    HasDerivAt (fun u => C + (3 * A) * u) (3 * A) t := by
  convert (hasDerivAt_const t C).add ((hasDerivAt_id t).const_mul (3 * A)) using 1
  ring

/-- Rohform der Ableitung (direkt aus Ketten- und konstanter Multiplikation). -/
theorem hasDerivAt_alphaSol_raw (C A : ℝ) (t : ℝ) (hpos : 0 < C + 3 * A * t) :
    HasDerivAt (alphaSol C A) ((3 : ℝ)⁻¹ * ((C + 3 * A * t)⁻¹ * (3 * A))) t := by
  refine ((Real.hasDerivAt_log (ne_of_gt hpos)).comp t (hasDerivAt_inner C A t)).const_mul (3 : ℝ)⁻¹

/-- exp(-3 alphaSol t) = 1 / (C + 3 A t) auf dem Definititsintervall. -/
theorem exp_neg_three_alphaSol (C A : ℝ) (t : ℝ) (hpos : 0 < C + 3 * A * t) :
    Real.exp (-3 * alphaSol C A t) = (C + 3 * A * t)⁻¹ := by
  unfold alphaSol
  have hmul : (-3 : ℝ) * ((3 : ℝ)⁻¹ * Real.log (C + 3 * A * t)) = - Real.log (C + 3 * A * t) := by
    ring
  rw [hmul, Real.exp_neg, Real.exp_log hpos]

/-- Die gefuehrte Differentialgleichung fuer das konkrete Feld:
alpha_dot(t) = A * exp(-3 alpha(t)). -/
theorem hasDerivAt_guided_alpha (C A : ℝ) (t : ℝ) (hpos : 0 < C + 3 * A * t) :
    HasDerivAt (alphaSol C A) (A * Real.exp (-3 * alphaSol C A t)) t := by
  convert hasDerivAt_alphaSol_raw C A t hpos using 1
  rw [exp_neg_three_alphaSol C A t hpos]
  field_simp [hpos.ne']

/-- Die gefuehrte Feldkomponente (proportionale Skalierung der Skalenfaktor-Bahn)
    erfuellt phi_dot = B * exp(-3 alpha). -/
theorem hasDerivAt_guided_phi (C A B : ℝ) (t : ℝ) (hpos : 0 < C + 3 * A * t) :
    HasDerivAt (fun u => (B / A) * alphaSol C A u)
      ((B / A) * (A * Real.exp (-3 * alphaSol C A t))) t :=
  (hasDerivAt_guided_alpha C A t hpos).const_mul (B / A)

end GuidedFlowExplicit

#print axioms GuidedFlowExplicit.hasDerivAt_inner
#print axioms GuidedFlowExplicit.hasDerivAt_alphaSol_raw
#print axioms GuidedFlowExplicit.exp_neg_three_alphaSol
#print axioms GuidedFlowExplicit.hasDerivAt_guided_alpha
#print axioms GuidedFlowExplicit.hasDerivAt_guided_phi
