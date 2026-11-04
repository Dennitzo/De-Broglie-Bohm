import Mathlib

open MeasureTheory

/-!
Kurvenaussage zur Periodenbedingung

Ersetzt den früheren Beleg Periods.integral_exact_zero, der nur das Integral
über das entartete Intervall [a,a] bewies. Formalisiert wird die tatsächliche
Aussage: Für eine glatte Kurve c : ℝ → ℝ und eine glatte Funktion F ist
d/dt F(c(t)) = F'(c(t))·c'(t)  (Kettenregel) und
∫_a^b (d/dt) F(c(t)) dt = F(c(b)) - F(c(a))  (Hauptsatz).
Erst aus c(a)=c(b) folgt das Verschwinden des geschlossenen Kurvenintegrals.
Die Parametergrenzen dürfen verschieden sein; die Kurve wird nicht durch ein
entartetes Intervall ersetzt.
-/

namespace PeriodsCurve

/-- Kettenregel entlang einer Kurve: die Ableitung von t ↦ F(c(t)) ist
    F'(c(t)) · c'(t). -/
theorem curve_chain_rule (c : ℝ → ℝ) (F : ℝ → ℝ) (t : ℝ)
    (hF : DifferentiableAt ℝ F (c t)) (hc : DifferentiableAt ℝ c t) :
    deriv (fun s : ℝ => F (c s)) t = deriv F (c t) * deriv c t := by
  have hcomp : HasDerivAt (fun s : ℝ => F (c s)) (deriv F (c t) * deriv c t) t := by
    exact hF.hasDerivAt.comp t hc.hasDerivAt
  exact hcomp.deriv

/-- Hauptsatz entlang einer Kurve: das Integral der Zeitableitung von F(c(t))
    über [a,b] ist die Differenz F(c(b)) - F(c(a)); die Grenzen dürfen
    verschieden sein und die Kurve wird nicht durch ein entartetes Intervall
    ersetzt. Die Stetigkeit der Zeitableitung auf [a,b] ist ausreichend für
    die Integrierbarkeit. -/
theorem curve_fundamental (c : ℝ → ℝ) (F : ℝ → ℝ) (a b : ℝ)
    (hF : Differentiable ℝ F) (hc : Differentiable ℝ c)
    (hcont : ContinuousOn (fun s : ℝ => deriv (fun s : ℝ => F (c s)) s) (Set.uIcc a b)) :
    ∫ t in a..b, deriv (fun s : ℝ => F (c s)) t = F (c b) - F (c a) := by
  have hcompOn : ∀ x ∈ Set.uIcc a b, DifferentiableAt ℝ (fun s : ℝ => F (c s)) x := by
    intro x hx
    have hF_at : DifferentiableAt ℝ F (c x) := hF.differentiableAt
    have hc_at : DifferentiableAt ℝ c x := hc.differentiableAt
    have hcomp_at : DifferentiableAt ℝ (fun s : ℝ => F (c s)) x := hF_at.comp x hc_at
    exact hcomp_at
  have hint : IntervalIntegrable (deriv (fun s : ℝ => F (c s))) volume a b :=
    hcont.intervalIntegrable (μ := volume)
  exact intervalIntegral.integral_deriv_eq_sub hcompOn hint

end PeriodsCurve
