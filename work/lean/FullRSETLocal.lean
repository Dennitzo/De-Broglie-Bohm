import Mathlib

/-! A local differential subtraction model with stated smooth-difference hypotheses.
This does not construct a Hadamard parametrix or sum angular/frequency modes. -/
namespace FullRSETLocal
noncomputable section
open Filter
open scoped Topology

def differentialExpression (a b c value first second : ℝ) : ℝ :=
  a*second+b*first+c*value

theorem common_subtraction_cancels (a b c f₀ f₁ f₂ g₀ g₁ g₂ h₀ h₁ h₂ : ℝ) :
    differentialExpression a b c (f₀-h₀) (f₁-h₁) (f₂-h₂)-
      differentialExpression a b c (g₀-h₀) (g₁-h₁) (g₂-h₂) =
    differentialExpression a b c (f₀-g₀) (f₁-g₁) (f₂-g₂) := by
  unfold differentialExpression
  ring

theorem actual_first_derivative_subtraction (f g : ℝ → ℝ) (x df dg : ℝ)
    (hf : HasDerivAt f df x) (hg : HasDerivAt g dg x) :
    HasDerivAt (fun t => f t-g t) (df-dg) x := hf.sub hg

theorem actual_second_derivative_subtraction (f g : ℝ → ℝ) (x ddf ddg : ℝ)
    (hf : HasDerivAt (deriv f) ddf x) (hg : HasDerivAt (deriv g) ddg x)
    (hf' : ∀ᶠ t in 𝓝 x, DifferentiableAt ℝ f t)
    (hg' : ∀ᶠ t in 𝓝 x, DifferentiableAt ℝ g t) :
    HasDerivAt (deriv (fun t => f t-g t)) (ddf-ddg) x := by
  apply (hf.sub hg).congr_of_eventuallyEq
  filter_upwards [hf', hg'] with t ht ht'
  exact deriv_sub ht ht'

theorem finite_coincidence_difference {ι : Type*} {l : Filter ι}
    (F₀ F₁ F₂ : ι → ℝ) (a b c f₀ f₁ f₂ : ℝ)
    (h₀ : Tendsto F₀ l (𝓝 f₀)) (h₁ : Tendsto F₁ l (𝓝 f₁))
    (h₂ : Tendsto F₂ l (𝓝 f₂)) :
    Tendsto (fun t => differentialExpression a b c (F₀ t) (F₁ t) (F₂ t)) l
      (𝓝 (differentialExpression a b c f₀ f₁ f₂)) := by
  exact ((h₂.const_mul a).add (h₁.const_mul b)).add (h₀.const_mul c)

theorem continuous_jet_has_finite_coincidence (F₀ F₁ F₂ : ℝ → ℝ) (x a b c : ℝ)
    (h₀ : ContinuousAt F₀ x) (h₁ : ContinuousAt F₁ x) (h₂ : ContinuousAt F₂ x) :
    Tendsto (fun t => differentialExpression a b c (F₀ t) (F₁ t) (F₂ t))
      (𝓝 x) (𝓝 (differentialExpression a b c (F₀ x) (F₁ x) (F₂ x))) := by
  exact finite_coincidence_difference F₀ F₁ F₂ a b c (F₀ x) (F₁ x) (F₂ x)
    h₀.tendsto h₁.tendsto h₂.tendsto

theorem common_finite_counterterm_preserves_state_difference (T₁ T₂ C : ℝ) :
    (T₁+C)-(T₂+C)=T₁-T₂ := by ring

end
end FullRSETLocal
