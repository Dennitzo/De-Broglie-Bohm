import Mathlib

/-! Finite-endpoint obstruction for the regular areal-radius chart.
All existence, coefficient regularity and interval hypotheses are explicit. -/
namespace SWaveEndpoint
noncomputable section
open Filter MeasureTheory
open scoped Topology

def coefficient (alpha r rho w h : ℝ) : ℝ :=
  r/(3*alpha)*(h-1-alpha*w^2-(2*r+9*alpha/r)*w+6*alpha*rho/r^2)

def primitive (P : ℝ → ℝ) (a x : ℝ) : ℝ := ∫ y in a..x, P y

theorem inverse_square_metric_hasDerivAt (u : ℝ → ℝ) (x P : ℝ)
    (hu : HasDerivAt u (P*u x) x) (hu0 : u x ≠ 0) :
    HasDerivAt (fun y => 1/(u y)^2) (-2*P*(1/(u x)^2)) x := by
  convert (hasDerivAt_const x (1 : ℝ)).div (hu.pow 2) (pow_ne_zero 2 hu0) using 1 <;>
    dsimp <;> field_simp <;> ring

theorem primitive_hasDerivAt (P : ℝ → ℝ) (a x : ℝ) (hP : Continuous P) :
    HasDerivAt (primitive P a) (P x) x := by
  exact intervalIntegral.integral_hasDerivAt_right (hP.intervalIntegrable a x)
    (hP.stronglyMeasurableAtFilter volume (𝓝 x)) hP.continuousAt

theorem integrating_factor_hasDerivAt (h A : ℝ → ℝ) (x P : ℝ)
    (hh : HasDerivAt h (-2*P*h x) x) (hA : HasDerivAt A P x) :
    HasDerivAt (fun y => h y*Real.exp (2*A y)) 0 x := by
  convert hh.mul ((hA.const_mul 2).exp) using 1 <;> ring

theorem integrating_factor_constant (h A : ℝ → ℝ) (a b x x0 : ℝ)
    (hx : x ∈ Set.Ioo a b) (hx0 : x0 ∈ Set.Ioo a b)
    (hderiv : ∀ y ∈ Set.Ioo a b, ∃ p,
      HasDerivAt h (-2*p*h y) y ∧ HasDerivAt A p y) :
    h x*Real.exp (2*A x)=h x0*Real.exp (2*A x0) := by
  have hz : ∀ y ∈ Set.Ioo a b,
      HasDerivAt (fun z => h z*Real.exp (2*A z)) 0 y := by
    intro y hy
    obtain ⟨p,hp,hA⟩ := hderiv y hy
    exact integrating_factor_hasDerivAt h A y p hp hA
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun y hy => (hz y hy).differentiableAt.differentiableWithinAt)
    (fun y hy => (hz y hy).deriv) hx hx0

theorem linear_ode_exact_solution (h P : ℝ → ℝ) (a b x x0 : ℝ)
    (hx : x ∈ Set.Ioo a b) (hx0 : x0 ∈ Set.Ioo a b) (hP : Continuous P)
    (hh : ∀ y ∈ Set.Ioo a b, HasDerivAt h (-2*P y*h y) y) :
    h x=h x0*Real.exp (-2*primitive P x0 x) := by
  have hc := integrating_factor_constant h (primitive P x0) a b x x0 hx hx0
    (fun y hy => ⟨P y,hh y hy,primitive_hasDerivAt P x0 y hP⟩)
  have hzero : primitive P x0 x0=0 := by simp [primitive]
  rw [hzero] at hc
  simp only [mul_zero, Real.exp_zero, mul_one] at hc
  have hE := Real.exp_pos (2*primitive P x0 x)
  apply (mul_right_cancel₀ (ne_of_gt hE))
  rw [hc]
  rw [mul_assoc, ← Real.exp_add]
  have he : -2*primitive P x0 x+2*primitive P x0 x=0 := by ring
  rw [he,Real.exp_zero,mul_one]

theorem positive_solution_no_zero (h P : ℝ → ℝ) (a b x x0 : ℝ)
    (hx : x ∈ Set.Ioo a b) (hx0 : x0 ∈ Set.Ioo a b) (hP : Continuous P)
    (hh : ∀ y ∈ Set.Ioo a b, HasDerivAt h (-2*P y*h y) y)
    (h0 : 0 < h x0) : 0 < h x := by
  rw [linear_ode_exact_solution h P a b x x0 hx hx0 hP hh]
  positivity

theorem finite_endpoint_positive_limit (h P : ℝ → ℝ) (a b x0 : ℝ)
    (hx0 : x0 ∈ Set.Ioo a b) (hP : Continuous P)
    (hh : ∀ y ∈ Set.Ioo a b, HasDerivAt h (-2*P y*h y) y)
    (h0 : 0 < h x0) :
    Tendsto h (𝓝[<] b) (𝓝 (h x0*Real.exp (-2*primitive P x0 b))) ∧
      0 < h x0*Real.exp (-2*primitive P x0 b) := by
  constructor
  · have hlim : Tendsto (fun x => h x0*Real.exp (-2*primitive P x0 x))
        (𝓝[<] b) (𝓝 (h x0*Real.exp (-2*primitive P x0 b))) := by
      have hc := ((((primitive_hasDerivAt P x0 b hP).const_mul (-2)).exp).const_mul (h x0)).continuousAt
      exact hc.tendsto.mono_left inf_le_left
    apply hlim.congr'
    filter_upwards [Ioo_mem_nhdsLT (lt_trans hx0.1 hx0.2)] with x hx
    exact (linear_ode_exact_solution h P a b x x0 hx hx0 hP hh).symm
  · positivity


theorem coefficient_tendsto {ι : Type*} {l : Filter ι}
    (r rho w h : ι → ℝ) (alpha rM rhoM wM hM : ℝ)
    (hr : Tendsto r l (𝓝 rM)) (hρ : Tendsto rho l (𝓝 rhoM))
    (hw : Tendsto w l (𝓝 wM)) (hh : Tendsto h l (𝓝 hM)) (hr0 : rM ≠ 0) :
    Tendsto (fun x => coefficient alpha (r x) (rho x) (w x) (h x))
      l (𝓝 (coefficient alpha rM rhoM wM hM)) := by
  have hc := (((hh.sub_const 1).sub ((hw.pow 2).const_mul alpha)).sub
    (((hr.const_mul 2).add ((tendsto_const_nhds (x := 9*alpha)).div hr hr0)).mul hw)).add
    ((hρ.const_mul (6*alpha)).div (hr.pow 2) (pow_ne_zero 2 hr0))
  simpa [coefficient] using (hr.div_const (3*alpha)).mul hc

theorem regular_coefficient_zero_derivative_limit (h rho w : ℝ → ℝ)
    (alpha b rhoM wM : ℝ) (hb : b ≠ 0)
    (hρ : Tendsto rho (𝓝[<] b) (𝓝 rhoM))
    (hw : Tendsto w (𝓝[<] b) (𝓝 wM)) (hh : Tendsto h (𝓝[<] b) (𝓝 0))
    (hode : ∀ᶠ x in 𝓝[<] b,
      HasDerivAt h (-2*coefficient alpha x (rho x) (w x) (h x)*h x) x) :
    Tendsto (deriv h) (𝓝[<] b) (𝓝 0) := by
  have hr : Tendsto (fun x : ℝ => x) (𝓝[<] b) (𝓝 b) := tendsto_id.mono_left inf_le_left
  have hp := coefficient_tendsto (fun x : ℝ => x) rho w h alpha b rhoM wM 0 hr hρ hw hh hb
  have hz := (hp.const_mul (-2)).mul hh
  simp only [mul_zero] at hz
  apply hz.congr'
  filter_upwards [hode] with x hx
  exact hx.deriv.symm

theorem regular_coefficient_zero_ratio_limit (h rho w : ℝ → ℝ)
    (alpha b rhoM wM : ℝ) (hb : b ≠ 0)
    (hρ : Tendsto rho (𝓝[<] b) (𝓝 rhoM))
    (hw : Tendsto w (𝓝[<] b) (𝓝 wM)) (hh : Tendsto h (𝓝[<] b) (𝓝 0))
    (hode : ∀ᶠ x in 𝓝[<] b,
      HasDerivAt h (-2*coefficient alpha x (rho x) (w x) (h x)*h x) x) :
    Tendsto (fun x => h x/(b-x)) (𝓝[<] b) (𝓝 0) := by
  have hd := regular_coefficient_zero_derivative_limit h rho w alpha b rhoM wM hb hρ hw hh hode
  have hdf : ∀ᶠ x in 𝓝[<] b, HasDerivAt h (deriv h x) x := by
    filter_upwards [hode] with x hx
    rw [hx.deriv]
    exact hx
  have hgg : ∀ᶠ x in 𝓝[<] b, HasDerivAt (fun y => b-y) (-1) x := by
    exact Eventually.of_forall (fun x => by simpa using (hasDerivAt_const x b).sub (hasDerivAt_id x))
  have hg0 : Tendsto (fun x : ℝ => b-x) (𝓝[<] b) (𝓝 0) := by
    have hi : Tendsto (fun x : ℝ => x) (𝓝[<] b) (𝓝 b) := tendsto_id.mono_left inf_le_left
    simpa using (tendsto_const_nhds (x := b)).sub hi
  have hdiv : Tendsto (fun x => deriv h x/(-1)) (𝓝[<] b) (𝓝 0) := by
    simpa using hd.div_const (-1)
  exact HasDerivAt.lhopital_zero_nhdsLT hdf hgg
    (Eventually.of_forall (fun _ => by norm_num)) hh hg0 hdiv

theorem regular_coefficient_excludes_simple_zero (h rho w : ℝ → ℝ)
    (alpha b rhoM wM L : ℝ) (hb : b ≠ 0) (hL : L ≠ 0)
    (hρ : Tendsto rho (𝓝[<] b) (𝓝 rhoM))
    (hw : Tendsto w (𝓝[<] b) (𝓝 wM)) (hh : Tendsto h (𝓝[<] b) (𝓝 0))
    (hode : ∀ᶠ x in 𝓝[<] b,
      HasDerivAt h (-2*coefficient alpha x (rho x) (w x) (h x)*h x) x)
    (hlinear : Tendsto (fun x => h x/(b-x)) (𝓝[<] b) (𝓝 L)) : False := by
  have hz := regular_coefficient_zero_ratio_limit h rho w alpha b rhoM wM hb hρ hw hh hode
  exact hL (tendsto_nhds_unique hlinear hz)

end
end SWaveEndpoint
