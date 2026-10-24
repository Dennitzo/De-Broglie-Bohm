import Mathlib

/-! Exact consequences of the stated conformal Dirac trace anomaly. -/
namespace AbsoluteSourceTrace
noncomputable section
open Filter Topology

def schwarzAnomaly (M r : ℝ) : ℝ := 7*M^2/(120*Real.pi^2*r^6)

theorem schwarz_anomaly_positive (M r : ℝ) (hM : M ≠ 0) (hr : 0 < r) :
    0 < schwarzAnomaly M r := by
  unfold schwarzAnomaly
  positivity

theorem schwarz_trace_obstruction (M r gamma : ℝ)
    (hM : M ≠ 0) (hr : 0 < r) (hg : 0 < gamma) :
    (0 : ℝ) ≠ gamma*schwarzAnomaly M r := by
  have hp := mul_pos hg (schwarz_anomaly_positive M r hM hr)
  linarith

theorem schwarz_trace_scaled (M r : ℝ) (hr : r ≠ 0) :
    r^6*schwarzAnomaly M r=7*M^2/(120*Real.pi^2) := by
  unfold schwarzAnomaly
  field_simp

theorem schwarz_anomaly_tendsto (M : ℝ) (hM : M ≠ 0) :
    Tendsto (schwarzAnomaly M) (𝓝[>] (0 : ℝ)) atTop := by
  have hi : Tendsto (fun r : ℝ => (r⁻¹)^6) (𝓝[>] (0 : ℝ)) atTop :=
    (tendsto_pow_atTop (by norm_num : (6 : ℕ) ≠ 0)).comp tendsto_inv_nhdsGT_zero
  have hc : 0 < 7*M^2/(120*Real.pi^2) := by positivity
  convert Tendsto.const_mul_atTop hc hi using 1
  ext r
  simp [schwarzAnomaly,div_eq_mul_inv,inv_pow]
  ring

def regularRicciScaled (M a R : ℝ) : ℝ := -2*a^2+6*M*a^2/R

theorem regular_ricci_scaled_limit (M a : ℝ) :
    Tendsto (regularRicciScaled M a) atTop (𝓝 (-2*a^2)) := by
  have hi : Tendsto (fun R : ℝ => R⁻¹) atTop (𝓝 0) := tendsto_inv_atTop_zero
  change Tendsto (fun R : ℝ => -2*a^2+6*M*a^2/R) atTop (𝓝 (-2*a^2))
  simpa [div_eq_mul_inv,neg_mul] using
    (tendsto_const_nhds (x := -(2*a^2))).add (hi.const_mul (6*M*a^2))

/- A trace solution would have a vanishing residual on every sufficiently
   large radius. This lemma states the general asymptotic contradiction. -/
theorem nonzero_limit_excludes_eventual_zero (f : ℝ → ℝ) (L : ℝ)
    (hf : Tendsto f atTop (𝓝 L)) (hL : L ≠ 0) :
    ¬ (∀ᶠ R in atTop, f R=0) := by
  intro hz
  have he : (fun _ : ℝ => (0 : ℝ)) =ᶠ[atTop] f := by
    filter_upwards [hz] with R hR using hR.symm
  have hh : Tendsto f atTop (𝓝 0) := tendsto_const_nhds.congr' he
  exact hL (tendsto_nhds_unique hf hh)

/-- Four-dimensional anomaly numerator in the areal inverse radius. -/
def anomalyPoly (M a u : ℝ) : ℝ :=
  336*M^2*u^6+224*M*a^2*u^7+(96*a^4-1008*M^2*a^2)*u^8
    -424*M*a^4*u^9+888*M^2*a^4*u^10

/-- Covariant wave operator applied to the Ricci scalar. -/
def boxRicciPoly (M a u : ℝ) : ℝ :=
  -24*a^2*u^6+184*M*a^2*u^7+(32*a^4-300*M^2*a^2)*u^8
    -230*M*a^4*u^9+360*M^2*a^4*u^10

def regularTrace (M a G b R : ℝ) : ℝ :=
  -2*a^2/R^4+6*M*a^2/R^5+
    G/(720*Real.pi)*(anomalyPoly M a R⁻¹+b*boxRicciPoly M a R⁻¹)

def scaledPoly (M a G b u : ℝ) : ℝ :=
  -2*a^2+6*M*a^2*u+G/(720*Real.pi)*
    (336*M^2*u^2+224*M*a^2*u^3+(96*a^4-1008*M^2*a^2)*u^4
      -424*M*a^4*u^5+888*M^2*a^4*u^6
      +b*(-24*a^2*u^2+184*M*a^2*u^3+(32*a^4-300*M^2*a^2)*u^4
        -230*M*a^4*u^5+360*M^2*a^4*u^6))

theorem regular_trace_scaled_identity (M a G b R : ℝ) (hR : R ≠ 0) :
    R^4*regularTrace M a G b R=scaledPoly M a G b R⁻¹ := by
  unfold regularTrace anomalyPoly boxRicciPoly scaledPoly
  field_simp

theorem scaled_polynomial_limit (M a G b : ℝ) :
    Tendsto (fun R : ℝ => scaledPoly M a G b R⁻¹) atTop (𝓝 (-2*a^2)) := by
  have hc : Continuous (scaledPoly M a G b) := by
    unfold scaledPoly
    fun_prop
  have hi : Tendsto (fun R : ℝ => R⁻¹) atTop (𝓝 0) := tendsto_inv_atTop_zero
  simpa [Function.comp_def,scaledPoly] using (hc.tendsto 0).comp hi

theorem regular_trace_scaled_limit (M a G b : ℝ) :
    Tendsto (fun R : ℝ => R^4*regularTrace M a G b R) atTop (𝓝 (-2*a^2)) := by
  apply (scaled_polynomial_limit M a G b).congr'
  filter_upwards [eventually_gt_atTop (0 : ℝ)] with R hR
  exact (regular_trace_scaled_identity M a G b R (ne_of_gt hR)).symm

theorem regular_trace_no_exact_tail (M a G b R0 : ℝ) (ha : a ≠ 0) :
    ¬ (∀ R : ℝ, R0 < R → regularTrace M a G b R=0) := by
  intro hsol
  have hlim := regular_trace_scaled_limit M a G b
  have hnon : -2*a^2 ≠ 0 := mul_ne_zero (by norm_num) (pow_ne_zero 2 ha)
  apply nonzero_limit_excludes_eventual_zero _ _ hlim hnon
  filter_upwards [eventually_gt_atTop R0] with R hR
  rw [hsol R hR,mul_zero]

theorem schwarz_reduction (M G b R : ℝ) :
    regularTrace M 0 G b R=7*G*M^2/(15*Real.pi*R^6) := by
  unfold regularTrace anomalyPoly boxRicciPoly
  simp [div_eq_mul_inv,inv_pow]
  ring

theorem schwarz_reduced_trace_positive (M G b R : ℝ)
    (hM : M ≠ 0) (hG : 0 < G) (hR : 0 < R) :
    0 < regularTrace M 0 G b R := by
  rw [schwarz_reduction]
  positivity

end
end AbsoluteSourceTrace
