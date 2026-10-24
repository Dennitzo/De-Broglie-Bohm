import Mathlib

/-! Fixed-point consequences of the explicitly specified invariant Dirac
   source. The digamma source formula is an external analytic input. -/
namespace AbsoluteSourceDeSitter
noncomputable section

def density (H value : ℝ) : ℝ := H^4*value/(960*Real.pi^2)
def conformalDensity (H : ℝ) : ℝ := density H 11

theorem conformal_density_nonnegative (H : ℝ) : 0 ≤ conformalDensity H := by
  unfold conformalDensity density
  positivity

theorem invariant_trace (rho : ℝ) : -rho+3*(-rho)= -4*rho := by ring

theorem conformal_trace (H : ℝ) : -4*conformalDensity H= -11*H^4/(240*Real.pi^2) := by
  unfold conformalDensity density
  ring

theorem fixedpoint_equation (F : ℝ → ℝ) (N G m x : ℝ) (hx : x ≠ 0)
    (hroot : (N*G*m^2)*F x=360*Real.pi*x^2) :
    3*(m/x)^2=8*Real.pi*G*N*density (m/x) (F x) := by
  have hpi := Real.pi_ne_zero
  unfold density
  calc
    3*(m/x)^2 = m^2/(120*Real.pi*x^4)*(360*Real.pi*x^2) := by field_simp; ring
    _ = m^2/(120*Real.pi*x^4)*((N*G*m^2)*F x) := by rw [hroot]
    _ = 8*Real.pi*G*N*((m/x)^4*F x/(960*Real.pi^2)) := by field_simp; ring

theorem fixedpoint_root_exists (F : ℝ → ℝ) (lam : ℝ) (hlam : 0 < lam)
    (hF : ContinuousOn F (Set.Icc 0 1)) (hzero : F 0=11) (hone : F 1<0) :
    ∃ x : ℝ, 0 < x ∧ x < 1 ∧ lam*F x=360*Real.pi*x^2 := by
  let q : ℝ → ℝ := fun x => lam*F x-360*Real.pi*x^2
  have hq : ContinuousOn q (Set.Icc 0 1) :=
    (continuousOn_const.mul hF).sub (continuousOn_const.mul (continuousOn_id.pow 2))
  have hq0 : 0 < q 0 := by dsimp [q]; rw [hzero]; nlinarith
  have hq1 : q 1 < 0 := by
    have hh := mul_neg_of_pos_of_neg hlam hone
    dsimp [q]
    nlinarith [Real.pi_pos]
  obtain ⟨x,hx,hroot⟩ := intermediate_value_Icc' (by norm_num : (0 : ℝ) ≤ 1) hq
    (show (0 : ℝ) ∈ Set.Icc (q 1) (q 0) from ⟨le_of_lt hq1,le_of_lt hq0⟩)
  refine ⟨x,?_,?_,?_⟩
  · have hn : x ≠ 0 := by intro he; subst x; linarith
    exact lt_of_le_of_ne hx.1 (Ne.symm hn)
  · have hn : x ≠ 1 := by intro he; subst x; linarith
    exact lt_of_le_of_ne hx.2 hn
  · dsimp [q] at hroot
    linarith

theorem fixedpoint_hubble_exceeds_mass (m x : ℝ) (hm : 0 < m)
    (hx : 0 < x) (hx1 : x < 1) : m < m/x := by
  apply (lt_div_iff₀ hx).2
  nlinarith

theorem conformal_nonzero_branch (N G H : ℝ) (hN : N ≠ 0) (hG : G ≠ 0)
    (hH : H ≠ 0) :
    3*H^2=8*Real.pi*G*N*conformalDensity H ↔ H^2=360*Real.pi/(11*N*G) := by
  have hpi := Real.pi_ne_zero
  unfold conformalDensity density
  field_simp
  constructor <;> intro h <;> nlinarith [sq_pos_of_ne_zero hH]

theorem invariant_source_counterterms_cancel_trace_difference
    (rho1 rho0 localTerm : ℝ) :
    (-4*rho1+localTerm)-(-4*rho0+localTerm)= -4*(rho1-rho0) := by ring

def psiPartial20 : ℝ := ∑ n ∈ Finset.range 20,
  1/(((n : ℝ)+1)*(((n : ℝ)+1)^2+1))

theorem euler_constant_bound_certificate : (751 : ℝ)/1296 < 29/50 := by norm_num

theorem log_two_lower_certificate : (61 : ℝ)/48-56/81=751/1296 := by norm_num

theorem finite_digamma_lower_bound : (21 : ℝ)+240*(29/50-psiPartial20)<0 := by
  norm_num [psiPartial20,Finset.sum_range_succ]

theorem endpoint_source_negative (gamma psiOne : ℝ)
    (hgamma : gamma ≤ 29/50) (hpsi : psiPartial20-gamma ≤ psiOne) :
    21-240*psiOne<0 := by
  have hpartial := finite_digamma_lower_bound
  linarith

end
end AbsoluteSourceDeSitter
