import Mathlib

/-! Componentwise covariance of the semiclassical equation under finite
   local changes of the renormalized stress prescription. -/
namespace AbsoluteSourceRenormalization
noncomputable section

def residual (z lam alpha beta gamma G g I J T : ℝ) : ℝ :=
  z*G+lam*g+alpha*I+beta*J-gamma*T

def finiteShift (m c0 c1 c2 c3 G g I J : ℝ) : ℝ :=
  c0*m^4*g+c1*m^2*G+c2*I+c3*J

theorem matching_coupling_shifts
    (z lam alpha beta gamma m c0 c1 c2 c3 G g I J T : ℝ) :
    residual (z+gamma*c1*m^2) (lam+gamma*c0*m^4)
      (alpha+gamma*c2) (beta+gamma*c3) gamma G g I J
      (T+finiteShift m c0 c1 c2 c3 G g I J)=
    residual z lam alpha beta gamma G g I J T := by
  unfold residual finiteShift
  ring

theorem equivalent_field_equations
    (z lam alpha beta gamma m c0 c1 c2 c3 G g I J T : ℝ) :
    (residual (z+gamma*c1*m^2) (lam+gamma*c0*m^4)
      (alpha+gamma*c2) (beta+gamma*c3) gamma G g I J
      (T+finiteShift m c0 c1 c2 c3 G g I J)=0) ↔
    residual z lam alpha beta gamma G g I J T=0 := by
  rw [matching_coupling_shifts]

theorem tensor_component_covariance {ι : Type*}
    (z lam alpha beta gamma m c0 c1 c2 c3 : ℝ)
    (G g I J T : ι → ℝ) :
    (fun i => residual (z+gamma*c1*m^2) (lam+gamma*c0*m^4)
      (alpha+gamma*c2) (beta+gamma*c3) gamma (G i) (g i) (I i) (J i)
      (T i+finiteShift m c0 c1 c2 c3 (G i) (g i) (I i) (J i))) =
    (fun i => residual z lam alpha beta gamma (G i) (g i) (I i) (J i) (T i)) := by
  funext i
  exact matching_coupling_shifts z lam alpha beta gamma m c0 c1 c2 c3
    (G i) (g i) (I i) (J i) (T i)

theorem state_difference_invariant (T1 T0 S : ℝ) :
    (T1+S)-(T0+S)=T1-T0 := by ring

theorem fixed_geometry_counterterm_mismatch
    (z lam alpha beta gamma G g I J T S : ℝ) :
    residual z lam alpha beta gamma G g I J (T+S)=
      residual z lam alpha beta gamma G g I J T-gamma*S := by
  unfold residual
  ring

end
end AbsoluteSourceRenormalization
