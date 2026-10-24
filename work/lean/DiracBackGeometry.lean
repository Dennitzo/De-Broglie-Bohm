import Mathlib

/-! Differential consequences of the specified reduced Einstein equations. -/
namespace DiracBackGeometry
noncomputable section
open Filter Topology

def massIntegral (B b : ℝ) : ℝ := B*(b^2+1)/2
def constraint (A a B b E gam : ℝ) : ℝ := 2*B*b*a+A*(b^2+1)-gam*E
def kretschmann (A a aa B b bb : ℝ) : ℝ :=
  4*((aa/A)^2+2*(bb/B)^2+2*(a*b/(A*B))^2+((1+b^2)/B^2)^2)

theorem radial_first_integral_hasDerivAt (B b : ℝ → ℝ) (t : ℝ)
    (hB : HasDerivAt B (b t) t) (hBn : B t ≠ 0)
    (hb : HasDerivAt b (-((b t)^2+1)/(2*B t)) t) :
    HasDerivAt (fun x => massIntegral (B x) (b x)) 0 t := by
  convert (hB.mul ((hb.pow 2).add_const 1)).div_const 2 using 1 <;>
    dsimp [massIntegral] <;> field_simp <;> ring

theorem radial_first_integral_constant (B b : ℝ → ℝ) (l u t t0 : ℝ)
    (hB : ∀ x ∈ Set.Ioo l u, HasDerivAt B (b x) x)
    (hBn : ∀ x ∈ Set.Ioo l u, B x ≠ 0)
    (hb : ∀ x ∈ Set.Ioo l u, HasDerivAt b (-((b x)^2+1)/(2*B x)) x)
    (ht : t ∈ Set.Ioo l u) (ht0 : t0 ∈ Set.Ioo l u) :
    massIntegral (B t) (b t)=massIntegral (B t0) (b t0) := by
  have hd := fun x hx => radial_first_integral_hasDerivAt B b x (hB x hx) (hBn x hx) (hb x hx)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo l u).isPreconnected
    (fun x hx => (hd x hx).differentiableAt.differentiableWithinAt)
    (fun x hx => (hd x hx).deriv) ht ht0

theorem radial_acceleration_negative (B b : ℝ) (hB : 0 < B) :
    - (b^2+1)/(2*B) < 0 := by
  exact div_neg_of_neg_of_pos (by nlinarith [sq_nonneg b]) (by linarith)

theorem radial_velocity_from_integral (B b C : ℝ) (hB : B ≠ 0)
    (hC : massIntegral B b=C) : b^2=2*C/B-1 := by
  unfold massIntegral at hC
  field_simp
  nlinarith

theorem einstein_constraint_hasDerivAt (A a B b E : ℝ → ℝ) (t gam kap S1 : ℝ)
    (hA : HasDerivAt A (a t) t) (hB : HasDerivAt B (b t) t)
    (hBn : B t ≠ 0)
    (hb : HasDerivAt b (-((b t)^2+1)/(2*B t)) t)
    (ha : HasDerivAt a
      (A t*((b t)^2+1)/(2*(B t)^2)-a t*b t/B t-gam*kap*S1/(2*(B t)^3)) t)
    (hE : HasDerivAt E (-kap*S1*b t/(B t)^2) t) :
    HasDerivAt (fun x => constraint (A x) (a x) (B x) (b x) (E x) gam) 0 t := by
  convert ((((hB.mul hb).mul ha).const_mul 2).add (hA.mul ((hb.pow 2).add_const 1))).sub
    (hE.const_mul gam) using 1 <;> dsimp [constraint] <;>
      (try funext x) <;> (try dsimp) <;> field_simp <;> ring

theorem constraint_linear_equation (A a B b E gam : ℝ)
    (hB : B ≠ 0) (hb : b ≠ 0) :
    constraint A a B b E gam=0 ↔ a=(gam*E-A*(b^2+1))/(2*B*b) := by
  unfold constraint
  constructor <;> intro h <;> field_simp at * <;> nlinarith

theorem kretschmann_lower_bound (A a aa B b C : ℝ) (hB : B ≠ 0)
    (hb : b^2+1=2*C/B) :
    24*C^2/B^6 ≤ kretschmann A a aa B b (-C/B^2) := by
  have h1 : ((-C/B^2)/B)^2=C^2/B^6 := by field_simp <;> ring
  have h2 : ((1+b^2)/B^2)^2=4*C^2/B^6 := by
    rw [show 1+b^2=2*C/B by linarith]
    field_simp
    ring
  unfold kretschmann
  rw [h1,h2]
  have he : 24*C^2/B^6=24*(C^2/B^6) := by ring
  have he2 : 4*C^2/B^6=4*(C^2/B^6) := by ring
  rw [he,he2]
  nlinarith [sq_nonneg (aa/A),sq_nonneg (a*b/(A*B))]

theorem curvature_lower_bound_tendsto (C : ℝ) (hC : C ≠ 0) :
    Tendsto (fun B : ℝ => 24*C^2/B^6) (𝓝[>] (0 : ℝ)) atTop := by
  have hi : Tendsto (fun B : ℝ => (B⁻¹)^6) (𝓝[>] (0 : ℝ)) atTop :=
    (tendsto_pow_atTop (by norm_num : (6 : ℕ) ≠ 0)).comp tendsto_inv_nhdsGT_zero
  simpa [div_eq_mul_inv,inv_pow] using
    Tendsto.const_mul_atTop (by positivity : 0 < 24*C^2) hi

theorem conditional_curvature_diverges (K : ℝ → ℝ) (C : ℝ) (hC : C ≠ 0)
    (hK : ∀ᶠ B in 𝓝[>] (0 : ℝ), 24*C^2/B^6 ≤ K B) :
    Tendsto K (𝓝[>] (0 : ℝ)) atTop :=
  tendsto_atTop_mono' _ hK (curvature_lower_bound_tendsto C hC)

end
end DiracBackGeometry
