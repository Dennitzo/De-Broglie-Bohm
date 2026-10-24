import Mathlib

/-! Finite Gaussian guidance curvature and its covariance geometry. -/
namespace FoundGuidance
noncomputable section

def curvatureX (a b c x y : ℝ) : ℝ := -b*(b*x+c*y)
def curvatureY (a b c x y : ℝ) : ℝ := b*(a*x+b*y)

theorem curvature_trace (b : ℝ) : -b^2+b^2=0 := by ring

theorem curvature_determinant (a b c : ℝ) :
    (-b^2)*b^2-(-b*c)*(a*b)=b^2*(a*c-b^2) := by ring

theorem curvature_determinant_positive (a b c : ℝ)
    (hb : b ≠ 0) (hd : 0 < a*c-b^2) :
    0 < (-b^2)*b^2-(-b*c)*(a*b) := by
  rw [curvature_determinant]
  exact mul_pos (sq_pos_of_ne_zero hb) hd

theorem curvature_squared_x (a b c x y : ℝ) :
    curvatureX a b c (curvatureX a b c x y) (curvatureY a b c x y)
      = -b^2*(a*c-b^2)*x := by
  dsimp [curvatureX,curvatureY]
  ring

theorem curvature_squared_y (a b c x y : ℝ) :
    curvatureY a b c (curvatureX a b c x y) (curvatureY a b c x y)
      = -b^2*(a*c-b^2)*y := by
  dsimp [curvatureX,curvatureY]
  ring

theorem curvature_zero_iff (a b c x y : ℝ)
    (hb : b ≠ 0) (hd : 0 < a*c-b^2) :
    (curvatureX a b c x y=0 ∧ curvatureY a b c x y=0) ↔ (x=0 ∧ y=0) := by
  constructor
  · rintro ⟨hx,hy⟩
    have hdet : -b^2*(a*c-b^2) ≠ 0 := mul_ne_zero (neg_ne_zero.mpr (pow_ne_zero 2 hb)) (ne_of_gt hd)
    have hxx := curvature_squared_x a b c x y
    have hyy := curvature_squared_y a b c x y
    rw [hx,hy] at hxx hyy
    simp only [curvatureX,curvatureY,mul_zero,add_zero,neg_zero] at hxx hyy
    exact ⟨(mul_eq_zero.mp hxx.symm).resolve_left hdet,
      (mul_eq_zero.mp hyy.symm).resolve_left hdet⟩
  · rintro ⟨rfl,rfl⟩
    simp [curvatureX,curvatureY]

theorem curvature_nonzero_iff (a b c x y : ℝ)
    (hb : b ≠ 0) (hd : 0 < a*c-b^2) :
    (curvatureX a b c x y ≠ 0 ∨ curvatureY a b c x y ≠ 0) ↔ (x ≠ 0 ∨ y ≠ 0) := by
  simpa only [not_and_or,not_not] using not_congr (curvature_zero_iff a b c x y hb hd)

theorem curvature_density_tangent (a b c x y : ℝ) :
    (a*x+b*y)*curvatureX a b c x y +
    (b*x+c*y)*curvatureY a b c x y = 0 := by
  dsimp [curvatureX,curvatureY]
  ring

theorem curvature_covariance_11 (a b c : ℝ) :
    2*((-b^2)*(c/(2*(a*c-b^2)))+(-b*c)*(-b/(2*(a*c-b^2))))=0 := by ring

theorem curvature_covariance_12 (a b c : ℝ) :
    (-b^2)*(-b/(2*(a*c-b^2)))+(-b*c)*(a/(2*(a*c-b^2)))+
    (c/(2*(a*c-b^2)))*(a*b)+(-b/(2*(a*c-b^2)))*b^2=0 := by ring

theorem curvature_covariance_22 (a b c : ℝ) :
    2*((a*b)*(-b/(2*(a*c-b^2)))+b^2*(a/(2*(a*c-b^2))))=0 := by ring

theorem curvature_example (x y : ℝ) :
    curvatureX 2 1 2 x y=-(x+2*y) ∧ curvatureY 2 1 2 x y=2*x+y := by
  dsimp [curvatureX,curvatureY]
  constructor <;> ring

theorem curvature_example_norm (x y : ℝ) :
    (curvatureX 2 1 2 x y)^2+(curvatureY 2 1 2 x y)^2
      = x^2+y^2+4*(x+y)^2 := by
  dsimp [curvatureX,curvatureY]
  ring

theorem curvature_example_norm_lower (x y : ℝ) :
    x^2+y^2 ≤ (curvatureX 2 1 2 x y)^2+(curvatureY 2 1 2 x y)^2 := by
  rw [curvature_example_norm]
  nlinarith [sq_nonneg (x+y)]

theorem curvature_no_entanglement (a c x y : ℝ) :
    curvatureX a 0 c x y=0 ∧ curvatureY a 0 c x y=0 := by
  simp [curvatureX,curvatureY]

theorem multitime_denominator_positive (a b c s t : ℝ)
    (ha : 0<a) (hc : 0<c) (hd : 0<a*c-b^2) :
    0 < (1-(a*c-b^2)*s*t)^2+(a*s+c*t)^2 := by
  have hr := sq_nonneg (1-(a*c-b^2)*s*t)
  have hi := sq_nonneg (a*s+c*t)
  by_contra hn
  have hnonpos := le_of_not_gt hn
  have hreal : 1-(a*c-b^2)*s*t=0 := by nlinarith
  have himag : a*s+c*t=0 := by nlinarith
  have hproduct : 0 < (a*c-b^2)*(s*t) := by nlinarith [hreal]
  have hst : 0 < s*t := pos_of_mul_pos_right hproduct (le_of_lt hd)
  have hacst := mul_pos (mul_pos ha hc) hst
  have hm := congrArg (fun z : ℝ => z*(c*t)) himag
  dsimp only at hm
  nlinarith [sq_nonneg (c*t)]

theorem conditional_covariance_cross (r v w v' : ℝ) (hw : w ≠ 0) :
    r*v+((v'-r*v)/w)*w=v' := by
  field_simp
  ring

theorem conditional_covariance_diagonal (r u v w u' v' : ℝ) (hw : w ≠ 0)
    (hscale : r^2*(u-v^2/w)=u'-v'^2/w) :
    r^2*u+2*r*((v'-r*v)/w)*v+((v'-r*v)/w)^2*w=u' := by
  calc
    _ = r^2*(u-v^2/w)+v'^2/w := by field_simp; ring
    _ = u' := by rw [hscale]; ring

theorem conditional_map_determinant (r s : ℝ) : r*1-s*0=r := by ring

end
end FoundGuidance
