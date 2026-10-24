import Mathlib

/-! One-dimensional Green branches. Boundary conditions and state existence are inputs. -/
namespace FullRSETGreen
noncomputable section

def wronskian (y z dy dz : ℝ) : ℝ := y*dz-dy*z
def leftBranch (y : ℝ → ℝ) (zY W x : ℝ) : ℝ := -y x*zY/W
def rightBranch (z : ℝ → ℝ) (yY W x : ℝ) : ℝ := -yY*z x/W

theorem wronskian_hasDerivAt (y z dy dz : ℝ → ℝ) (x V : ℝ)
    (hy : HasDerivAt y (dy x) x) (hz : HasDerivAt z (dz x) x)
    (hdy : HasDerivAt dy (V*y x) x) (hdz : HasDerivAt dz (V*z x) x) :
    HasDerivAt (fun t => wronskian (y t) (z t) (dy t) (dz t)) 0 x := by
  convert (hy.mul hdz).sub (hdy.mul hz) using 1 <;>
    dsimp [wronskian] <;> (try funext t) <;> (try dsimp) <;> ring

theorem wronskian_constant_on_interval (y z dy dz V : ℝ → ℝ) (a b x x₀ : ℝ)
    (hy : ∀ t ∈ Set.Ioo a b, HasDerivAt y (dy t) t)
    (hz : ∀ t ∈ Set.Ioo a b, HasDerivAt z (dz t) t)
    (hdy : ∀ t ∈ Set.Ioo a b, HasDerivAt dy (V t*y t) t)
    (hdz : ∀ t ∈ Set.Ioo a b, HasDerivAt dz (V t*z t) t)
    (hx : x ∈ Set.Ioo a b) (hx₀ : x₀ ∈ Set.Ioo a b) :
    wronskian (y x) (z x) (dy x) (dz x) =
      wronskian (y x₀) (z x₀) (dy x₀) (dz x₀) := by
  have hw : ∀ t ∈ Set.Ioo a b,
      HasDerivAt (fun s => wronskian (y s) (z s) (dy s) (dz s)) 0 t := by
    intro t ht
    exact wronskian_hasDerivAt y z dy dz t (V t) (hy t ht) (hz t ht) (hdy t ht) (hdz t ht)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun t ht => (hw t ht).differentiableAt.differentiableWithinAt)
    (fun t ht => (hw t ht).deriv) hx hx₀

theorem change_left_solution_preserves_wronskian (y z dy dz c : ℝ) :
    wronskian (y+c*z) z (dy+c*dz) dz = wronskian y z dy dz := by
  unfold wronskian
  ring

theorem green_left_hasDerivAt (y : ℝ → ℝ) (x dy zY W : ℝ)
    (hy : HasDerivAt y dy x) :
    HasDerivAt (leftBranch y zY W) (-dy*zY/W) x := by
  exact (hy.neg.mul_const zY).div_const W

theorem green_right_hasDerivAt (z : ℝ → ℝ) (x dz yY W : ℝ)
    (hz : HasDerivAt z dz x) :
    HasDerivAt (rightBranch z yY W) (-yY*dz/W) x := by
  exact (hz.const_mul (-yY)).div_const W

theorem green_branches_match (y z : ℝ → ℝ) (x W : ℝ) :
    leftBranch y (z x) W x = rightBranch z (y x) W x := rfl

theorem green_derivative_jump (y z dy dz : ℝ)
    (hW : wronskian y z dy dz ≠ 0) :
    (-y*dz/wronskian y z dy dz)-(-dy*z/wronskian y z dy dz) = -1 := by
  field_simp
  unfold wronskian
  ring

theorem green_left_homogeneous_equation (y dy : ℝ → ℝ) (x V zY W : ℝ)
    (hy : HasDerivAt dy (V*y x) x) :
    HasDerivAt (fun t => -dy t*zY/W) (V*leftBranch y zY W x) x := by
  convert (hy.neg.mul_const zY).div_const W using 1
  unfold leftBranch
  ring

theorem green_boundary_change (yX yY zX zY W c : ℝ) :
    -(yX+c*zX)*zY/W-(-yX*zY/W) = -c*zX*zY/W ∧
    -(yY+c*zY)*zX/W-(-yY*zX/W) = -c*zX*zY/W := by
  constructor <;> ring

theorem green_boundary_change_nonzero (zX zY W c : ℝ)
    (hX : zX ≠ 0) (hY : zY ≠ 0) (hW : W ≠ 0) (hc : c ≠ 0) :
    -c*zX*zY/W ≠ 0 := by
  exact div_ne_zero (mul_ne_zero (mul_ne_zero (neg_ne_zero.mpr hc) hX) hY) hW

theorem boundary_difference_derivative (z : ℝ → ℝ) (x dz c zY W : ℝ)
    (hz : HasDerivAt z dz x) :
    HasDerivAt (fun t => -c*z t*zY/W) (-c*dz*zY/W) x := by
  exact ((hz.const_mul (-c)).mul_const zY).div_const W

def robinCorrection (zX zY zeta gamma : ℝ) : ℝ := zX*zY/(zeta+gamma)

theorem robin_left_initial_data (y q : ℝ → ℝ) (x zeta gamma : ℝ)
    (hy : HasDerivAt y 1 x) (hq : HasDerivAt q (-zeta) x)
    (hy₀ : y x=0) (hq₀ : q x=1) :
    (q x+(zeta+gamma)*y x=1) ∧
    HasDerivAt (fun t => q t+(zeta+gamma)*y t) gamma x := by
  constructor
  · rw [hy₀,hq₀]; ring
  · convert hq.add (hy.const_mul (zeta+gamma)) using 1 <;> ring

theorem robin_initial_wronskian (zeta gamma : ℝ) :
    wronskian 0 1 1 (-zeta) = -1 ∧
    wronskian 1 1 gamma (-zeta) = -(zeta+gamma) := by
  unfold wronskian
  constructor <;> ring

theorem robin_minus_dirichlet_green (yX qX qY zeta gamma : ℝ)
    (h : zeta+gamma ≠ 0) :
    -(qX+(zeta+gamma)*yX)*qY/(-(zeta+gamma)) - (-yX*qY/(-1)) =
      robinCorrection qX qY zeta gamma := by
  unfold robinCorrection
  field_simp
  ring

theorem robin_diagonal_positive (z zeta gamma : ℝ)
    (hz : z ≠ 0) (hzet : 0 < zeta) (hgam : 0 ≤ gamma) :
    0 < robinCorrection z z zeta gamma := by
  unfold robinCorrection
  exact div_pos (mul_self_pos.mpr hz) (by positivity)

theorem robin_parameter_resolvent_difference (zX zY zeta gamma₁ gamma₂ : ℝ)
    (h₁ : zeta+gamma₁ ≠ 0) (h₂ : zeta+gamma₂ ≠ 0) :
    robinCorrection zX zY zeta gamma₁-robinCorrection zX zY zeta gamma₂ =
      (gamma₂-gamma₁)*zX*zY/((zeta+gamma₁)*(zeta+gamma₂)) := by
  unfold robinCorrection
  field_simp
  ring

theorem robin_diagonal_antitone (z zeta gamma₁ gamma₂ : ℝ)
    (hzet : 0 < zeta) (h₁ : 0 ≤ gamma₁) (h₁₂ : gamma₁ ≤ gamma₂) :
    robinCorrection z z zeta gamma₂ ≤ robinCorrection z z zeta gamma₁ := by
  unfold robinCorrection
  apply div_le_div_of_nonneg_left (mul_self_nonneg z) (by positivity)
  linarith

theorem robin_parameter_hasDerivAt (zX zY zeta gamma : ℝ)
    (h : zeta+gamma ≠ 0) :
    HasDerivAt (fun g => robinCorrection zX zY zeta g)
      (-zX*zY/(zeta+gamma)^2) gamma := by
  have hd := (hasDerivAt_const gamma zeta).add (hasDerivAt_id gamma)
  convert (hasDerivAt_const gamma (zX*zY)).div hd h using 1 <;>
    dsimp [robinCorrection] <;> (try funext g) <;> (try dsimp) <;> ring

theorem robin_strict_diagonal_dependence (z zeta gamma₁ gamma₂ : ℝ)
    (hz : z ≠ 0) (hzet : 0 < zeta) (h₁ : 0 ≤ gamma₁) (h₁₂ : gamma₁ < gamma₂) :
    robinCorrection z z zeta gamma₂ < robinCorrection z z zeta gamma₁ := by
  unfold robinCorrection
  apply div_lt_div_of_pos_left (mul_self_pos.mpr hz) (by positivity)
  linarith

end
end FullRSETGreen
