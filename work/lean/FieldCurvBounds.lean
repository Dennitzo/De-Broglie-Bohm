import Mathlib

/-! Pointwise bilinear bounds and actual radial coupling derivatives. -/
namespace FieldCurvBounds
noncomputable section
def q (a b c d : ℝ) : ℝ := a^2+b^2+c^2+d^2
def overlapR (a b c d : ℝ) : ℝ := a*c+b*d
def overlapI (a b c d : ℝ) : ℝ := a*d-b*c
def potential (mu w a b c d : ℝ) : ℝ := 2*mu*overlapR a b c d-2*w*overlapI a b c d
def radius (a r : ℝ) : ℝ := Real.sqrt (r^2+a^2)

theorem potential_norm_identity (mu w a b c d : ℝ) :
    q (mu*c-w*d) (mu*d+w*c) (mu*a+w*b) (mu*b-w*a)=
      (mu^2+w^2)*q a b c d := by
  unfold q
  ring

theorem potential_gap_identity (mu w a b c d : ℝ) :
    (mu^2+w^2)*(q a b c d)^2-(potential mu w a b c d)^2=
      (mu^2+w^2)*(a^2+b^2-c^2-d^2)^2+
      4*(mu*overlapI a b c d+w*overlapR a b c d)^2 := by
  unfold q potential overlapR overlapI
  ring

theorem potential_abs_bound (mu w a b c d Omega : ℝ)
    (hO : 0 ≤ Omega) (hc : mu^2+w^2 ≤ Omega^2) :
    |potential mu w a b c d| ≤ Omega*q a b c d := by
  have hg := potential_gap_identity mu w a b c d
  have hq : 0 ≤ q a b c d := by unfold q; positivity
  have he : (potential mu w a b c d)^2 ≤ (mu^2+w^2)*(q a b c d)^2 := by
    have hnn : 0 ≤ (mu^2+w^2)*(a^2+b^2-c^2-d^2)^2+
        4*(mu*overlapI a b c d+w*overlapR a b c d)^2 := by positivity
    linarith
  have hc' := mul_le_mul_of_nonneg_right hc (sq_nonneg (q a b c d))
  apply (sq_le_sq₀ (abs_nonneg _) (mul_nonneg hO hq)).1
  simp only [sq_abs]
  nlinarith

theorem four_cauchy_identity (a b c d e f g h : ℝ) :
    q a b c d*q e f g h-(a*e+b*f+c*g+d*h)^2=
    (a*f-b*e)^2+(a*g-c*e)^2+(a*h-d*e)^2+
      (b*g-c*f)^2+(b*h-d*f)^2+(c*h-d*g)^2 := by
  unfold q
  ring

theorem four_cauchy_bound (a b c d e f g h A B : ℝ)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hu : q a b c d ≤ A^2)
    (hx : q e f g h ≤ B^2) : |a*e+b*f+c*g+d*h| ≤ A*B := by
  have hq : 0 ≤ q a b c d := by unfold q; positivity
  have hqx : 0 ≤ q e f g h := by unfold q; positivity
  have hg := four_cauchy_identity a b c d e f g h
  have hsq : (a*e+b*f+c*g+d*h)^2 ≤ q a b c d*q e f g h := by
    have hnn : 0 ≤ (a*f-b*e)^2+(a*g-c*e)^2+(a*h-d*e)^2+
      (b*g-c*f)^2+(b*h-d*f)^2+(c*h-d*g)^2 := by positivity
    linarith
  have hm := mul_le_mul hu hx hqx (sq_nonneg A)
  apply (sq_le_sq₀ (abs_nonneg _) (mul_nonneg hA hB)).1
  simp only [sq_abs]
  nlinarith

theorem radial_kinetic_abs_bound (a b c d ar br cr dr A B : ℝ)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hu : q a b c d ≤ A^2)
    (hx : q ar br cr dr ≤ B^2) : |a*br-b*ar-c*dr+d*cr| ≤ A*B := by
  have hh := four_cauchy_bound a b c d br (-ar) (-dr) cr A B hA hB hu
    (show q br (-ar) (-dr) cr ≤ B^2 by convert hx using 1 <;> unfold q <;> ring)
  convert hh using 1 <;> congr 1 <;> ring

theorem momentum_abs_bound (a b c d ar br cr dr A B : ℝ)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hu : q a b c d ≤ A^2)
    (hx : q ar br cr dr ≤ B^2) : |a*br-b*ar+c*dr-d*cr| ≤ A*B := by
  have hh := four_cauchy_bound a b c d br (-ar) dr (-cr) A B hA hB hu
    (show q br (-ar) dr (-cr) ≤ B^2 by convert hx using 1 <;> unfold q <;> ring)
  convert hh using 1 <;> congr 1 <;> ring

theorem overlap_twice_abs_bound (a b c d : ℝ) : |2*overlapI a b c d| ≤ q a b c d := by
  have hh := potential_abs_bound 0 1 a b c d 1 (by norm_num) (by norm_num)
  simpa [potential,abs_mul,abs_neg] using hh

theorem radius_lower (a r : ℝ) (ha : 0 < a) : a ≤ radius a r := by
  have hs := Real.sq_sqrt (show 0 ≤ r^2+a^2 by positivity)
  have hp := Real.sqrt_nonneg (r^2+a^2)
  unfold radius
  nlinarith [sq_nonneg r]

theorem radius_hasDerivAt (a r : ℝ) (ha : 0 < a) :
    HasDerivAt (radius a) (r/radius a r) r := by
  convert (((hasDerivAt_id r).pow 2).add_const (a^2)).sqrt
    (show r^2+a^2 ≠ 0 by positivity) using 1 <;>
    dsimp [radius] <;> field_simp <;> ring

theorem regular_coupling_hasDerivAt (a k r : ℝ) (ha : 0 < a) :
    HasDerivAt (fun x => k/radius a x) (-k*r/(radius a r)^3) r := by
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  convert (hasDerivAt_const r k).div (radius_hasDerivAt a r ha) (ne_of_gt hR) using 1
  dsimp
  field_simp
  ring

theorem regular_coupling_sq_bound (a k r : ℝ) (ha : 0 < a) :
    (k/radius a r)^2 ≤ k^2/a^2 := by
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  rw [div_pow]
  apply div_le_div_of_nonneg_left (sq_nonneg k) (sq_pos_of_pos ha)
  exact pow_le_pow_left₀ (le_of_lt ha) (radius_lower a r ha) 2

theorem singular_coupling_hasDerivAt (k r : ℝ) (hr : r ≠ 0) :
    HasDerivAt (fun x : ℝ => k/x) (-k/r^2) r := by
  convert (hasDerivAt_const r k).div (hasDerivAt_id r) hr using 1 <;> simp <;> ring

theorem singular_coupling_arbitrarily_large (k L : ℝ) (hk : 0 < k) (hL : 0 < L) :
    ∃ r : ℝ, 0 < r ∧ k/r > L := by
  refine ⟨k/(2*L),by positivity,?_⟩
  have hk0 : k ≠ 0 := ne_of_gt hk
  have hL0 : L ≠ 0 := ne_of_gt hL
  have he : k/(k/(2*L))=2*L := by field_simp
  rw [he]
  linarith

theorem singular_coupling_tendsto (k : ℝ) (hk : 0 < k) :
    Filter.Tendsto (fun r : ℝ => k/r) (nhdsWithin 0 (Set.Ioi 0)) Filter.atTop := by
  simpa only [div_eq_mul_inv] using
    (Filter.Tendsto.const_mul_atTop hk (tendsto_inv_nhdsGT_zero :
      Filter.Tendsto (fun r : ℝ => r⁻¹) (nhdsWithin 0 (Set.Ioi 0)) Filter.atTop))

theorem regular_coupling_gradient_bound (a k r : ℝ) (ha : 0 < a) :
    |(-k*r)/(radius a r)^3| ≤ |k|/a^2 := by
  have hR : 0 < radius a r := lt_of_lt_of_le ha (radius_lower a r ha)
  have hsq : (radius a r)^2=r^2+a^2 := Real.sq_sqrt (by positivity)
  have habs : |r| ≤ radius a r := by
    apply (sq_le_sq₀ (abs_nonneg _) (le_of_lt hR)).1
    simp only [sq_abs]
    nlinarith [sq_nonneg a]
  rw [abs_div,abs_mul,abs_neg,abs_of_pos (pow_pos hR 3)]
  calc
    |k| * |r|/(radius a r)^3 ≤ |k| * radius a r/(radius a r)^3 := by gcongr
    _ = |k|/(radius a r)^2 := by field_simp
    _ ≤ |k|/a^2 := by
      apply div_le_div_of_nonneg_left (abs_nonneg k) (sq_pos_of_pos ha)
      exact pow_le_pow_left₀ (le_of_lt ha) (radius_lower a r ha) 2

theorem normal_energy_stress_bound (P W R a A B Omega : ℝ)
    (ha : 0 < a) (hR : a ≤ R) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hO : 0 ≤ Omega) (hP : |P| ≤ A*B) (hW : |W| ≤ Omega*A^2) :
    |(P+W)/(4*Real.pi*R^2)| ≤ (A*B+Omega*A^2)/(4*Real.pi*a^2) := by
  have hRp : 0 < R := lt_of_lt_of_le ha hR
  rw [abs_div,abs_of_pos (by positivity : 0 < 4*Real.pi*R^2)]
  calc
    |P+W|/(4*Real.pi*R^2) ≤ (A*B+Omega*A^2)/(4*Real.pi*R^2) := by
      gcongr
      exact le_trans (abs_add_le P W) (add_le_add hP hW)
    _ ≤ (A*B+Omega*A^2)/(4*Real.pi*a^2) := by
      apply div_le_div_of_nonneg_left (by positivity) (by positivity)
      gcongr

theorem radial_stress_bound (P R a A B : ℝ)
    (ha : 0 < a) (hR : a ≤ R) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hP : |P| ≤ A*B) : |P/(4*Real.pi*R^2)| ≤ A*B/(4*Real.pi*a^2) := by
  have hh := normal_energy_stress_bound P 0 R a A B 0 ha hR hA hB (by norm_num) hP (by simp)
  simpa using hh

theorem angular_stress_bound (k s R a A : ℝ)
    (ha : 0 < a) (hR : a ≤ R) (hs : |s| ≤ A^2) :
    |(-k*s)/(8*Real.pi*R^3)| ≤ |k| *A^2/(8*Real.pi*a^3) := by
  have hRp : 0 < R := lt_of_lt_of_le ha hR
  rw [abs_div,abs_mul,abs_neg,abs_of_pos (by positivity : 0 < 8*Real.pi*R^3)]
  calc
    |k| *|s|/(8*Real.pi*R^3) ≤ |k| *A^2/(8*Real.pi*R^3) := by gcongr
    _ ≤ |k| *A^2/(8*Real.pi*a^3) := by
      apply div_le_div_of_nonneg_left (by positivity) (by positivity)
      gcongr

end
end FieldCurvBounds
