import Mathlib

/-! Algebraic interior barriers and conditional curvature growth.
Invariance and global existence of a solution are separate analytic inputs. -/
namespace SWaveInterior
noncomputable section
open Filter
open scoped Topology

def F (s z k b : ℝ) : ℝ :=
  z*(1-2*k-b)/3+k^2/3-4*k-6*k/z+2*s*k^2/z+6*s*k^2/z^2

def J (s z k b : ℝ) : ℝ :=
  z*(1-k)*(b-1+2*k)/3+k^3/3-13*k^2/3+3*k-5*k^2/z+
    2*s*k^2*(k-1)/z+6*s*k^3/z^2

theorem interior_decay_bound (s z k b : ℝ) (hz : 0 < z)
    (hs : s ≤ z) (hk : (3/4 : ℝ) ≤ k) (hk1 : k ≤ 1) (hb : 0 ≤ b) :
    F s z k b ≤ -z/6-5/4 := by
  have hk0 : 0 ≤ k := by linarith
  have hk2 : k^2 ≤ k := by nlinarith [mul_nonneg hk0 (sub_nonneg.mpr hk1)]
  have hsz : s/z ≤ 1 := (div_le_one hz).2 hs
  have ha : (s/z)*k^2 ≤ k^2 := mul_le_of_le_one_left (sq_nonneg k) hsz
  have hb1 : 0 ≤ z*(k-3/4) := mul_nonneg (le_of_lt hz) (sub_nonneg.mpr hk)
  have hb2 : 0 ≤ z*b := mul_nonneg (le_of_lt hz) hb
  have hc : 0 ≤ (k-(s/z)*k^2)/z := div_nonneg (by linarith) (le_of_lt hz)
  have hi : F s z k b+z/6+5/4 =
      -(2/3 : ℝ)*z*(k-3/4)-z*b/3+(7/3 : ℝ)*(k^2-k)-
      (5/3 : ℝ)*(k-3/4)+2*((s/z)*k^2-k^2)-6*((k-(s/z)*k^2)/z) := by
    unfold F
    field_simp
    <;> ring
  linarith

theorem lower_barrier (s z b : ℝ) (hs0 : 0 ≤ s) (hs : s ≤ z)
    (hz : 16 ≤ z) (hb : 0 ≤ b) :
    (125/768 : ℝ) ≤ J s z (3/4) b := by
  have hz0 : 0 < z := by linarith
  have hsz : s/z ≤ 1 := (div_le_one hz0).2 hs
  have hzp : 0 ≤ z*b := mul_nonneg (le_of_lt hz0) hb
  have hsp : 0 ≤ s/z^2 := div_nonneg hs0 (sq_nonneg z)
  have hinv : 1/z ≤ (1/16 : ℝ) := (div_le_div_iff₀ hz0 (by norm_num)).2 (by linarith)
  have hi : J s z (3/4) b = z/24+z*b/12-3/64-45/(16*z)-9*s/(32*z)+81*s/(32*z^2) := by
    unfold J
    ring
  rw [hi]
  have ht1 : s/z ≤ 1 := hsz
  have ht2 : 0 ≤ s/z^2 := hsp
  have he : z/24+z*b/12-3/64-45/(16*z)-9*s/(32*z)+81*s/(32*z^2) =
      z/24+z*b/12-3/64-(45/16 : ℝ)*(1/z)-(9/32 : ℝ)*(s/z)+(81/32 : ℝ)*(s/z^2) := by
    ring
  rw [he]
  linarith

theorem upper_barrier (s z b : ℝ) (hz : 0 < z) (hs : s ≤ z) :
    J s z 1 b ≤ -1+1/z := by
  have hsz : s/z^2 ≤ 1/z := by
    apply (div_le_div_iff₀ (sq_pos_of_pos hz) hz).2
    nlinarith [mul_le_mul_of_nonneg_left hs (le_of_lt hz)]
  have hi : J s z 1 b = -1-5/z+6*(s/z^2) := by
    unfold J
    ring
  rw [hi, show 5/z=5*(1/z) by ring]
  linarith

theorem radius_gap_hasDerivAt (z : ℝ → ℝ) (s k : ℝ)
    (hz : HasDerivAt z (2*k) s) (hk : (3/4 : ℝ) ≤ k) :
    HasDerivAt (fun t => z t-t) (2*k-1) s ∧ (1/2 : ℝ) ≤ 2*k-1 := by
  constructor
  · exact hz.sub (hasDerivAt_id s)
  · linarith


theorem normalized_curvature_lower_bound (alpha z b k K : ℝ)
    (ha : 0 < alpha) (hz : 0 < z) (hb : 0 < b) (hk : (3/4 : ℝ) ≤ k)
    (hK : 8*k^2/(alpha^2*z^2*b^2) ≤ K) :
    9/(2*alpha^2*z^2*b^2) ≤ K := by
  have hd : 0 < alpha^2*z^2*b^2 := by positivity
  have hn : (9/2 : ℝ) ≤ 8*k^2 := by nlinarith [sq_nonneg (k-3/4)]
  have hh := div_le_div_of_nonneg_right hn (le_of_lt hd)
  have he : (9/2 : ℝ)/(alpha^2*z^2*b^2)=9/(2*alpha^2*z^2*b^2) := by ring
  rw [he] at hh
  exact hh.trans hK

def curvatureLowerBound (alpha b0 z0 t : ℝ) : ℝ :=
  9*Real.exp (2*(z0/3+5/2)*t+t^2/2)/(2*alpha^2*b0^2*(z0+2*t)^2)

theorem radius_ratio_tendsto (z0 : ℝ) (hz0 : 0 ≤ z0) :
    Tendsto (fun t : ℝ => t/(z0+2*t)) atTop (𝓝 (1/2 : ℝ)) := by
  have hq : Tendsto (fun t : ℝ => z0/t) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop tendsto_id
  have hd : Tendsto (fun t : ℝ => z0/t+2) atTop (𝓝 (2 : ℝ)) := by
    simpa using hq.add_const 2
  have hh : Tendsto (fun t : ℝ => 1/(z0/t+2)) atTop (𝓝 (1/2 : ℝ)) :=
    tendsto_const_nhds.div hd (by norm_num)
  apply hh.congr'
  filter_upwards [eventually_gt_atTop (0 : ℝ)] with t ht
  have ht0 : t ≠ 0 := ne_of_gt ht
  have hn : z0+2*t ≠ 0 := by positivity
  field_simp
  <;> ring

theorem simple_curvature_lower_bound_tendsto (alpha b0 z0 : ℝ)
    (ha : 0 < alpha) (hb : 0 < b0) (hz0 : 0 ≤ z0) :
    Tendsto (fun t : ℝ => 9*Real.exp (2*t)/(2*alpha^2*b0^2*(z0+2*t)^2)) atTop atTop := by
  have hr := radius_ratio_tendsto z0 hz0
  have hfac : Tendsto (fun t : ℝ => (18/(alpha^2*b0^2))*(t/(z0+2*t))^2)
      atTop (𝓝 ((18/(alpha^2*b0^2))*(1/2 : ℝ)^2)) :=
    (hr.pow 2).const_mul _
  have hx : Tendsto (fun t : ℝ => 2*t) atTop atTop :=
    Tendsto.const_mul_atTop (by norm_num) tendsto_id
  have he : Tendsto (fun t : ℝ => Real.exp (2*t)/(2*t)^2) atTop atTop :=
    (Real.tendsto_exp_div_pow_atTop 2).comp hx
  have hp : 0 < (18/(alpha^2*b0^2))*(1/2 : ℝ)^2 := by positivity
  apply (Tendsto.pos_mul_atTop hp hfac he).congr'
  filter_upwards [eventually_gt_atTop (0 : ℝ)] with t ht
  have ht0 : t ≠ 0 := ne_of_gt ht
  have hn : z0+2*t ≠ 0 := by positivity
  have ha0 := ne_of_gt ha
  have hb0 := ne_of_gt hb
  field_simp
  <;> ring

theorem curvature_lower_bound_tendsto (alpha b0 z0 : ℝ)
    (ha : 0 < alpha) (hb : 0 < b0) (hz0 : 0 ≤ z0) :
    Tendsto (curvatureLowerBound alpha b0 z0) atTop atTop := by
  apply tendsto_atTop_mono' atTop _ (simple_curvature_lower_bound_tendsto alpha b0 z0 ha hb hz0)
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with t ht
  have hexp : 2*t ≤ 2*(z0/3+5/2)*t+t^2/2 := by
    nlinarith [mul_nonneg hz0 ht, sq_nonneg t]
  have he := Real.exp_le_exp.mpr hexp
  unfold curvatureLowerBound
  exact div_le_div_of_nonneg_right (by linarith)
    (by positivity : 0 ≤ 2*alpha^2*b0^2*(z0+2*t)^2)

theorem conditional_curvature_diverges (K : ℝ → ℝ) (alpha b0 z0 : ℝ)
    (ha : 0 < alpha) (hb : 0 < b0) (hz0 : 0 ≤ z0)
    (hK : ∀ᶠ t in atTop, curvatureLowerBound alpha b0 z0 t ≤ K t) :
    Tendsto K atTop atTop :=
  tendsto_atTop_mono' atTop hK (curvature_lower_bound_tendsto alpha b0 z0 ha hb hz0)

end
end SWaveInterior
