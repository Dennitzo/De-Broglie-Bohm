import Mathlib

/-! A centered massless finite Fourier mode on fixed de Sitter spacetime.
The exact width is verified with real derivatives and its force equation.
H and k are positive constants; a is a positive scale factor. -/

namespace GravityDeSitter

noncomputable section

def variance (H k a : ℝ) : ℝ := H^2/(2*k^3)+1/(2*k*a^2)
def width (H k a : ℝ) : ℝ := Real.sqrt (variance H k a)

def modeQuantumPotential (a sigma q : ℝ) : ℝ :=
  1/(4*a^3*sigma^2)-q^2/(8*a^3*sigma^4)

theorem quantum_potential_hasDerivAt (a sigma q : ℝ) (ha : a ≠ 0) (hs : sigma ≠ 0) :
    HasDerivAt (modeQuantumPotential a sigma) (-q/(4*a^3*sigma^4)) q := by
  have h := (hasDerivAt_const q (1/(4*a^3*sigma^2))).sub
    (((hasDerivAt_id q).pow 2).div_const (8*a^3*sigma^4))
  convert h using 1 <;> dsimp [modeQuantumPotential] <;> field_simp <;> ring

theorem quantum_acceleration_identity (a sigma q : ℝ) (ha : a ≠ 0) (hs : sigma ≠ 0) :
    -deriv (modeQuantumPotential a sigma) q/a^3=q/(4*a^6*sigma^4) := by
  rw [(quantum_potential_hasDerivAt a sigma q ha hs).deriv]
  field_simp
  <;> ring

theorem variance_positive (H k a : ℝ) (hH : 0 < H) (hk : 0 < k) (ha : 0 < a) :
    0 < variance H k a := by unfold variance; positivity

theorem width_positive (H k a : ℝ) (hH : 0 < H) (hk : 0 < k) (ha : 0 < a) :
    0 < width H k a := Real.sqrt_pos.mpr (variance_positive H k a hH hk ha)

theorem width_squared (H k a : ℝ) (hH : 0 < H) (hk : 0 < k) (ha : 0 < a) :
    (width H k a)^2=H^2/(2*k^3)+1/(2*k*a^2) := by
  exact Real.sq_sqrt (le_of_lt (variance_positive H k a hH hk ha))

theorem variance_hasDerivAt (a : ℝ → ℝ) (t H k : ℝ)
    (ha : HasDerivAt a (H*a t) t) (ha0 : a t ≠ 0) (hk : k ≠ 0) :
    HasDerivAt (fun u => variance H k (a u)) (-H/(k*(a t)^2)) t := by
  have hd := (ha.pow 2).const_mul (2*k)
  have hdn : 2*k*(a t)^2 ≠ 0 := by positivity
  have h := ((hasDerivAt_const t 1).div hd hdn).const_add (H^2/(2*k^3))
  convert h using 1 <;> dsimp [variance] <;> field_simp <;> ring

theorem width_hasDerivAt (a : ℝ → ℝ) (t H k : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hH : 0 < H) (hk : 0 < k) (ha0 : 0 < a t) :
    HasDerivAt (fun u => width H k (a u)) (-H/(2*k*(a t)^2*width H k (a t))) t := by
  have hv := variance_hasDerivAt a t H k ha (ne_of_gt ha0) (ne_of_gt hk)
  have hn := ne_of_gt (variance_positive H k (a t) hH hk ha0)
  have hw := ne_of_gt (width_positive H k (a t) hH hk ha0)
  convert hv.sqrt hn using 1 <;> dsimp [width] <;> field_simp <;> ring

theorem width_velocity_hasDerivAt (a sigma : ℝ → ℝ) (t H k : ℝ)
    (ha : HasDerivAt a (H*a t) t)
    (hsigma : HasDerivAt sigma (-H/(2*k*(a t)^2*sigma t)) t)
    (ha0 : a t ≠ 0) (hs : sigma t ≠ 0) (hk : k ≠ 0) :
    HasDerivAt (fun u => -H/(2*k*(a u)^2*sigma u))
      (H^2/(k*(a t)^2*sigma t)-H^2/(4*k^2*(a t)^4*(sigma t)^3)) t := by
  have hd := ((ha.pow 2).const_mul (2*k)).mul hsigma
  have hdn : 2*k*(a t)^2*sigma t ≠ 0 := by positivity
  convert (hasDerivAt_const t (-H)).div hd hdn using 1 <;> dsimp <;> field_simp <;> ring

theorem width_ermakov_algebra (H k a sigma : ℝ)
    (ha : a ≠ 0) (hk : k ≠ 0) (hs : sigma ≠ 0)
    (hr : sigma^2=H^2/(2*k^3)+1/(2*k*a^2)) :
    (H^2/(k*a^2*sigma)-H^2/(4*k^2*a^4*sigma^3))+
      3*H*(-H/(2*k*a^2*sigma))+(k^2/a^2)*sigma=1/(4*a^6*sigma^3) := by
  have hp : 2*k^3*a^2*sigma^2-H^2*a^2-k^2=0 := by
    rw [hr]
    field_simp
    <;> ring
  have hpoly : 4*k^4*a^4*sigma^4-2*H^2*k*a^4*sigma^2-H^2*a^2-k^2=0 := by
    linear_combination (2*k*a^2*sigma^2+1)*hp
  field_simp
  nlinarith [hpoly]

theorem de_sitter_width_ermakov (a : ℝ → ℝ) (H k t : ℝ)
    (hH : 0 < H) (hk : 0 < k) (ha0 : ∀ u, 0 < a u)
    (ha : ∀ u, HasDerivAt a (H*a u) u) :
    deriv (fun u => deriv (fun z => width H k (a z)) u) t+
      3*H*deriv (fun u => width H k (a u)) t+
      (k^2/(a t)^2)*width H k (a t)=1/(4*(a t)^6*(width H k (a t))^3) := by
  have hd : ∀ u, deriv (fun z => width H k (a z)) u = -H/(2*k*(a u)^2*width H k (a u)) := by
    intro u
    exact (width_hasDerivAt a u H k (ha u) hH hk (ha0 u)).deriv
  simp_rw [hd]
  rw [(width_velocity_hasDerivAt a (fun u => width H k (a u)) t H k (ha t)
    (width_hasDerivAt a t H k (ha t) hH hk (ha0 t)) (ne_of_gt (ha0 t))
    (ne_of_gt (width_positive H k (a t) hH hk (ha0 t))) (ne_of_gt hk)).deriv]
  exact width_ermakov_algebra H k (a t) (width H k (a t)) (ne_of_gt (ha0 t))
    (ne_of_gt hk) (ne_of_gt (width_positive H k (a t) hH hk (ha0 t)))
    (width_squared H k (a t) hH hk (ha0 t))

theorem frozen_variance_lower_bound (H k a : ℝ)
    (hH : 0 < H) (hk : 0 < k) (ha : 0 < a) :
    H^2/(2*k^3) < (width H k a)^2 := by
  rw [width_squared H k a hH hk ha]
  have h : 0 < 1/(2*k*a^2) := by positivity
  linarith

theorem width_freeze_error_bound (H k a sigma sigmaInf : ℝ)
    (ha : 0 < a) (hk : 0 < k) (hs : 0 ≤ sigma) (hsi : 0 ≤ sigmaInf)
    (hr : sigma^2=H^2/(2*k^3)+1/(2*k*a^2)) (hi : sigmaInf^2=H^2/(2*k^3)) :
    0 ≤ sigma-sigmaInf ∧ (sigma-sigmaInf)^2 ≤ 1/(2*k*a^2) := by
  have hpos : 0 < 1/(2*k*a^2) := by positivity
  have hord : sigmaInf ≤ sigma := by nlinarith
  constructor
  · linarith
  · nlinarith [mul_nonneg hsi (by linarith : 0 ≤ sigma-sigmaInf)]

end

end GravityDeSitter
