import Mathlib

/-! Finite multipole identities, not a renormalization prescription. -/
namespace FourDimModes
noncomputable section
open Filter
open scoped Topology BigOperators

def angularEigenvalue (l : ℕ) : ℝ := (l : ℝ)*(l+1)
def potential (rho r rxx : ℝ) (l : ℕ) : ℝ :=
  rxx/r+Real.exp (2*rho)*angularEigenvalue l/r^2

theorem multipole_reduced_wave_operator (rho r rx rxx f fx fxx ftt : ℝ)
    (l : ℕ) (hr : r ≠ 0) :
    -(r*ftt)+(rxx*f+2*rx*fx+r*fxx)-potential rho r rxx l*(r*f) =
      r*(-ftt+fxx+2*rx/r*fx-Real.exp (2*rho)*angularEigenvalue l/r^2*f) := by
  unfold potential
  field_simp
  <;> ring

theorem angular_barrier_nonnegative (rho r : ℝ) (l : ℕ) :
    0 ≤ Real.exp (2*rho)*angularEigenvalue l/r^2 := by
  unfold angularEigenvalue
  positivity

theorem angular_barrier_strict (rho r : ℝ) (l : ℕ) (hl : 0 < l) (hr : r ≠ 0) :
    0 < Real.exp (2*rho)*angularEigenvalue l/r^2 := by
  have hlr : (0 : ℝ) < l := by exact_mod_cast hl
  unfold angularEigenvalue
  exact div_pos (mul_pos (Real.exp_pos _) (mul_pos hlr (by positivity)))
    (sq_pos_of_ne_zero hr)

theorem adjacent_multipole_difference (rho r rxx : ℝ) (l : ℕ) :
    potential rho r rxx (l+1)-potential rho r rxx l =
      2*(l+1)*Real.exp (2*rho)/r^2 := by
  unfold potential angularEigenvalue
  push_cast
  ring

theorem potential_ge_swave (rho r rxx : ℝ) (l : ℕ) :
    potential rho r rxx 0 ≤ potential rho r rxx l := by
  have hp := angular_barrier_nonnegative rho r l
  simpa [potential, angularEigenvalue] using add_le_add_left hp (rxx/r)

theorem degeneracy_sum (L : ℕ) :
    (∑ l ∈ Finset.range (L+1), (2*l+1))=(L+1)^2 := by
  induction L with
  | zero => norm_num
  | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      ring

theorem weighted_degeneracy_sum (L : ℕ) :
    2*(∑ l ∈ Finset.range (L+1), (2*l+1)*l*(l+1))=L*(L+1)^2*(L+2) := by
  induction L with
  | zero => norm_num
  | succ n ih =>
      rw [Finset.sum_range_succ]
      nlinarith [ih]

theorem degeneracy_sum_tendsto :
    Tendsto (fun L : ℕ => ∑ l ∈ Finset.range (L+1), (2*l+1)) atTop atTop := by
  apply tendsto_atTop_mono (fun L => ?_) tendsto_id
  rw [degeneracy_sum]
  change L ≤ (L+1)^2
  nlinarith


def generalPotential (f r rxx xi R : ℝ) (l : ℕ) : ℝ :=
  rxx/r+f*(angularEigenvalue l/r^2+xi*R)

theorem nonminimal_reduced_wave_operator (fMetric r rx rxx f fx fxx ftt xi R : ℝ)
    (l : ℕ) (hr : r ≠ 0) :
    -(r*ftt)+(rxx*f+2*rx*fx+r*fxx)-generalPotential fMetric r rxx xi R l*(r*f) =
      r*(-ftt+fxx+2*rx/r*fx-fMetric*(angularEigenvalue l/r^2+xi*R)*f) := by
  unfold generalPotential
  field_simp
  <;> ring

theorem fixed_multipole_endpoint_limit {ι : Type*} {L : Filter ι}
    (x f r rxx R : ι → ℝ) (xi : ℝ) (l : ℕ)
    (hRad : Tendsto (fun t => (x t)^2*(rxx t/r t)) L (𝓝 0))
    (hAng : Tendsto (fun t => (x t)^2*f t/(r t)^2) L (𝓝 0))
    (hCurv : Tendsto (fun t => (x t)^2*f t*R t) L (𝓝 (-(1/8 : ℝ)))) :
    Tendsto (fun t => (x t)^2*generalPotential (f t) (r t) (rxx t) xi (R t) l)
      L (𝓝 (-xi/8)) := by
  have hh := (hRad.add (hAng.mul_const (angularEigenvalue l))).add (hCurv.mul_const xi)
  have he : 0+0*angularEigenvalue l+(-(1/8 : ℝ))*xi = -xi/8 := by ring
  rw [he] at hh
  convert hh using 1
  funext t
  unfold generalPotential
  ring

end
end FourDimModes
