import Mathlib

/-! Reduced Schwarzschild curvature and specified scalar trace anomaly.
Curvature identification and field-theoretic derivation are external
analytical inputs; the following statements prove their consequences. -/

namespace BlackHoleAnomaly

noncomputable section

def kretschmann (m r : ℝ) : ℝ := 48*m^2/r^6
def scalarTrace (m r : ℝ) : ℝ := kretschmann m r/(2880*Real.pi^2)
def leadingRicciShift (G m r : ℝ) : ℝ := -8*Real.pi*G*scalarTrace m r

theorem kretschmann_positive (m r : ℝ) (hm : 0 < m) (hr : 0 < r) :
    0 < kretschmann m r := by
  unfold kretschmann
  positivity

theorem scalar_trace_positive (m r : ℝ) (hm : 0 < m) (hr : 0 < r) :
    0 < scalarTrace m r := by
  unfold scalarTrace
  exact div_pos (kretschmann_positive m r hm hr) (by positivity)

theorem scalar_trace_reduced (m r : ℝ) :
    scalarTrace m r=m^2/(60*Real.pi^2*r^6) := by
  unfold scalarTrace kretschmann
  ring

theorem leading_ricci_shift_reduced (G m r : ℝ) :
    leadingRicciShift G m r= -2*G*m^2/(15*Real.pi*r^6) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  unfold leadingRicciShift
  rw [scalar_trace_reduced]
  by_cases hr : r=0
  · simp [hr]
  · field_simp
    <;> ring

theorem leading_ricci_shift_negative (G m r : ℝ)
    (hG : 0 < G) (hm : 0 < m) (hr : 0 < r) :
    leadingRicciShift G m r < 0 := by
  have hT := scalar_trace_positive m r hm hr
  unfold leadingRicciShift
  have hp := mul_pos (mul_pos (mul_pos (by norm_num : (0:ℝ)<8) Real.pi_pos) hG) hT
  nlinarith

theorem ricci_flat_trace_obstruction (G m r : ℝ)
    (hG : 0 < G) (hm : 0 < m) (hr : 0 < r) :
    ¬ ((0:ℝ) = -8*Real.pi*G*scalarTrace m r) := by
  have hn := leading_ricci_shift_negative G m r hG hm hr
  unfold leadingRicciShift at hn
  linarith

theorem zero_mass_trace (r : ℝ) : scalarTrace 0 r=0 := by
  simp [scalarTrace, kretschmann]

theorem schwarzschild_horizon_trace (m : ℝ) (hm : m ≠ 0) :
    scalarTrace m (2*m)=1/(3840*Real.pi^2*m^4) := by
  have hp : Real.pi ≠ 0 := ne_of_gt Real.pi_pos
  rw [scalar_trace_reduced]
  field_simp
  <;> ring

end

end BlackHoleAnomaly
