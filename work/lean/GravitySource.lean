import Mathlib

/-! Actual continuity and Friedmann-constraint derivatives for a specified
FLRW truncation. The Einstein equations used as premises are explicit;
no curvature-tensor or general Bianchi proof is claimed here. -/

namespace GravitySource

noncomputable section

def rho (W q v : ℝ) : ℝ := (v^2+W*q^2)/2
def pressure (W q v : ℝ) : ℝ := v^2/2-W*q^2/6

theorem mode_energy_hasDerivAt (W q v : ℝ → ℝ) (t H acc : ℝ)
    (hW : HasDerivAt W (-2*H*W t) t) (hq : HasDerivAt q (v t) t)
    (hv : HasDerivAt v acc t) :
    HasDerivAt (fun u => rho (W u) (q u) (v u))
      (v t*(acc+3*H*v t+W t*q t)-3*H*(rho (W t) (q t) (v t)+pressure (W t) (q t) (v t))) t := by
  have h := ((hv.pow 2).add (hW.mul (hq.pow 2))).div_const 2
  convert h using 1 <;> dsimp [rho, pressure] <;> ring

theorem mode_source_balance (W q v : ℝ → ℝ) (t H acc : ℝ)
    (hW : HasDerivAt W (-2*H*W t) t) (hq : HasDerivAt q (v t) t)
    (hv : HasDerivAt v acc t) :
    deriv (fun u => rho (W u) (q u) (v u)) t+
      3*H*(rho (W t) (q t) (v t)+pressure (W t) (q t) (v t))=
      v t*(acc+3*H*v t+W t*q t) := by
  rw [(mode_energy_hasDerivAt W q v t H acc hW hq hv).deriv]
  ring

theorem quantum_source_work (W q v : ℝ → ℝ) (t H acc F : ℝ)
    (hW : HasDerivAt W (-2*H*W t) t) (hq : HasDerivAt q (v t) t)
    (hv : HasDerivAt v acc t) (hF : acc+3*H*v t+W t*q t=F) :
    deriv (fun u => rho (W u) (q u) (v u)) t+
      3*H*(rho (W t) (q t) (v t)+pressure (W t) (q t) (v t))=v t*F := by
  rw [mode_source_balance W q v t H acc hW hq hv, hF]

theorem de_sitter_quantum_work_identity (H k a sigma q : ℝ)
    (ha : a ≠ 0) (hk : k ≠ 0) (hs : sigma ≠ 0) :
    (-H*q/(2*k*a^2*sigma^2))*(q/(4*a^6*sigma^4))=
      -H*q^2/(8*k*a^8*sigma^6) := by
  field_simp
  <;> ring

theorem de_sitter_quantum_work_negative (H k a sigma q : ℝ)
    (hH : 0 < H) (hk : 0 < k) (ha : 0 < a) (hs : 0 < sigma) (hq : q ≠ 0) :
    (-H*q/(2*k*a^2*sigma^2))*(q/(4*a^6*sigma^4)) < 0 := by
  rw [de_sitter_quantum_work_identity H k a sigma q (ne_of_gt ha) (ne_of_gt hk) (ne_of_gt hs)]
  have hn : 0 < H*q^2 := by positivity
  have hd : 0 < 8*k*a^8*sigma^6 := by positivity
  exact div_neg_of_neg_of_pos (by nlinarith) hd

def quantumCorrection (a sigma : ℝ) : ℝ := 1/(8*a^6*sigma^2)

theorem quantum_correction_hasDerivAt (a sigma : ℝ → ℝ) (t H sd : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hsigma : HasDerivAt sigma sd t)
    (ha0 : a t ≠ 0) (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => quantumCorrection (a u) (sigma u))
      (-6*H*quantumCorrection (a t) (sigma t)-sd/(4*(a t)^6*(sigma t)^3)) t := by
  have hd := ((ha.pow 6).const_mul 8).mul (hsigma.pow 2)
  have hdn : 8*(a t)^6*(sigma t)^2 ≠ 0 := by positivity
  convert (hasDerivAt_const t 1).div hd hdn using 1 <;>
    dsimp [quantumCorrection] <;> field_simp <;> ring

theorem gaussian_bare_mean_source (W sigma sd : ℝ → ℝ) (t H a : ℝ)
    (hW : HasDerivAt W (-2*H*W t) t) (hsigma : HasDerivAt sigma (sd t) t)
    (hsd : HasDerivAt sd (-3*H*sd t-W t*sigma t+1/(4*a^6*(sigma t)^3)) t) :
    deriv (fun u => rho (W u) (sigma u) (sd u)) t+
      3*H*(rho (W t) (sigma t) (sd t)+pressure (W t) (sigma t) (sd t))=
      sd t/(4*a^6*(sigma t)^3) := by
  rw [mode_source_balance W sigma sd t H _ hW hsigma hsd]
  ring

theorem operator_mean_source_hasDerivAt (a W sigma sd : ℝ → ℝ) (t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hW : HasDerivAt W (-2*H*W t) t)
    (hsigma : HasDerivAt sigma (sd t) t)
    (hsd : HasDerivAt sd (-3*H*sd t-W t*sigma t+1/(4*(a t)^6*(sigma t)^3)) t)
    (ha0 : a t ≠ 0) (hs : sigma t ≠ 0) :
    HasDerivAt (fun u => rho (W u) (sigma u) (sd u)+quantumCorrection (a u) (sigma u))
      (-3*H*((rho (W t) (sigma t) (sd t)+quantumCorrection (a t) (sigma t))+
        (pressure (W t) (sigma t) (sd t)+quantumCorrection (a t) (sigma t)))) t := by
  have hr := mode_energy_hasDerivAt W sigma sd t H _ hW hsigma hsd
  have hd := quantum_correction_hasDerivAt a sigma t H (sd t) ha hsigma ha0 hs
  convert hr.add hd using 1 <;> dsimp <;> ring

theorem operator_mean_source_conserved (a W sigma sd : ℝ → ℝ) (t H : ℝ)
    (ha : HasDerivAt a (H*a t) t) (hW : HasDerivAt W (-2*H*W t) t)
    (hsigma : HasDerivAt sigma (sd t) t)
    (hsd : HasDerivAt sd (-3*H*sd t-W t*sigma t+1/(4*(a t)^6*(sigma t)^3)) t)
    (ha0 : a t ≠ 0) (hs : sigma t ≠ 0) :
    deriv (fun u => rho (W u) (sigma u) (sd u)+quantumCorrection (a u) (sigma u)) t+
      3*H*((rho (W t) (sigma t) (sd t)+quantumCorrection (a t) (sigma t))+
        (pressure (W t) (sigma t) (sd t)+quantumCorrection (a t) (sigma t)))=0 := by
  rw [(operator_mean_source_hasDerivAt a W sigma sd t H ha hW hsigma hsd ha0 hs).deriv]
  ring

def constraint (c Lambda H density : ℝ) : ℝ := H^2-c*density-Lambda/3

theorem friedmann_constraint_hasDerivAt (H density : ℝ → ℝ) (t c Lambda p drho : ℝ)
    (hH : HasDerivAt H (-(3*c/2)*(density t+p)) t)
    (hrho : HasDerivAt density drho t) :
    HasDerivAt (fun u => constraint c Lambda (H u) (density u))
      (-c*(drho+3*H t*(density t+p))) t := by
  have h := ((hH.pow 2).sub (hrho.const_mul c)).sub_const (Lambda/3)
  convert h using 1 <;> dsimp [constraint] <;> ring

theorem friedmann_constraint_zero_derivative (H density : ℝ → ℝ) (t c Lambda p : ℝ)
    (hH : HasDerivAt H (-(3*c/2)*(density t+p)) t)
    (hrho : HasDerivAt density (-3*H t*(density t+p)) t) :
    HasDerivAt (fun u => constraint c Lambda (H u) (density u)) 0 t := by
  convert friedmann_constraint_hasDerivAt H density t c Lambda p _ hH hrho using 1 <;> ring

theorem friedmann_constraint_preserved (H density pressure : ℝ → ℝ) (c Lambda t t0 : ℝ)
    (hH : ∀ u, HasDerivAt H (-(3*c/2)*(density u+pressure u)) u)
    (hrho : ∀ u, HasDerivAt density (-3*H u*(density u+pressure u)) u) :
    constraint c Lambda (H t) (density t)=constraint c Lambda (H t0) (density t0) := by
  have hz : ∀ u, HasDerivAt (fun z => constraint c Lambda (H z) (density z)) 0 u := by
    intro u
    exact friedmann_constraint_zero_derivative H density u c Lambda (pressure u) (hH u) (hrho u)
  exact is_const_of_deriv_eq_zero (fun u => (hz u).differentiableAt)
    (fun u => (hz u).deriv) t t0

end

end GravitySource
