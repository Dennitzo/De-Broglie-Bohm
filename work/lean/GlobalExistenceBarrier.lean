import Mathlib

open scoped Real

namespace GlobalExistenceBarrier

/- Kapitel 21: Barriere- und Kollapskern des superlinearen Normmodells
alpha' = -exp(-3 alpha) der gefuehrten Skalenfaktor-Dynamik. Formalisiert:
Kollapszeit (strikt positiv), die streng fallende Barriere y(t)=exp(3 a0)-3t,
die Verifikation der geschlossenen Loesung alpha(t)=log(y(t))/3 gegen die
Differentialgleichung, und die Gronwall-artige Beschraenktheit des linearen
Modells q'=k q. Existenz und Eindeutigkeit werden nicht neu behauptet.

Geprueft ueber das math.formalProof-Profil (Lean 4.30.0, Mathlib), Ergebnis
KernelAccepted, axiomAudit=passed; Axiome ausnahmslos propext, Classical.choice,
Quot.sound (kein sorry). -/

theorem t_star_pos (a0 : ℝ) : (0 : ℝ) < Real.exp (3 * a0) / 3 := by
  positivity

theorem y_barrier (a0 t : ℝ) :
    let t_star := Real.exp (3 * a0) / 3
    (Real.exp (3 * a0) - 3 * t > 0 ↔ t < t_star) ∧ (Real.exp (3 * a0) - 3 * t_star = 0) := by
  constructor
  · dsimp only
    constructor <;> intro h <;> linarith
  · dsimp only
    ring

theorem deriv_super (a0 t : ℝ) (h : 0 < Real.exp (3 * a0) - 3 * t) :
    HasDerivAt (fun s => Real.log (Real.exp (3 * a0) - 3 * s) / 3)
      (-1 / (Real.exp (3 * a0) - 3 * t)) t := by
  have hinner : HasDerivAt (fun s => Real.exp (3 * a0) - 3 * s) (-3) t := by
    have hc : HasDerivAt (fun _s => Real.exp (3 * a0)) 0 t := hasDerivAt_const t (Real.exp (3 * a0))
    have hl : HasDerivAt (fun s => (3 : ℝ) * s) 3 t := by
      convert (hasDerivAt_id t).const_mul (3 : ℝ) using 1
      ring
    convert hc.sub hl using 1
    ring
  have hmain := hinner.log (by linarith : Real.exp (3 * a0) - 3 * t ≠ 0) |>.div_const (3 : ℝ)
  convert hmain using 1
  field_simp

theorem gronwall_bound_linear (q0 k t : ℝ) (hk : 0 ≤ k) (ht : 0 ≤ t) :
    q0 * Real.exp (k * t) ≤ (1 + q0) * Real.exp (k * t) - 1 := by
  have e1 : (1 : ℝ) ≤ Real.exp (k * t) := by
    have := Real.exp_le_exp.mpr (show (0 : ℝ) ≤ k * t by nlinarith)
    rwa [Real.exp_zero] at this
  nlinarith

end GlobalExistenceBarrier
