import Mathlib

/-! Specified static spherical Einstein-Polyakov model. The definitions
below are reduced source and geometric equations. General tensor and
quantum-state constructions are separate analytical inputs. -/

namespace BHBackSource

noncomputable section

def constraint (alpha r u v : ℝ) : ℝ := u^2+2*r*u*v+alpha*v^2-1
def uRate (alpha r u v : ℝ) : ℝ := v*((r^2+alpha)*u+alpha*r*v)/(r^2-alpha)
def vRate (alpha r u v : ℝ) : ℝ := -r*v*(2*u+r*v)/(r^2-alpha)
def scaledEnergy (alpha r v dv : ℝ) : ℝ := alpha*(2*dv+v^2)/r^2
def scaledRadialPressure (alpha r v : ℝ) : ℝ := -alpha*v^2/r^2

theorem radial_pressure_hasDerivAt (r v : ℝ → ℝ) (t alpha u dv : ℝ)
    (hr : HasDerivAt r u t) (hv : HasDerivAt v dv t) (hr0 : r t ≠ 0) :
    HasDerivAt (fun x => scaledRadialPressure alpha (r x) (v x))
      (-2*alpha*v t*dv/(r t)^2+2*alpha*(v t)^2*u/(r t)^3) t := by
  have hn : (r t)^2 ≠ 0 := pow_ne_zero 2 hr0
  convert ((hv.pow 2).const_mul (-alpha)).div (hr.pow 2) hn using 1 <;>
    dsimp [scaledRadialPressure] <;> field_simp <;> ring

theorem anisotropic_source_conservation (r v : ℝ → ℝ) (t alpha u dv : ℝ)
    (hr : HasDerivAt r u t) (hv : HasDerivAt v dv t) (hr0 : r t ≠ 0) :
    deriv (fun x => scaledRadialPressure alpha (r x) (v x)) t+
      v t*(scaledEnergy alpha (r t) (v t) dv+scaledRadialPressure alpha (r t) (v t))+
      2*u/(r t)*scaledRadialPressure alpha (r t) (v t)=0 := by
  rw [(radial_pressure_hasDerivAt r v t alpha u dv hr hv hr0).deriv]
  unfold scaledEnergy scaledRadialPressure
  field_simp
  <;> ring

theorem angular_einstein_identity (alpha r u v : ℝ) (hd : r^2-alpha ≠ 0) :
    r*(vRate alpha r u v+v^2)+uRate alpha r u v+u*v=0 := by
  unfold uRate vRate
  field_simp
  <;> ring

theorem energy_einstein_residual (alpha r u v : ℝ) (hd : r^2-alpha ≠ 0) :
    1-u^2-2*r*uRate alpha r u v-alpha*(2*vRate alpha r u v+v^2)=
      -constraint alpha r u v := by
  unfold uRate vRate constraint
  field_simp
  <;> ring

theorem radial_einstein_residual (alpha r u v : ℝ) :
    u^2-1+2*r*u*v+alpha*v^2=constraint alpha r u v := by
  unfold constraint
  ring

theorem einstein_equations_on_constraint (alpha r u v : ℝ)
    (hd : r^2-alpha ≠ 0) (hC : constraint alpha r u v=0) :
    1-u^2-2*r*uRate alpha r u v=alpha*(2*vRate alpha r u v+v^2) ∧
    u^2-1+2*r*u*v= -alpha*v^2 ∧
    r*(vRate alpha r u v+v^2)+uRate alpha r u v+u*v=0 := by
  have he := energy_einstein_residual alpha r u v hd
  have hp := radial_einstein_residual alpha r u v
  constructor
  · linarith
  constructor
  · linarith
  · exact angular_einstein_identity alpha r u v hd

theorem constraint_hasDerivAt (r u v : ℝ → ℝ) (t alpha : ℝ)
    (hd : (r t)^2-alpha ≠ 0)
    (hr : HasDerivAt r (u t) t)
    (hu : HasDerivAt u (uRate alpha (r t) (u t) (v t)) t)
    (hv : HasDerivAt v (vRate alpha (r t) (u t) (v t)) t) :
    HasDerivAt (fun x => constraint alpha (r x) (u x) (v x)) 0 t := by
  have hall := (((hu.pow 2).add (((hr.mul hu).mul hv).const_mul 2)).add
    ((hv.pow 2).const_mul alpha)).sub_const 1
  convert hall using 1 <;> dsimp [constraint, uRate, vRate] <;>
    (try funext x) <;> field_simp <;> ring

theorem constraint_preserved_on_interval (r u v : ℝ → ℝ) (alpha a b t t0 : ℝ)
    (ht : t ∈ Set.Ioo a b) (ht0 : t0 ∈ Set.Ioo a b)
    (hd : ∀ x ∈ Set.Ioo a b, (r x)^2-alpha ≠ 0)
    (hr : ∀ x ∈ Set.Ioo a b, HasDerivAt r (u x) x)
    (hu : ∀ x ∈ Set.Ioo a b, HasDerivAt u (uRate alpha (r x) (u x) (v x)) x)
    (hv : ∀ x ∈ Set.Ioo a b, HasDerivAt v (vRate alpha (r x) (u x) (v x)) x) :
    constraint alpha (r t) (u t) (v t)=constraint alpha (r t0) (u t0) (v t0) := by
  have hz : ∀ x ∈ Set.Ioo a b,
      HasDerivAt (fun z => constraint alpha (r z) (u z) (v z)) 0 x := by
    intro x hx
    exact constraint_hasDerivAt r u v x alpha (hd x hx) (hr x hx) (hu x hx) (hv x hx)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun x hx => (hz x hx).differentiableAt.differentiableWithinAt)
    (fun x hx => (hz x hx).deriv) ht ht0

theorem throat_constraint (alpha r v : ℝ) (ha : alpha ≠ 0)
    (hC : constraint alpha r 0 v=0) : v^2=1/alpha := by
  unfold constraint at hC
  apply (eq_div_iff ha).2
  nlinarith

theorem throat_radius_rate (alpha r v : ℝ) (ha : alpha ≠ 0)
    (hd : r^2-alpha ≠ 0) (hC : constraint alpha r 0 v=0) :
    uRate alpha r 0 v=r/(r^2-alpha) := by
  have hv := throat_constraint alpha r v ha hC
  unfold uRate
  have hh : v*((r^2+alpha)*0+alpha*r*v)=r := by
    have hp := congrArg (fun x : ℝ => alpha*r*x) hv
    dsimp at hp
    field_simp at hp
    nlinarith [hp]
  rw [hh]

theorem throat_flare_positive (alpha r v : ℝ) (ha : 0 < alpha)
    (hr : 0 < r) (hd : alpha < r^2) (hC : constraint alpha r 0 v=0) :
    0 < uRate alpha r 0 v := by
  rw [throat_radius_rate alpha r v (ne_of_gt ha) (ne_of_gt (sub_pos.mpr hd)) hC]
  exact div_pos hr (sub_pos.mpr hd)

theorem throat_null_energy_value (alpha r v : ℝ)
    (ha : alpha ≠ 0) (hr : r ≠ 0) (hd : r^2-alpha ≠ 0)
    (hC : constraint alpha r 0 v=0) :
    scaledEnergy alpha r v (vRate alpha r 0 v)+scaledRadialPressure alpha r v=
      -2/(r^2-alpha) := by
  have hv := throat_constraint alpha r v ha hC
  unfold scaledEnergy scaledRadialPressure vRate
  have hp : alpha*v^2=1 := by
    unfold constraint at hC
    nlinarith
  field_simp
  nlinarith [congrArg (fun x : ℝ => r^2*x) hp]

theorem throat_null_energy_negative (alpha r v : ℝ)
    (ha : 0 < alpha) (hr : 0 < r) (hd : alpha < r^2)
    (hC : constraint alpha r 0 v=0) :
    scaledEnergy alpha r v (vRate alpha r 0 v)+scaledRadialPressure alpha r v < 0 := by
  rw [throat_null_energy_value alpha r v (ne_of_gt ha) (ne_of_gt hr)
    (ne_of_gt (sub_pos.mpr hd)) hC]
  exact div_neg_of_neg_of_pos (by norm_num) (sub_pos.mpr hd)

end

end BHBackSource
