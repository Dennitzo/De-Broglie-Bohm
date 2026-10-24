import Mathlib

/-! Necessary trace conditions for the specified exact OR family.
The geometric endpoint limits are explicit hypotheses, not a construction of a full RSET. -/
namespace FullRSETTrace
noncomputable section
open Filter
open scoped Topology

def scalarR (A B C S : ℝ) : ℝ := -2*A-4*B-4*C+2*S
def kretschmann (A B C S : ℝ) : ℝ := 4*(A^2+2*B^2+2*C^2+S^2)
def ricciSquared (A B C S : ℝ) : ℝ := (A+2*B)^2+(-A-2*C)^2+2*(-B-C+S)^2
def anomalyQ (A B C S : ℝ) : ℝ := kretschmann A B C S-ricciSquared A B C S

theorem radial_scalar_divergence_actual_derivative (rho u dR : ℝ → ℝ)
    (r drho H ddR : ℝ) (hr : r ≠ 0) (hu : 0 < u r)
    (hρ : HasDerivAt rho drho r) (hw : HasDerivAt u H r)
    (hR : HasDerivAt dR ddR r) :
    Real.sqrt (u r)/(r^2*Real.exp (rho r))*
      deriv (fun x => x^2*Real.exp (rho x)*Real.sqrt (u x)*dR x) r =
    u r*ddR+(2*u r/r+u r*drho+H/2)*dR r := by
  have hs := hw.sqrt (ne_of_gt hu)
  have hf := ((((hasDerivAt_id r).pow 2).mul hρ.exp).mul hs).mul hR
  change HasDerivAt (fun x => x^2*Real.exp (rho x)*Real.sqrt (u x)*dR x) _ r at hf
  rw [hf.deriv]
  have hsq : (Real.sqrt (u r))^2=u r := Real.sq_sqrt (le_of_lt hu)
  have hn : Real.sqrt (u r) ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hu)
  dsimp
  field_simp
  rw [hsq]
  ring

theorem curvature_weight_scaling (w A B C S : ℝ) :
    w*scalarR A B C S=scalarR (w*A) (w*B) (w*C) (w*S) ∧
    w^2*anomalyQ A B C S=anomalyQ (w*A) (w*B) (w*C) (w*S) := by
  unfold scalarR anomalyQ kretschmann ricciSquared
  constructor <;> ring

theorem leading_radial_curvature_limits {ι : Type*} {l : Filter ι}
    (A B C S : ι → ℝ) (a : ℝ)
    (hA : Tendsto A l (𝓝 a)) (hB : Tendsto B l (𝓝 0))
    (hC : Tendsto C l (𝓝 0)) (hS : Tendsto S l (𝓝 0)) :
    Tendsto (fun t => scalarR (A t) (B t) (C t) (S t)) l (𝓝 (-2*a)) ∧
    Tendsto (fun t => anomalyQ (A t) (B t) (C t) (S t)) l (𝓝 (2*a^2)) := by
  have hR := (((hA.const_mul (-2)).sub (hB.const_mul 4)).sub (hC.const_mul 4)).add (hS.const_mul 2)
  have hK := (((hA.pow 2).add ((hB.pow 2).const_mul 2)).add
    ((hC.pow 2).const_mul 2)).add (hS.pow 2)
  have hRic := (((hA.add (hB.const_mul 2)).pow 2).add
    ((hA.neg.sub (hC.const_mul 2)).pow 2)).add
    ((((hB.neg.sub hC).add hS).pow 2).const_mul 2)
  constructor
  · simpa [scalarR] using hR
  · convert (hK.const_mul 4).sub hRic using 1 <;>
      dsimp [anomalyQ,kretschmann,ricciSquared] <;> ring

theorem weighted_scalar_vanishes (R : ℝ → ℝ) (r₅ : ℝ)
    (hR : Tendsto (fun w => w^5*R w) (𝓝[>] (0 : ℝ)) (𝓝 r₅)) :
    Tendsto (fun w => w^10*R w) (𝓝[>] (0 : ℝ)) (𝓝 0) ∧
    Tendsto (fun w => w^9*R w) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hi : Tendsto (fun w : ℝ => w) (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
    tendsto_id.mono_left inf_le_left
  constructor
  · convert (hi.pow 5).mul hR using 1 <;> (try norm_num) <;> (try funext w) <;> ring
  · convert (hi.pow 4).mul hR using 1 <;> (try norm_num) <;> (try funext w) <;> ring

theorem necessary_normalized_trace {ι : Type*} {l : Filter ι} [NeBot l]
    (weight R Q boxR : ι → ℝ) (q b A B : ℝ) (hq : q ≠ 0)
    (hR : Tendsto (fun t => weight t*R t) l (𝓝 0))
    (hQ : Tendsto (fun t => weight t*Q t) l (𝓝 A))
    (hBox : Tendsto (fun t => weight t*boxR t) l (𝓝 B))
    (hTrace : ∀ᶠ t in l, R t+q*(Q t+b*boxR t)=0) : A+b*B=0 := by
  have hh := hR.add ((hQ.add (hBox.const_mul b)).const_mul q)
  have hz : Tendsto (fun t => weight t*R t+q*(weight t*Q t+b*(weight t*boxR t))) l (𝓝 0) := by
    apply tendsto_const_nhds.congr'
    filter_upwards [hTrace] with t ht
    calc 0 = weight t*(R t+q*(Q t+b*boxR t)) := by rw [ht]; ring
      _ = _ := by ring
  have he := tendsto_nhds_unique hh hz
  have hp : q*(A+b*B)=0 := by simpa using he
  exact (mul_eq_zero.mp hp).resolve_left hq

theorem exact_OR_leading_trace_forces_scheme (R Q boxR : ℝ → ℝ) (q b a r₅ : ℝ)
    (hq : q ≠ 0) (ha : a ≠ 0)
    (hR : Tendsto (fun w => w^5*R w) (𝓝[>] (0 : ℝ)) (𝓝 r₅))
    (hQ : Tendsto (fun w => w^10*Q w) (𝓝[>] (0 : ℝ)) (𝓝 (32*a^4/81)))
    (hBox : Tendsto (fun w => w^10*boxR w) (𝓝[>] (0 : ℝ)) (𝓝 (-920*a^4/27)))
    (hTrace : ∀ᶠ w in 𝓝[>] (0 : ℝ), R w+q*(Q w+b*boxR w)=0) :
    b=4/345 := by
  have he := necessary_normalized_trace (fun w => w^10) R Q boxR q b
    (32*a^4/81) (-920*a^4/27) hq (weighted_scalar_vanishes R r₅ hR).1 hQ hBox hTrace
  have hp : a^4 ≠ 0 := pow_ne_zero 4 ha
  have hm : a^4*(4-345*b)=0 := by nlinarith [he]
  have hb := (mul_eq_zero.mp hm).resolve_left hp
  linarith

theorem canonical_trace_excludes_exact_OR (R Q boxR : ℝ → ℝ) (q a r₅ : ℝ)
    (hq : q ≠ 0) (ha : a ≠ 0)
    (hR : Tendsto (fun w => w^5*R w) (𝓝[>] (0 : ℝ)) (𝓝 r₅))
    (hQ : Tendsto (fun w => w^10*Q w) (𝓝[>] (0 : ℝ)) (𝓝 (32*a^4/81)))
    (hBox : Tendsto (fun w => w^10*boxR w) (𝓝[>] (0 : ℝ)) (𝓝 (-920*a^4/27))) :
    ¬ (∀ᶠ w in 𝓝[>] (0 : ℝ), R w+q*(Q w+0*boxR w)=0) := by
  intro ht
  have he := exact_OR_leading_trace_forces_scheme R Q boxR q 0 a r₅ hq ha hR hQ hBox ht
  norm_num at he

theorem conformal_subleading_coefficient (L r : ℝ) :
    4*L^4*(1/60 : ℝ)^3*(84*(1/60 : ℝ)-1)/(243*r^12) = L^4/(32805000*r^12) := by
  ring

theorem conformal_subleading_coefficient_positive (L r : ℝ) (hL : 0 < L) (hr : 0 < r) :
    0 < L^4/(32805000*r^12) := by positivity

theorem every_fixed_scheme_excludes_exact_conformal_OR
    (R Q boxR : ℝ → ℝ) (q b a r₅ L rD : ℝ)
    (hq : q ≠ 0) (ha : a ≠ 0) (hL : 0 < L) (hrD : 0 < rD)
    (hR : Tendsto (fun w => w^5*R w) (𝓝[>] (0 : ℝ)) (𝓝 r₅))
    (hQ : Tendsto (fun w => w^10*Q w) (𝓝[>] (0 : ℝ)) (𝓝 (32*a^4/81)))
    (hBox : Tendsto (fun w => w^10*boxR w) (𝓝[>] (0 : ℝ)) (𝓝 (-920*a^4/27)))
    (hNext : Tendsto (fun w => w^9*(Q w+(4/345 : ℝ)*boxR w))
      (𝓝[>] (0 : ℝ)) (𝓝 (L^4/(32805000*rD^12)))) :
    ¬ (∀ᶠ w in 𝓝[>] (0 : ℝ), R w+q*(Q w+b*boxR w)=0) := by
  intro ht
  have hb := exact_OR_leading_trace_forces_scheme R Q boxR q b a r₅ hq ha hR hQ hBox ht
  subst b
  have hh := (weighted_scalar_vanishes R r₅ hR).2.add (hNext.const_mul q)
  have hz : Tendsto (fun w => w^9*R w+q*(w^9*(Q w+(4/345 : ℝ)*boxR w)))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    apply tendsto_const_nhds.congr'
    filter_upwards [ht] with w hw
    calc 0 = w^9*(R w+q*(Q w+(4/345 : ℝ)*boxR w)) := by rw [hw]; ring
      _ = _ := by ring
  have he := tendsto_nhds_unique hh hz
  have hp : L^4/(32805000*rD^12) ≠ 0 := ne_of_gt (conformal_subleading_coefficient_positive L rD hL hrD)
  exact (mul_ne_zero hq hp) (by simpa using he)

end
end FullRSETTrace
