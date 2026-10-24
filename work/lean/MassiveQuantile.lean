import Mathlib

/-! Quantiles with absorbed boundary mass. Joint differentiability is explicit. -/
namespace MassiveQuantile
noncomputable section

def absorbed (boundaryFlux : ℝ → ℝ) (t : ℝ) : ℝ := -∫ s in (0 : ℝ)..t, boundaryFlux s
def cumulative (boundaryFlux : ℝ → ℝ) (density : ℝ → ℝ → ℝ) (a t r : ℝ) : ℝ :=
  absorbed boundaryFlux t+∫ x in a..r, density t x
def gradient (jt density : ℝ) : (ℝ × ℝ) →L[ℝ] ℝ :=
  (-jt) • ContinuousLinearMap.fst ℝ ℝ ℝ + density • ContinuousLinearMap.snd ℝ ℝ ℝ

theorem absorbed_mass_hasDerivAt (j : ℝ → ℝ) (t : ℝ) (hj : Continuous j) :
    HasDerivAt (absorbed j) (-j t) t := by
  exact (intervalIntegral.integral_hasDerivAt_right (hj.intervalIntegrable 0 t)
    (hj.stronglyMeasurableAtFilter MeasureTheory.volume (nhds t)) hj.continuousAt).neg

theorem cumulative_radial_hasDerivAt (j : ℝ → ℝ) (q : ℝ → ℝ → ℝ)
    (a t r : ℝ) (hq : Continuous (q t)) :
    HasDerivAt (cumulative j q a t) (q t r) r := by
  exact (intervalIntegral.integral_hasDerivAt_right (hq.intervalIntegrable a r)
    (hq.stronglyMeasurableAtFilter MeasureTheory.volume (nhds r)) hq.continuousAt).const_add (absorbed j t)

theorem cumulative_time_hasDerivAt (j : ℝ → ℝ) (q : ℝ → ℝ → ℝ)
    (a t r jr : ℝ) (hj : Continuous j)
    (hbalance : HasDerivAt (fun s => ∫ x in a..r, q s x) (j t-jr) t) :
    HasDerivAt (fun s => cumulative j q a s r) (-jr) t := by
  convert (absorbed_mass_hasDerivAt j t hj).add hbalance using 1 <;>
    dsimp [cumulative] <;> (try funext s) <;> (try dsimp) <;> ring

theorem cumulative_joint_derivative (A : ℝ → ℝ) (I : ℝ × ℝ → ℝ)
    (t r ja jr q : ℝ) (hA : HasDerivAt A (-ja) t)
    (hI : HasFDerivAt I ((ja-jr) • ContinuousLinearMap.fst ℝ ℝ ℝ+
      q • ContinuousLinearMap.snd ℝ ℝ ℝ) (t,r)) :
    HasFDerivAt (fun p : ℝ × ℝ => A p.1+I p) (gradient jr q) (t,r) := by
  have ha := hA.hasFDerivAt.comp (t,r) (ContinuousLinearMap.fst ℝ ℝ ℝ).hasFDerivAt
  have hall := ha.add hI
  have hd : (ContinuousLinearMap.toSpanSingleton ℝ (-ja)).comp
      (ContinuousLinearMap.fst ℝ ℝ ℝ)+
      ((ja-jr) • ContinuousLinearMap.fst ℝ ℝ ℝ+
        q • ContinuousLinearMap.snd ℝ ℝ ℝ)=gradient jr q := by
    apply ContinuousLinearMap.ext
    intro p
    simp [gradient]
    ring
  rw [hd] at hall
  exact hall

theorem quantile_along_guidance_hasDerivAt (F : ℝ × ℝ → ℝ) (r : ℝ → ℝ)
    (t jr q : ℝ) (hq : q ≠ 0)
    (hF : HasFDerivAt F (gradient jr q) (t,r t))
    (hr : HasDerivAt r (jr/q) t) :
    HasDerivAt (fun s => F (s,r s)) 0 t := by
  have hc := (hasDerivAt_id t).prodMk hr
  convert hF.comp_hasDerivAt t hc using 1
  simp [gradient]
  field_simp
  ring

theorem quantile_constant_on_interval (F : ℝ × ℝ → ℝ) (r j q : ℝ → ℝ)
    (a b t t₀ : ℝ)
    (hq : ∀ s ∈ Set.Ioo a b, q s ≠ 0)
    (hF : ∀ s ∈ Set.Ioo a b, HasFDerivAt F (gradient (j s) (q s)) (s,r s))
    (hr : ∀ s ∈ Set.Ioo a b, HasDerivAt r (j s/q s) s)
    (ht : t ∈ Set.Ioo a b) (ht₀ : t₀ ∈ Set.Ioo a b) :
    F (t,r t)=F (t₀,r t₀) := by
  have hd : ∀ s ∈ Set.Ioo a b, HasDerivAt (fun u => F (u,r u)) 0 s := by
    intro s hs
    exact quantile_along_guidance_hasDerivAt F r s (j s) (q s) (hq s hs) (hF s hs) (hr s hs)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun s hs => (hd s hs).differentiableAt.differentiableWithinAt)
    (fun s hs => (hd s hs).deriv) ht ht₀

theorem quantile_constancy_forces_guidance (F : ℝ × ℝ → ℝ) (r : ℝ → ℝ)
    (t jr q speed : ℝ) (hq : q ≠ 0)
    (hF : HasFDerivAt F (gradient jr q) (t,r t))
    (hr : HasDerivAt r speed t)
    (hconstant : HasDerivAt (fun s => F (s,r s)) 0 t) : speed=jr/q := by
  have hc := hF.comp_hasDerivAt t ((hasDerivAt_id t).prodMk hr)
  have he := hc.unique hconstant
  simp [gradient] at he
  apply (eq_div_iff hq).2
  linarith

theorem absorbed_bulk_outer_mass_derivative (A N B : ℝ → ℝ) (t ja jb : ℝ)
    (hA : HasDerivAt A (-ja) t) (hN : HasDerivAt N (ja-jb) t)
    (hB : HasDerivAt B jb t) : HasDerivAt (fun s => A s+N s+B s) 0 t := by
  convert (hA.add hN).add hB using 1 <;> dsimp <;> ring

theorem absorbed_bulk_outer_mass_constant (A N B ja jb : ℝ → ℝ)
    (a b t t₀ : ℝ)
    (hA : ∀ s ∈ Set.Ioo a b, HasDerivAt A (-ja s) s)
    (hN : ∀ s ∈ Set.Ioo a b, HasDerivAt N (ja s-jb s) s)
    (hB : ∀ s ∈ Set.Ioo a b, HasDerivAt B (jb s) s)
    (ht : t ∈ Set.Ioo a b) (ht₀ : t₀ ∈ Set.Ioo a b) :
    A t+N t+B t=A t₀+N t₀+B t₀ := by
  have hd : ∀ s ∈ Set.Ioo a b, HasDerivAt (fun u => A u+N u+B u) 0 s := by
    intro s hs
    exact absorbed_bulk_outer_mass_derivative A N B s (ja s) (jb s) (hA s hs) (hN s hs) (hB s hs)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun s hs => (hd s hs).differentiableAt.differentiableWithinAt)
    (fun s hs => (hd s hs).deriv) ht ht₀

end
end MassiveQuantile
