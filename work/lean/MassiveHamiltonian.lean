import Mathlib

/-! Classical massive PG Hamiltonian benchmark; separate from Bohmian guidance. -/
namespace MassiveHamiltonian
noncomputable section

def energy (mu p : ℝ) : ℝ := Real.sqrt (p^2+mu^2)
def hamiltonian (mu v p : ℝ) : ℝ := energy mu p-v*p
def velocity (mu v p : ℝ) : ℝ := p/energy mu p-v

theorem energy_positive (mu p : ℝ) (hmu : mu ≠ 0) : 0 < energy mu p := by
  unfold energy
  exact Real.sqrt_pos.2 (by nlinarith [sq_nonneg p,sq_pos_of_ne_zero hmu])

theorem classical_hamiltonian_hasDerivAt (r p v : ℝ → ℝ) (t mu vp : ℝ)
    (hmu : mu ≠ 0)
    (hr : HasDerivAt r (velocity mu (v (r t)) (p t)) t)
    (hp : HasDerivAt p (vp*p t) t) (hv : HasDerivAt v vp (r t)) :
    HasDerivAt (fun s => hamiltonian mu (v (r s)) (p s)) 0 t := by
  have hn : (p t)^2+mu^2 ≠ 0 := by nlinarith [sq_nonneg (p t),sq_pos_of_ne_zero hmu]
  have hE := ((hp.pow 2).add_const (mu^2)).sqrt hn
  have hH := hE.sub ((hv.comp t hr).mul hp)
  have he : energy mu (p t) ≠ 0 := ne_of_gt (energy_positive mu (p t) hmu)
  change Real.sqrt ((p t)^2+mu^2) ≠ 0 at he
  convert hH using 1 <;> dsimp [hamiltonian,velocity,energy] <;>
    (try funext s) <;> (try dsimp) <;> field_simp <;> ring

theorem classical_energy_constant (r p v vp : ℝ → ℝ) (mu a b t t₀ : ℝ)
    (hmu : mu ≠ 0)
    (hr : ∀ s ∈ Set.Ioo a b, HasDerivAt r (velocity mu (v (r s)) (p s)) s)
    (hp : ∀ s ∈ Set.Ioo a b, HasDerivAt p (vp s*p s) s)
    (hv : ∀ s ∈ Set.Ioo a b, HasDerivAt v (vp s) (r s))
    (ht : t ∈ Set.Ioo a b) (ht₀ : t₀ ∈ Set.Ioo a b) :
    hamiltonian mu (v (r t)) (p t)=hamiltonian mu (v (r t₀)) (p t₀) := by
  have hd : ∀ s ∈ Set.Ioo a b,
      HasDerivAt (fun u => hamiltonian mu (v (r u)) (p u)) 0 s := by
    intro s hs
    exact classical_hamiltonian_hasDerivAt r p v s mu (vp s) hmu (hr s hs) (hp s hs) (hv s hs)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun s hs => (hd s hs).differentiableAt.differentiableWithinAt)
    (fun s hs => (hd s hs).deriv) ht ht₀

theorem classical_timelike_identity (mu v p : ℝ) (hmu : mu ≠ 0) :
    1-(velocity mu v p+v)^2=mu^2/(p^2+mu^2) := by
  have hsq : (energy mu p)^2=p^2+mu^2 := by
    exact Real.sq_sqrt (by positivity)
  have he : energy mu p ≠ 0 := ne_of_gt (energy_positive mu p hmu)
  unfold velocity
  rw [sub_add_cancel,div_pow,hsq]
  have hn : p^2+mu^2 ≠ 0 := by nlinarith [sq_nonneg p,sq_pos_of_ne_zero hmu]
  field_simp
  ring

theorem classical_massive_strictly_timelike (mu v p : ℝ) (hmu : mu ≠ 0) :
    0 < 1-(velocity mu v p+v)^2 := by
  rw [classical_timelike_identity mu v p hmu]
  exact div_pos (sq_pos_of_ne_zero hmu) (by nlinarith [sq_nonneg p,sq_pos_of_ne_zero hmu])

end
end MassiveHamiltonian
