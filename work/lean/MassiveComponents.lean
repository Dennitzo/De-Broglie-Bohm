import Mathlib

/-! Real components of the specified massive Dirac equation in PG coordinates. -/
namespace MassiveComponents
noncomputable section

def componentDensity (a b : ℝ) : ℝ := a^2+b^2
def density (a b c e : ℝ) : ℝ := componentDensity a b+componentDensity c e
def flux (v a b c e : ℝ) : ℝ := (1-v)*componentDensity a b+(-1-v)*componentDensity c e
def transfer (mu a b c e : ℝ) : ℝ := 2*mu*(a*e-b*c)

theorem complex_mass_term_components (mu c e : ℝ) :
    (-Complex.I*(mu : ℂ)*((c : ℂ)+(e : ℂ)*Complex.I)).re=mu*e ∧
    (-Complex.I*(mu : ℂ)*((c : ℂ)+(e : ℂ)*Complex.I)).im= -mu*c := by
  simp [Complex.mul_re,Complex.mul_im]

theorem component_density_hasDerivAt (a b : ℝ → ℝ) (x da db : ℝ)
    (ha : HasDerivAt a da x) (hb : HasDerivAt b db x) :
    HasDerivAt (fun t => componentDensity (a t) (b t)) (2*(a x*da+b x*db)) x := by
  convert (ha.pow 2).add (hb.pow 2) using 1 <;>
    dsimp [componentDensity] <;> (try funext t) <;> (try dsimp) <;> ring

theorem upper_density_time_derivative (a b : ℝ → ℝ) (t v vp mu ar br c e : ℝ)
    (ha : HasDerivAt a ((v-1)*ar+vp/2*a t+mu*e) t)
    (hb : HasDerivAt b ((v-1)*br+vp/2*b t-mu*c) t) :
    HasDerivAt (fun s => componentDensity (a s) (b s))
      (2*(v-1)*(a t*ar+b t*br)+vp*componentDensity (a t) (b t)+transfer mu (a t) (b t) c e) t := by
  convert component_density_hasDerivAt a b t _ _ ha hb using 1
  unfold componentDensity transfer
  ring

theorem lower_density_time_derivative (c e : ℝ → ℝ) (t v vp mu cr er a b : ℝ)
    (hc : HasDerivAt c ((v+1)*cr+vp/2*c t+mu*b) t)
    (he : HasDerivAt e ((v+1)*er+vp/2*e t-mu*a) t) :
    HasDerivAt (fun s => componentDensity (c s) (e s))
      (2*(v+1)*(c t*cr+e t*er)+vp*componentDensity (c t) (e t)-transfer mu a b (c t) (e t)) t := by
  convert component_density_hasDerivAt c e t _ _ hc he using 1
  unfold componentDensity transfer
  ring

theorem upper_flux_space_derivative (v a b : ℝ → ℝ) (r vp ar br : ℝ)
    (hv : HasDerivAt v vp r) (ha : HasDerivAt a ar r) (hb : HasDerivAt b br r) :
    HasDerivAt (fun x => (1-v x)*componentDensity (a x) (b x))
      (-vp*componentDensity (a r) (b r)+2*(1-v r)*(a r*ar+b r*br)) r := by
  convert ((hasDerivAt_const r (1 : ℝ)).sub hv).mul
    (component_density_hasDerivAt a b r ar br ha hb) using 1 <;> dsimp <;> ring

theorem lower_flux_space_derivative (v c e : ℝ → ℝ) (r vp cr er : ℝ)
    (hv : HasDerivAt v vp r) (hc : HasDerivAt c cr r) (he : HasDerivAt e er r) :
    HasDerivAt (fun x => (-1-v x)*componentDensity (c x) (e x))
      (-vp*componentDensity (c r) (e r)-2*(1+v r)*(c r*cr+e r*er)) r := by
  convert ((hasDerivAt_const r (-1 : ℝ)).sub hv).mul
    (component_density_hasDerivAt c e r cr er hc he) using 1 <;> dsimp <;> ring

theorem upper_mass_exchange (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu ar br : ℝ) (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hat : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+mu*e t r) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r-mu*c t r) t) :
    deriv (fun s => componentDensity (a s r) (b s r)) t+
      deriv (fun x => (1-v x)*componentDensity (a t x) (b t x)) r =
      transfer mu (a t r) (b t r) (c t r) (e t r) := by
  rw [(upper_density_time_derivative (fun s => a s r) (fun s => b s r)
    t (v r) vp mu ar br (c t r) (e t r) hat hbt).deriv,
    (upper_flux_space_derivative v (a t) (b t) r vp ar br hv har hbr).deriv]
  ring

theorem lower_mass_exchange (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu cr er : ℝ) (hv : HasDerivAt v vp r)
    (hcr : HasDerivAt (c t) cr r) (her : HasDerivAt (e t) er r)
    (hct : HasDerivAt (fun s => c s r) ((v r+1)*cr+vp/2*c t r+mu*b t r) t)
    (het : HasDerivAt (fun s => e s r) ((v r+1)*er+vp/2*e t r-mu*a t r) t) :
    deriv (fun s => componentDensity (c s r) (e s r)) t+
      deriv (fun x => (-1-v x)*componentDensity (c t x) (e t x)) r =
      -transfer mu (a t r) (b t r) (c t r) (e t r) := by
  rw [(lower_density_time_derivative (fun s => c s r) (fun s => e s r)
    t (v r) vp mu cr er (a t r) (b t r) hct het).deriv,
    (lower_flux_space_derivative v (c t) (e t) r vp cr er hv hcr her).deriv]
  ring

theorem massive_dirac_continuity (a b c e : ℝ → ℝ → ℝ) (v : ℝ → ℝ)
    (t r vp mu ar br cr er : ℝ) (hv : HasDerivAt v vp r)
    (har : HasDerivAt (a t) ar r) (hbr : HasDerivAt (b t) br r)
    (hcr : HasDerivAt (c t) cr r) (her : HasDerivAt (e t) er r)
    (hat : HasDerivAt (fun s => a s r) ((v r-1)*ar+vp/2*a t r+mu*e t r) t)
    (hbt : HasDerivAt (fun s => b s r) ((v r-1)*br+vp/2*b t r-mu*c t r) t)
    (hct : HasDerivAt (fun s => c s r) ((v r+1)*cr+vp/2*c t r+mu*b t r) t)
    (het : HasDerivAt (fun s => e s r) ((v r+1)*er+vp/2*e t r-mu*a t r) t) :
    deriv (fun s => density (a s r) (b s r) (c s r) (e s r)) t+
      deriv (fun x => flux (v x) (a t x) (b t x) (c t x) (e t x)) r=0 := by
  have htu := upper_density_time_derivative (fun s => a s r) (fun s => b s r)
    t (v r) vp mu ar br (c t r) (e t r) hat hbt
  have htd := lower_density_time_derivative (fun s => c s r) (fun s => e s r)
    t (v r) vp mu cr er (a t r) (b t r) hct het
  have hru := upper_flux_space_derivative v (a t) (b t) r vp ar br hv har hbr
  have hrd := lower_flux_space_derivative v (c t) (e t) r vp cr er hv hcr her
  have ht := htu.add htd
  have hr := hru.add hrd
  change HasDerivAt (fun s => density (a s r) (b s r) (c s r) (e s r)) _ t at ht
  change HasDerivAt (fun x => flux (v x) (a t x) (b t x) (c t x) (e t x)) _ r at hr
  rw [ht.deriv,hr.deriv]
  ring

theorem mass_transfer_bound (mu a b c e : ℝ) :
    (transfer mu a b c e)^2 ≤ 4*mu^2*componentDensity a b*componentDensity c e := by
  unfold transfer componentDensity
  nlinarith [mul_nonneg (sq_nonneg mu) (sq_nonneg (a*c+b*e))]

end
end MassiveComponents
