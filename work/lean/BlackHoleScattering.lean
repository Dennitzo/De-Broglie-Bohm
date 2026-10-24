import Mathlib

/-! Real-potential radial current and scattering-flux identities.
Asymptotic amplitudes are explicit matching premises; this module does
not derive the quantum Hawking spectrum or an asymptotic existence theorem. -/

namespace BlackHoleScattering

noncomputable section

def radialPotential (m r ell : ℝ) : ℝ :=
  (1-2*m/r)*(ell*(ell+1)/r^2+2*m/r^3)
def realCurrent (x y dx dy : ℝ) : ℝ := x*dy-y*dx
def complexCurrent (psi dpsi : ℂ) : ℝ := (star psi*dpsi).im

theorem exterior_lapse_positive (m r : ℝ) (hm : 0 < m) (hr : 2*m < r) :
    0 < 1-2*m/r := by
  have hr0 : 0 < r := by linarith
  have hdiv : 2*m/r < 1 := (div_lt_one hr0).2 hr
  linarith

theorem radial_potential_positive (m r ell : ℝ)
    (hm : 0 < m) (hr : 2*m < r) (hell : 0 ≤ ell) :
    0 < radialPotential m r ell := by
  have hr0 : 0 < r := by linarith
  have hl := exterior_lapse_positive m r hm hr
  unfold radialPotential
  apply mul_pos hl
  have h1 : 0 ≤ ell*(ell+1)/r^2 := by positivity
  have h2 : 0 < 2*m/r^3 := by positivity
  linarith

theorem curvature_coupling_vanishes (xi : ℝ) : xi*(0:ℝ)=0 := by ring

theorem current_real_imaginary_bridge (psi dpsi : ℂ) :
    complexCurrent psi dpsi=realCurrent psi.re psi.im dpsi.re dpsi.im := by
  simp [complexCurrent, realCurrent, Complex.mul_im]
  ring

theorem current_hasDerivAt (x y dx dy : ℝ → ℝ) (t V omega : ℝ)
    (hx : HasDerivAt x (dx t) t) (hy : HasDerivAt y (dy t) t)
    (hdx : HasDerivAt dx ((V-omega^2)*x t) t)
    (hdy : HasDerivAt dy ((V-omega^2)*y t) t) :
    HasDerivAt (fun u => realCurrent (x u) (y u) (dx u) (dy u)) 0 t := by
  convert (hx.mul hdy).sub (hy.mul hdx) using 1 <;> dsimp [realCurrent] <;> ring

theorem current_constant_on_interval (x y dx dy V : ℝ → ℝ) (omega a b t t0 : ℝ)
    (ht : t ∈ Set.Ioo a b) (ht0 : t0 ∈ Set.Ioo a b)
    (hx : ∀ u ∈ Set.Ioo a b, HasDerivAt x (dx u) u)
    (hy : ∀ u ∈ Set.Ioo a b, HasDerivAt y (dy u) u)
    (hdx : ∀ u ∈ Set.Ioo a b, HasDerivAt dx ((V u-omega^2)*x u) u)
    (hdy : ∀ u ∈ Set.Ioo a b, HasDerivAt dy ((V u-omega^2)*y u) u) :
    realCurrent (x t) (y t) (dx t) (dy t)=realCurrent (x t0) (y t0) (dx t0) (dy t0) := by
  have hz : ∀ u ∈ Set.Ioo a b,
      HasDerivAt (fun z => realCurrent (x z) (y z) (dx z) (dy z)) 0 u := by
    intro u hu
    exact current_hasDerivAt x y dx dy u (V u) omega (hx u hu) (hy u hu) (hdx u hu) (hdy u hu)
  exact isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
    (fun u hu => (hz u hu).differentiableAt.differentiableWithinAt)
    (fun u hu => (hz u hu).deriv) ht ht0

theorem two_wave_current (incoming outgoing : ℂ) (omega : ℝ) :
    complexCurrent (incoming+outgoing)
      ((-Complex.I)*(omega:ℂ)*incoming+Complex.I*(omega:ℂ)*outgoing)=
      omega*(Complex.normSq outgoing-Complex.normSq incoming) := by
  simp [complexCurrent, Complex.mul_im, Complex.mul_re, Complex.normSq_apply]
  ring

theorem ingoing_wave_current (amplitude : ℂ) (omega : ℝ) :
    complexCurrent amplitude ((-Complex.I)*(omega:ℂ)*amplitude)=
      -omega*Complex.normSq amplitude := by
  have hh := two_wave_current amplitude 0 omega
  simpa using hh

theorem amplitude_balance (incoming outgoing transmitted : ℂ) (omega : ℝ)
    (hw : omega ≠ 0)
    (hflux : omega*(Complex.normSq outgoing-Complex.normSq incoming)=
      -omega*Complex.normSq transmitted) :
    Complex.normSq incoming=Complex.normSq outgoing+Complex.normSq transmitted := by
  have hh : omega*(Complex.normSq incoming-Complex.normSq outgoing-Complex.normSq transmitted)=0 := by
    nlinarith [hflux]
  have hz := (mul_eq_zero.mp hh).resolve_left hw
  linarith

def reflection (incoming outgoing : ℂ) : ℝ := Complex.normSq outgoing/Complex.normSq incoming
def transmission (incoming transmitted : ℂ) : ℝ := Complex.normSq transmitted/Complex.normSq incoming

theorem reflection_transmission_sum (incoming outgoing transmitted : ℂ)
    (hi : Complex.normSq incoming ≠ 0)
    (hbalance : Complex.normSq incoming=Complex.normSq outgoing+Complex.normSq transmitted) :
    reflection incoming outgoing+transmission incoming transmitted=1 := by
  unfold reflection transmission
  rw [← add_div, ← hbalance]
  exact div_self hi

theorem reflection_transmission_bounds (incoming outgoing transmitted : ℂ)
    (hi : Complex.normSq incoming ≠ 0)
    (hbalance : Complex.normSq incoming=Complex.normSq outgoing+Complex.normSq transmitted) :
    0 ≤ reflection incoming outgoing ∧ reflection incoming outgoing ≤ 1 ∧
      0 ≤ transmission incoming transmitted ∧ transmission incoming transmitted ≤ 1 := by
  have hr : 0 ≤ reflection incoming outgoing :=
    div_nonneg (Complex.normSq_nonneg outgoing) (Complex.normSq_nonneg incoming)
  have ht : 0 ≤ transmission incoming transmitted :=
    div_nonneg (Complex.normSq_nonneg transmitted) (Complex.normSq_nonneg incoming)
  have hs := reflection_transmission_sum incoming outgoing transmitted hi hbalance
  constructor
  · exact hr
  constructor
  · linarith
  constructor
  · exact ht
  · linarith

theorem quasistatic_mass_cube_hasDerivAt (m : ℝ → ℝ) (t G p : ℝ)
    (hm0 : m t ≠ 0) (hm : HasDerivAt m (-G*p/(m t)^2) t) :
    HasDerivAt (fun u => (m u)^3) (-3*G*p) t := by
  convert hm.pow 3 using 1 <;> dsimp <;> field_simp <;> ring

theorem quasistatic_mass_cube_on_interval (m : ℝ → ℝ) (G p a b t t0 : ℝ)
    (ht : t ∈ Set.Ioo a b) (ht0 : t0 ∈ Set.Ioo a b)
    (hm0 : ∀ u ∈ Set.Ioo a b, m u ≠ 0)
    (hm : ∀ u ∈ Set.Ioo a b, HasDerivAt m (-G*p/(m u)^2) u) :
    (m t)^3=(m t0)^3-3*G*p*(t-t0) := by
  have hz : ∀ u ∈ Set.Ioo a b,
      HasDerivAt (fun z => (m z)^3+3*G*p*z) 0 u := by
    intro u hu
    convert (quasistatic_mass_cube_hasDerivAt m u G p (hm0 u hu) (hm u hu)).add
      ((hasDerivAt_id u).const_mul (3*G*p)) using 1 <;> dsimp <;> ring
  have hc : (m t)^3+3*G*p*t=(m t0)^3+3*G*p*t0 :=
    isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo a b).isPreconnected
      (fun u hu => (hz u hu).differentiableAt.differentiableWithinAt)
      (fun u hu => (hz u hu).deriv) ht ht0
  linarith

theorem quasistatic_surface_gravity_hasDerivAt (m : ℝ → ℝ) (t G p : ℝ)
    (hm0 : m t ≠ 0) (hm : HasDerivAt m (-G*p/(m t)^2) t) :
    HasDerivAt (fun u => 1/(4*m u)) (G*p/(4*(m t)^4)) t := by
  have hn : 4*m t ≠ 0 := mul_ne_zero (by norm_num) hm0
  convert (hasDerivAt_const t 1).div (hm.const_mul 4) hn using 1 <;>
    dsimp <;> field_simp <;> ring

theorem quasistatic_adiabaticity_ratio (m : ℝ → ℝ) (t G p : ℝ)
    (hG : 0 ≤ G) (hp : 0 ≤ p) (hm0 : m t ≠ 0)
    (hm : HasDerivAt m (-G*p/(m t)^2) t) :
    |deriv (fun u => 1/(4*m u)) t|/(1/(4*m t))^2=4*G*p/(m t)^2 := by
  rw [(quasistatic_surface_gravity_hasDerivAt m t G p hm0 hm).deriv]
  rw [abs_of_nonneg (by positivity : 0 ≤ G*p/(4*(m t)^4))]
  field_simp
  <;> ring

end

end BlackHoleScattering
