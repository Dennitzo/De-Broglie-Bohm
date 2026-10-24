import Mathlib

/-! Exact FLRW quantum-source exchange, with actual derivative hypotheses. -/
namespace UnifiedWard
noncomputable section

def rhoQ (b u q : ℝ) := b*u^2*q
def pressureQ (b u q qa : ℝ) := b*u^2*(q-qa/3)
def rhoF (b v : ℝ) := v^2/(2*b)

theorem volume_derivative (a : ℝ → ℝ) (t h : ℝ)
    (ha : HasDerivAt a h t) :
    HasDerivAt (fun s => Real.exp (-3*a s))
      (-3*h*Real.exp (-3*a t)) t := by
  convert (ha.const_mul (-3)).exp using 1 <;> ring

theorem quantum_energy_derivative (u q : ℝ → ℝ) (t b h qa qp v : ℝ)
    (hu : HasDerivAt u (-3*h*u t) t)
    (hq : HasDerivAt q (qa*h+qp*v) t) :
    HasDerivAt (fun s => rhoQ b (u s) (q s))
      (b*(u t)^2*(qa*h+qp*v-6*h*q t)) t := by
  convert (((hu.pow 2).const_mul b).mul hq) using 1 <;>
    dsimp [rhoQ] <;> ring

theorem quantum_exchange (u q : ℝ → ℝ) (t b h qa qp v : ℝ)
    (hu : HasDerivAt u (-3*h*u t) t)
    (hq : HasDerivAt q (qa*h+qp*v) t) :
    deriv (fun s => rhoQ b (u s) (q s)) t+
      3*h*(rhoQ b (u t) (q t)+pressureQ b (u t) (q t) qa)=
      b*(u t)^2*qp*v := by
  rw [(quantum_energy_derivative u q t b h qa qp v hu hq).deriv]
  dsimp [rhoQ, pressureQ]
  ring

theorem scalar_exchange (v : ℝ → ℝ) (t b u h qp : ℝ) (hb : b ≠ 0)
    (hv : HasDerivAt v (-3*h*v t-b^2*u^2*qp) t) :
    deriv (fun s => rhoF b (v s)) t+3*h*(2*rhoF b (v t))=
      -b*u^2*qp*v t := by
  have hd := (hv.pow 2).div_const (2*b)
  change HasDerivAt (fun s => rhoF b (v s)) _ t at hd
  rw [hd.deriv]
  dsimp [rhoF]
  field_simp
  <;> ring

theorem total_exchange (u q v : ℝ → ℝ) (t b h qa qp : ℝ) (hb : b ≠ 0)
    (hu : HasDerivAt u (-3*h*u t) t)
    (hq : HasDerivAt q (qa*h+qp*v t) t)
    (hv : HasDerivAt v (-3*h*v t-b^2*(u t)^2*qp) t) :
    deriv (fun s => rhoF b (v s)) t + deriv (fun s => rhoQ b (u s) (q s)) t+
      3*h*(2*rhoF b (v t)+rhoQ b (u t) (q t)+pressureQ b (u t) (q t) qa)=0 := by
  have h1 := scalar_exchange v t b (u t) h qp hb hv
  have h2 := quantum_exchange u q t b h qa qp (v t) hu hq
  linarith

theorem friedmann (b u sa sp q : ℝ)
    (hj : (-sa^2+sp^2)/2+q=0) :
    (-b*u*sa)^2=(b*u*sp)^2+2*b^2*u^2*q := by
  linear_combination -2*b^2*u^2*hj

theorem friedmann_source (b u h v q : ℝ) (hb : b ≠ 0)
    (hf : h^2=v^2+2*b^2*u^2*q) :
    h^2=2*b*(rhoF b v+rhoQ b u q) := by
  dsimp [rhoF, rhoQ]
  field_simp
  nlinarith [hf]

theorem guidance_acceleration (u sa : ℝ → ℝ) (t b h v saa sap sp qa : ℝ)
    (hu : HasDerivAt u (-3*h*u t) t)
    (hs : HasDerivAt sa (saa*h+sap*v) t)
    (hh : h= -b*u t*sa t) (hv : v=b*u t*sp)
    (hj : -sa t*saa+sp*sap+qa=0) :
    HasDerivAt (fun s => -b*u s*sa s)
      (-3*h^2+b^2*(u t)^2*qa) t := by
  have hd := (hu.const_mul (-b)).mul hs
  convert hd using 1
  rw [hh, hv]
  have he : qa=sa t*saa-sp*sap := by linarith
  rw [he]
  ring

theorem scalar_guidance_acceleration (u sp : ℝ → ℝ) (t b h v spa spp sa qp : ℝ)
    (hu : HasDerivAt u (-3*h*u t) t)
    (hs : HasDerivAt sp (spa*h+spp*v) t)
    (hh : h= -b*u t*sa) (hv : v=b*u t*sp t)
    (hj : -sa*spa+sp t*spp+qp=0) :
    HasDerivAt (fun s => b*u s*sp s)
      (-3*h*v-b^2*(u t)^2*qp) t := by
  have hd := (hu.const_mul b).mul hs
  convert hd using 1
  rw [hh, hv]
  have he : qp=sa*spa-sp t*spp := by linarith
  rw [he]
  ring

theorem friedmann_constraint_propagation (h v u q : ℝ → ℝ) (t b qa qp : ℝ)
    (hh : HasDerivAt h (-3*(h t)^2+b^2*(u t)^2*qa) t)
    (hv : HasDerivAt v (-3*h t*v t-b^2*(u t)^2*qp) t)
    (hu : HasDerivAt u (-3*h t*u t) t)
    (hq : HasDerivAt q (qa*h t+qp*v t) t) :
    HasDerivAt (fun s => (h s)^2-(v s)^2-2*b^2*(u s)^2*q s)
      (-6*h t*((h t)^2-(v t)^2-2*b^2*(u t)^2*q t)) t := by
  convert ((hh.pow 2).sub (hv.pow 2)).sub (((hu.pow 2).const_mul (2*b^2)).mul hq) using 1
  simp only [Pi.pow_apply]
  ring

theorem weighted_constraint_conserved (a c : ℝ → ℝ) (t h : ℝ)
    (ha : HasDerivAt a h t) (hc : HasDerivAt c (-6*h*c t) t) :
    HasDerivAt (fun s => Real.exp (6*a s)*c s) 0 t := by
  convert ((ha.const_mul 6).exp).mul hc using 1
  ring

theorem raychaudhuri (b u h v q qa hd : ℝ) (hb : b ≠ 0)
    (hf : h^2=v^2+2*b^2*u^2*q)
    (ha : hd= -3*h^2+b^2*u^2*qa) :
    hd= -3*b*(2*rhoF b v+rhoQ b u q+pressureQ b u q qa) := by
  dsimp [rhoF, rhoQ, pressureQ]
  field_simp
  nlinarith [hf, ha]

theorem bounce_requires_negative_quantum_density (b u v q : ℝ)
    (hb : 0 < b) (hu : 0 < u) (hv : v ≠ 0)
    (hf : 0=v^2+2*b^2*u^2*q) : q < 0 := by
  have hv2 : 0 < v^2 := sq_pos_of_ne_zero hv
  by_contra h
  have hq : 0 ≤ q := le_of_not_gt h
  have hprod : 0 ≤ 2*b^2*u^2*q := by positivity
  linarith

end
end UnifiedWard
