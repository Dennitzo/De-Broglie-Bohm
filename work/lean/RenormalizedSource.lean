import Mathlib

/-! Finite FLRW source identities for specified trace coefficients.
The source is defined explicitly; a field-theoretic construction is not
part of this module. Time derivatives are genuine HasDerivAt premises. -/

namespace RenormalizedSource

noncomputable section

def curvature (H H1 : ℝ) : ℝ := 6*(H1+2*H^2)
def curvatureVelocity (H H1 H2 : ℝ) : ℝ := 6*(H2+4*H*H1)
def curvatureAcceleration (H H1 H2 H3 : ℝ) : ℝ := 6*(H3+4*H1^2+4*H*H2)
def boxCurvature (H H1 H2 H3 : ℝ) : ℝ :=
  -curvatureAcceleration H H1 H2 H3-3*H*curvatureVelocity H H1 H2

def traceSource (A b H H1 H2 H3 : ℝ) : ℝ :=
  -4*A*(H^4+H^2*H1)+b*boxCurvature H H1 H2 H3

def density (A b C a H H1 H2 : ℝ) : ℝ :=
  A*H^4+b*(H*curvatureVelocity H H1 H2+H^2*curvature H H1-
    (curvature H H1)^2/12)+C/a^4

def pressure (A b C a H H1 H2 H3 : ℝ) : ℝ :=
  (traceSource A b H H1 H2 H3+density A b C a H H1 H2)/3

def densityRate (A b C a H H1 H2 H3 : ℝ) : ℝ :=
  4*A*H^3*H1+b*(H1*curvatureVelocity H H1 H2+
    H*curvatureAcceleration H H1 H2 H3+2*H*H1*curvature H H1+
    H^2*curvatureVelocity H H1 H2-
    curvature H H1*curvatureVelocity H H1 H2/6)-4*H*C/a^4

theorem curvature_hasDerivAt (H H1 : ℝ → ℝ) (t dH1 : ℝ)
    (hH : HasDerivAt H (H1 t) t) (hH1 : HasDerivAt H1 dH1 t) :
    HasDerivAt (fun u => curvature (H u) (H1 u))
      (curvatureVelocity (H t) (H1 t) dH1) t := by
  convert (hH1.add ((hH.pow 2).const_mul 2)).const_mul 6 using 1 <;>
    dsimp [curvature, curvatureVelocity] <;> ring

theorem curvature_velocity_hasDerivAt (H H1 H2 : ℝ → ℝ) (t dH2 : ℝ)
    (hH : HasDerivAt H (H1 t) t) (hH1 : HasDerivAt H1 (H2 t) t)
    (hH2 : HasDerivAt H2 dH2 t) :
    HasDerivAt (fun u => curvatureVelocity (H u) (H1 u) (H2 u))
      (curvatureAcceleration (H t) (H1 t) (H2 t) dH2) t := by
  convert (hH2.add ((hH.mul hH1).const_mul 4)).const_mul 6 using 1 <;>
    dsimp [curvatureVelocity, curvatureAcceleration] <;> ring

theorem radiation_density_hasDerivAt (a : ℝ → ℝ) (t H C : ℝ)
    (ha : HasDerivAt a (H*a t) t) (ha0 : a t ≠ 0) :
    HasDerivAt (fun u => C/(a u)^4) (-4*H*C/(a t)^4) t := by
  have hn : (a t)^4 ≠ 0 := pow_ne_zero 4 ha0
  convert (hasDerivAt_const t C).div (ha.pow 4) hn using 1 <;>
    dsimp <;> field_simp <;> ring

theorem density_hasDerivAt (a H H1 H2 : ℝ → ℝ) (t A b C H3 : ℝ)
    (ha : HasDerivAt a (H t*a t) t) (ha0 : a t ≠ 0)
    (hH : HasDerivAt H (H1 t) t) (hH1 : HasDerivAt H1 (H2 t) t)
    (hH2 : HasDerivAt H2 H3 t) :
    HasDerivAt (fun u => density A b C (a u) (H u) (H1 u) (H2 u))
      (densityRate A b C (a t) (H t) (H1 t) (H2 t) H3) t := by
  have hR := curvature_hasDerivAt H H1 t (H2 t) hH hH1
  have hRd := curvature_velocity_hasDerivAt H H1 H2 t H3 hH hH1 hH2
  have hlocal := (((hH.mul hRd).add ((hH.pow 2).mul hR)).sub
    ((hR.pow 2).div_const 12)).const_mul b
  have hall := (((hH.pow 4).const_mul A).add hlocal).add
    (radiation_density_hasDerivAt a t (H t) C ha ha0)
  convert hall using 1 <;> dsimp [density, densityRate] <;> ring

theorem source_trace_identity (A b C a H H1 H2 H3 : ℝ) :
    -density A b C a H H1 H2+3*pressure A b C a H H1 H2 H3=
      traceSource A b H H1 H2 H3 := by
  unfold pressure
  ring

theorem source_continuity_algebra (A b C a H H1 H2 H3 : ℝ) :
    densityRate A b C a H H1 H2 H3+
      3*H*(density A b C a H H1 H2+pressure A b C a H H1 H2 H3)=0 := by
  unfold pressure densityRate density traceSource boxCurvature
    curvature curvatureVelocity curvatureAcceleration
  ring

theorem source_conserved_hasDerivAt (a H H1 H2 : ℝ → ℝ) (t A b C H3 : ℝ)
    (ha : HasDerivAt a (H t*a t) t) (ha0 : a t ≠ 0)
    (hH : HasDerivAt H (H1 t) t) (hH1 : HasDerivAt H1 (H2 t) t)
    (hH2 : HasDerivAt H2 H3 t) :
    HasDerivAt (fun u => density A b C (a u) (H u) (H1 u) (H2 u))
      (-3*H t*(density A b C (a t) (H t) (H1 t) (H2 t)+
        pressure A b C (a t) (H t) (H1 t) (H2 t) H3)) t := by
  have hh := density_hasDerivAt a H H1 H2 t A b C H3 ha ha0 hH hH1 hH2
  have he := source_continuity_algebra A b C (a t) (H t) (H1 t) (H2 t) H3
  convert hh using 1 <;> linarith

theorem source_continuity (a H H1 H2 : ℝ → ℝ) (t A b C H3 : ℝ)
    (ha : HasDerivAt a (H t*a t) t) (ha0 : a t ≠ 0)
    (hH : HasDerivAt H (H1 t) t) (hH1 : HasDerivAt H1 (H2 t) t)
    (hH2 : HasDerivAt H2 H3 t) :
    deriv (fun u => density A b C (a u) (H u) (H1 u) (H2 u)) t+
      3*H t*(density A b C (a t) (H t) (H1 t) (H2 t)+
        pressure A b C (a t) (H t) (H1 t) (H2 t) H3)=0 := by
  rw [(source_conserved_hasDerivAt a H H1 H2 t A b C H3 ha ha0 hH hH1 hH2).deriv]
  ring

theorem zero_counterterm_density (A C a H H1 H2 : ℝ) :
    density A 0 C a H H1 H2=A*H^4+C/a^4 := by
  unfold density
  ring

theorem zero_counterterm_density_hasDerivAt (a H : ℝ → ℝ) (t A C H1 : ℝ)
    (ha : HasDerivAt a (H t*a t) t) (ha0 : a t ≠ 0)
    (hH : HasDerivAt H H1 t) :
    HasDerivAt (fun u => A*(H u)^4+C/(a u)^4)
      (4*A*(H t)^3*H1-4*H t*C/(a t)^4) t := by
  convert ((hH.pow 4).const_mul A).add
    (radiation_density_hasDerivAt a t (H t) C ha ha0) using 1 <;> dsimp <;> ring

theorem zero_counterterm_null_energy (A C a H H1 H2 H3 : ℝ) :
    density A 0 C a H H1 H2+pressure A 0 C a H H1 H2 H3=
      (-4*A*H^2*H1+4*C/a^4)/3 := by
  unfold pressure density traceSource
  ring

theorem de_sitter_density (A b a H : ℝ) :
    density A b 0 a H 0 0=A*H^4 := by
  unfold density curvature curvatureVelocity
  ring

theorem de_sitter_equation_of_state (A b a H : ℝ) :
    pressure A b 0 a H 0 0 0 = -density A b 0 a H 0 0 := by
  unfold pressure density traceSource boxCurvature curvature
    curvatureVelocity curvatureAcceleration
  ring

theorem de_sitter_density_positive (A b a H : ℝ) (hA : 0 < A) (hH : H ≠ 0) :
    0 < density A b 0 a H 0 0 := by
  rw [de_sitter_density]
  have hh : 0 < H^4 := by positivity
  exact mul_pos hA hh

def constraintResidual (beta lam D a h r s : ℝ) : ℝ :=
  h^2-h^4-beta*(h*s+h^2*r-r^2/12)-lam-D/a^4

theorem constraint_residual_hasDerivAt (a h r s : ℝ → ℝ) (t beta lam D : ℝ)
    (hb : beta ≠ 0) (ha0 : a t ≠ 0)
    (ha : HasDerivAt a (h t*a t) t)
    (hh : HasDerivAt h (r t/6-2*(h t)^2) t)
    (hr : HasDerivAt r (s t) t)
    (hs : HasDerivAt s
      ((r t-12*lam+12*(h t)^4-2*(h t)^2*r t)/(3*beta)-3*h t*s t) t) :
    HasDerivAt (fun u => constraintResidual beta lam D (a u) (h u) (r u) (s u))
      (-4*h t*constraintResidual beta lam D (a t) (h t) (r t) (s t)) t := by
  have hlocal := (((hh.mul hs).add ((hh.pow 2).mul hr)).sub
    ((hr.pow 2).div_const 12)).const_mul beta
  have hall := (((((hh.pow 2).sub (hh.pow 4)).sub hlocal).sub_const lam).sub
    (radiation_density_hasDerivAt a t (h t) D ha ha0))
  convert hall using 1 <;> dsimp [constraintResidual] <;> field_simp <;> ring

theorem scaled_constraint_hasDerivAt (a E : ℝ → ℝ) (t h : ℝ)
    (ha : HasDerivAt a (h*a t) t) (hE : HasDerivAt E (-4*h*E t) t) :
    HasDerivAt (fun u => E u*(a u)^4) 0 t := by
  convert hE.mul (ha.pow 4) using 1 <;> dsimp <;> ring

theorem constraint_zero_preserved_on_interval (a E h : ℝ → ℝ) (l r t t0 : ℝ)
    (ht : t ∈ Set.Ioo l r) (ht0 : t0 ∈ Set.Ioo l r) (ha0 : a t ≠ 0)
    (ha : ∀ u ∈ Set.Ioo l r, HasDerivAt a (h u*a u) u)
    (hE : ∀ u ∈ Set.Ioo l r, HasDerivAt E (-4*h u*E u) u)
    (hinit : E t0=0) : E t=0 := by
  have hz : ∀ u ∈ Set.Ioo l r, HasDerivAt (fun z => E z*(a z)^4) 0 u := by
    intro u hu
    exact scaled_constraint_hasDerivAt a E u (h u) (ha u hu) (hE u hu)
  have hc : E t*(a t)^4=E t0*(a t0)^4 :=
    isOpen_Ioo.is_const_of_deriv_eq_zero (convex_Ioo l r).isPreconnected
      (fun u hu => (hz u hu).differentiableAt.differentiableWithinAt)
      (fun u hu => (hz u hu).deriv) ht ht0
  rw [hinit, zero_mul] at hc
  exact (mul_eq_zero.mp hc).resolve_right (pow_ne_zero 4 ha0)

end

end RenormalizedSource
