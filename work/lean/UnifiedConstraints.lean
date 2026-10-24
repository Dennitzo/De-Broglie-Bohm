import Mathlib

/-! Polar constraint identities and weak Hamilton-Jacobi graph closure. -/
namespace UnifiedConstraints
noncomputable section

theorem polar_hamilton_jacobi (r raa rpp sa sp : ℝ) (hr : r ≠ 0)
    (hreal : raa-rpp-r*(sa^2-sp^2)=0) :
    (-sa^2+sp^2)/2+(raa-rpp)/(2*r)=0 := by
  field_simp
  nlinarith [hreal]

theorem polar_continuity (r ra rp sa sp saa spp : ℝ)
    (him : 2*ra*sa+r*saa-2*rp*sp-r*spp=0) :
    -(2*r*ra*sa+r^2*saa)+(2*r*rp*sp+r^2*spp)=0 := by
  nlinarith [mul_eq_zero_of_left him r]

theorem continuity_derivatives (r sa sp : ℝ → ℝ → ℝ)
    (a p ra rp saa spp : ℝ)
    (hra : HasDerivAt (fun x => r x p) ra a)
    (hrp : HasDerivAt (r a) rp p)
    (hsa : HasDerivAt (fun x => sa x p) saa a)
    (hsp : HasDerivAt (sp a) spp p)
    (him : 2*ra*sa a p+r a p*saa-2*rp*sp a p-r a p*spp=0) :
    -deriv (fun x => (r x p)^2*sa x p) a+
      deriv (fun x => (r a x)^2*sp a x) p=0 := by
  have hA := (hra.pow 2).mul hsa
  have hP := (hrp.pow 2).mul hsp
  change HasDerivAt (fun x => (r x p)^2*sa x p) _ a at hA
  change HasDerivAt (fun x => (r a x)^2*sp a x) _ p at hP
  rw [hA.deriv, hP.deriv]
  have hc := polar_continuity (r a p) ra rp (sa a p) (sp a p) saa spp him
  simp only [Pi.pow_apply, Nat.reduceSub, pow_one, Nat.cast_ofNat]
  exact hc

theorem quantum_force_tangent (sa sp saa sap sppa spp qaa qpp : ℝ)
    (hshared : sap=sppa)
    (ha : -sa*saa+sp*sppa+qaa=0)
    (hp : -sa*sap+sp*spp+qpp=0) :
    (saa*(-sa)+sap*sp= -qaa) ∧
      (sppa*(-sa)+spp*sp= -qpp) := by
  subst sppa
  constructor <;> nlinarith

theorem graph_weak_poisson_closure
    (hq1 hq2 hp1 hp2 kq1 kq2 kp1 kp2 s11 s12 s21 s22 : ℝ)
    (hschwarz : s12=s21)
    (hh1 : hq1+hp1*s11+hp2*s21=0)
    (hh2 : hq2+hp1*s12+hp2*s22=0)
    (hk1 : kq1+kp1*s11+kp2*s21=0)
    (hk2 : kq2+kp1*s12+kp2*s22=0) :
    hq1*kp1+hq2*kp2-hp1*kq1-hp2*kq2=0 := by
  have e1 : hq1= -hp1*s11-hp2*s21 := by linarith
  have e2 : hq2= -hp1*s12-hp2*s22 := by linarith
  have e3 : kq1= -kp1*s11-kp2*s21 := by linarith
  have e4 : kq2= -kp1*s12-kp2*s22 := by linarith
  rw [e1,e2,e3,e4,hschwarz]
  ring

theorem weak_graph_counterexample_bracket (q p : ℝ) :
    deriv (fun _ : ℝ => p) q * deriv (fun y : ℝ => q*y) p-
      deriv (fun y : ℝ => y) p * deriv (fun x : ℝ => x*p) q= -p := by
  have h1 := (hasDerivAt_id p).const_mul q
  have h2 := (hasDerivAt_id q).mul_const p
  have h3 := hasDerivAt_id p
  change HasDerivAt (fun y : ℝ => q*y) _ p at h1
  change HasDerivAt (fun x : ℝ => x*p) _ q at h2
  change HasDerivAt (fun y : ℝ => y) _ p at h3
  rw [(hasDerivAt_const q p).deriv, h1.deriv, h2.deriv, h3.deriv]
  ring

theorem weak_graph_counterexample_off_graph :
    deriv (fun _ : ℝ => (1:ℝ)) 0 * deriv (fun y : ℝ => 0*y) 1-
      deriv (fun y : ℝ => y) 1 * deriv (fun x : ℝ => x*1) 0 ≠ 0 := by
  rw [weak_graph_counterexample_bracket]
  norm_num

end
end UnifiedConstraints
