import Mathlib

/-! Finite dimensional homogeneous spinor dynamics and the reduced anisotropic source. -/
namespace DiracBackMatter
noncomputable section

def normSq (a b c d : ℝ) : ℝ := a^2+b^2+c^2+d^2
def sOne (a b c d : ℝ) : ℝ := 2*(a*c+b*d)
def sTwo (a b c d : ℝ) : ℝ := 2*(a*d-b*c)
def sThree (a b c d : ℝ) : ℝ := a^2+b^2-c^2-d^2
def energy (mu kap B S1 S3 : ℝ) : ℝ := mu*S3+kap*S1/B
def density (mu kap A B S1 S3 : ℝ) : ℝ := energy mu kap B S1 S3/(A*B^2)
def angularPressure (kap A B S1 : ℝ) : ℝ := kap*S1/(2*A*B^3)

theorem spinor_bloch_norm (a b c d : ℝ) :
    (sOne a b c d)^2+(sTwo a b c d)^2+(sThree a b c d)^2=(normSq a b c d)^2 := by
  unfold sOne sTwo sThree normSq
  ring

theorem spinor_norm_hasDerivAt (a b c d : ℝ → ℝ) (t mu ell : ℝ)
    (ha : HasDerivAt a (mu*b t+ell*d t) t)
    (hb : HasDerivAt b (-mu*a t-ell*c t) t)
    (hc : HasDerivAt c (ell*b t-mu*d t) t)
    (hd : HasDerivAt d (-ell*a t+mu*c t) t) :
    HasDerivAt (fun x => normSq (a x) (b x) (c x) (d x)) 0 t := by
  convert (((ha.pow 2).add (hb.pow 2)).add (hc.pow 2)).add (hd.pow 2) using 1 <;>
    dsimp [normSq] <;> ring

theorem spinor_sOne_hasDerivAt (a b c d : ℝ → ℝ) (t mu ell : ℝ)
    (ha : HasDerivAt a (mu*b t+ell*d t) t)
    (hb : HasDerivAt b (-mu*a t-ell*c t) t)
    (hc : HasDerivAt c (ell*b t-mu*d t) t)
    (hd : HasDerivAt d (-ell*a t+mu*c t) t) :
    HasDerivAt (fun x => sOne (a x) (b x) (c x) (d x))
      (-2*mu*sTwo (a t) (b t) (c t) (d t)) t := by
  convert ((ha.mul hc).add (hb.mul hd)).const_mul 2 using 1 <;>
    dsimp [sOne,sTwo] <;> ring

theorem spinor_sTwo_hasDerivAt (a b c d : ℝ → ℝ) (t mu ell : ℝ)
    (ha : HasDerivAt a (mu*b t+ell*d t) t)
    (hb : HasDerivAt b (-mu*a t-ell*c t) t)
    (hc : HasDerivAt c (ell*b t-mu*d t) t)
    (hd : HasDerivAt d (-ell*a t+mu*c t) t) :
    HasDerivAt (fun x => sTwo (a x) (b x) (c x) (d x))
      (2*mu*sOne (a t) (b t) (c t) (d t)-2*ell*sThree (a t) (b t) (c t) (d t)) t := by
  convert ((ha.mul hd).sub (hb.mul hc)).const_mul 2 using 1 <;>
    dsimp [sOne,sTwo,sThree] <;> ring

theorem spinor_sThree_hasDerivAt (a b c d : ℝ → ℝ) (t mu ell : ℝ)
    (ha : HasDerivAt a (mu*b t+ell*d t) t)
    (hb : HasDerivAt b (-mu*a t-ell*c t) t)
    (hc : HasDerivAt c (ell*b t-mu*d t) t)
    (hd : HasDerivAt d (-ell*a t+mu*c t) t) :
    HasDerivAt (fun x => sThree (a x) (b x) (c x) (d x))
      (2*ell*sTwo (a t) (b t) (c t) (d t)) t := by
  convert (((ha.pow 2).add (hb.pow 2)).sub (hc.pow 2)).sub (hd.pow 2) using 1 <;>
    dsimp [sTwo,sThree] <;> ring

theorem bloch_norm_hasDerivAt (S1 S2 S3 : ℝ → ℝ) (t mu ell : ℝ)
    (h1 : HasDerivAt S1 (-2*mu*S2 t) t)
    (h2 : HasDerivAt S2 (2*mu*S1 t-2*ell*S3 t) t)
    (h3 : HasDerivAt S3 (2*ell*S2 t) t) :
    HasDerivAt (fun x => (S1 x)^2+(S2 x)^2+(S3 x)^2) 0 t := by
  convert ((h1.pow 2).add (h2.pow 2)).add (h3.pow 2) using 1 <;> ring

theorem mode_energy_hasDerivAt (B S1 S3 : ℝ → ℝ) (t mu kap b S2 : ℝ)
    (hB : HasDerivAt B b t) (hBn : B t ≠ 0)
    (h1 : HasDerivAt S1 (-2*mu*S2) t)
    (h3 : HasDerivAt S3 (2*kap/B t*S2) t) :
    HasDerivAt (fun x => energy mu kap (B x) (S1 x) (S3 x))
      (-kap*S1 t*b/(B t)^2) t := by
  convert (h3.const_mul mu).add ((h1.const_mul kap).div hB hBn) using 1 <;>
    dsimp [energy] <;> field_simp <;> ring

theorem density_hasDerivAt (A B E : ℝ → ℝ) (t a b kap S1 : ℝ)
    (hA : HasDerivAt A a t) (hB : HasDerivAt B b t)
    (hE : HasDerivAt E (-kap*S1*b/(B t)^2) t)
    (hAn : A t ≠ 0) (hBn : B t ≠ 0) :
    HasDerivAt (fun x => E x/(A x*(B x)^2))
      (-(a/A t)*(E t/(A t*(B t)^2))-
        2*(b/B t)*(E t/(A t*(B t)^2)+angularPressure kap (A t) (B t) S1)) t := by
  convert hE.div (hA.mul (hB.pow 2)) (by exact mul_ne_zero hAn (pow_ne_zero 2 hBn)) using 1 <;>
    dsimp [angularPressure] <;> field_simp <;> ring

theorem source_trace (mu kap A B S1 S3 : ℝ) (hA : A ≠ 0) (hB : B ≠ 0) :
    -density mu kap A B S1 S3+2*angularPressure kap A B S1 = -mu*S3/(A*B^2) := by
  unfold density energy angularPressure
  field_simp
  ring

theorem longitudinal_hamiltonian_anticommutator (mu ell : ℂ) :
    (!![0,-Complex.I;Complex.I,0] : Matrix (Fin 2) (Fin 2) ℂ) * !![mu,ell;ell,-mu] +
      !![mu,ell;ell,-mu] * !![0,-Complex.I;Complex.I,0] = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply,Fin.sum_univ_two] <;> ring

theorem single_mode_tangent_norm (S1 S2 S3 theta : ℝ)
    (h : S1^2+S2^2+S3^2=1) :
    -1+S2^2+S1^2*(Real.sin theta)^2 = -S3^2-S1^2*(Real.cos theta)^2 := by
  nlinarith [mul_nonneg (sq_nonneg S1) (sq_nonneg (Real.cos theta)),
    show S1^2*((Real.sin theta)^2+(Real.cos theta)^2)=S1^2 by
      rw [Real.sin_sq_add_cos_sq]; ring]

theorem single_mode_causal (S1 S2 S3 theta : ℝ)
    (h : S1^2+S2^2+S3^2=1) : -1+S2^2+S1^2*(Real.sin theta)^2 ≤ 0 := by
  rw [single_mode_tangent_norm S1 S2 S3 theta h]
  nlinarith [sq_nonneg S3,mul_nonneg (sq_nonneg S1) (sq_nonneg (Real.cos theta))]

theorem single_mode_timelike_if_S3_nonzero (S1 S2 S3 theta : ℝ)
    (h : S1^2+S2^2+S3^2=1) (h3 : S3 ≠ 0) :
    -1+S2^2+S1^2*(Real.sin theta)^2 < 0 := by
  rw [single_mode_tangent_norm S1 S2 S3 theta h]
  nlinarith [sq_pos_of_ne_zero h3,mul_nonneg (sq_nonneg S1) (sq_nonneg (Real.cos theta))]

end
end DiracBackMatter
