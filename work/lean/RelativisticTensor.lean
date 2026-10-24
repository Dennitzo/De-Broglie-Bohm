import Mathlib

/-! Separate tensor-factor actions for free Fourier modes and the normalized
Schmidt coefficient matrix. Spatial mode orthonormality is an analytic
input when these coefficients are used for wave packets. -/

namespace RelativisticTensor

noncomputable section

abbrev TwoMatrix := Matrix (Fin 2) (Fin 2) ℂ

def actLeft (H P : TwoMatrix) : TwoMatrix :=
  fun i j => ∑ k : Fin 2, H i k * P k j

def actRight (K P : TwoMatrix) : TwoMatrix :=
  fun i j => ∑ k : Fin 2, K j k * P i k

theorem separate_tensor_actions_commute (H K P : TwoMatrix) :
    actLeft H (actRight K P) = actRight K (actLeft H P) := by
  funext i j
  simp only [actLeft, actRight, Fin.sum_univ_two]
  ring

def diracMode (p m : ℝ) : TwoMatrix := !![(p : ℂ), (m : ℂ); (m : ℂ), -(p : ℂ)]

theorem free_dirac_modes_commute (p q m : ℝ) (P : TwoMatrix) :
    actLeft (diracMode p m) (actRight (diracMode q m) P) =
      actRight (diracMode q m) (actLeft (diracMode p m) P) := by
  exact separate_tensor_actions_commute _ _ _

theorem dirac_mode_square (p m : ℝ) :
    (diracMode p m) * (diracMode p m) =
      !![((p^2+m^2 : ℝ) : ℂ), 0; 0, ((p^2+m^2 : ℝ) : ℂ)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, diracMode, Fin.sum_univ_two] <;> ring

def bellCoefficients (r : ℝ) : TwoMatrix := !![(r : ℂ), 0; 0, Complex.I*(r : ℂ)]

def reduced (P : TwoMatrix) : TwoMatrix :=
  fun i j => ∑ k : Fin 2, P i k * star (P j k)

def coefficientDet (P : TwoMatrix) : ℂ := P 0 0 * P 1 1 - P 0 1 * P 1 0

theorem product_coefficients_zero_determinant (u v : Fin 2 → ℂ) :
    coefficientDet (fun i j => u i * v j) = 0 := by
  unfold coefficientDet
  ring

theorem bell_determinant (r : ℝ) (hr : r*r=(1:ℝ)/2) :
    coefficientDet (bellCoefficients r) = Complex.I/2 := by
  have hc : (r:ℂ)*(r:ℂ)=(1:ℂ)/2 := by
    simpa using congrArg (fun x : ℝ => (x : ℂ)) hr
  simp only [coefficientDet, bellCoefficients, Matrix.of_apply]
  change (r:ℂ)*(Complex.I*(r:ℂ))-0*0=Complex.I/2
  linear_combination Complex.I * hc

theorem bell_entangled (r : ℝ) (hr : r*r=(1:ℝ)/2) :
    ¬ ∃ u v : Fin 2 → ℂ, bellCoefficients r = fun i j => u i * v j := by
  rintro ⟨u,v,h⟩
  have hp := product_coefficients_zero_determinant u v
  rw [← h, bell_determinant r hr] at hp
  norm_num [Complex.ext_iff] at hp

theorem bell_reduced_density (r : ℝ) (hr : r*r=(1:ℝ)/2) :
    reduced (bellCoefficients r) = !![(1:ℂ)/2, 0; 0, (1:ℂ)/2] := by
  have hc : (r:ℂ)*(r:ℂ)=(1:ℂ)/2 := by
    simpa using congrArg (fun x : ℝ => (x : ℂ)) hr
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [reduced, bellCoefficients, Fin.sum_univ_two, hc]
  calc
    _ = (r:ℂ)*(r:ℂ) := by linear_combination -(r:ℂ)^2 * Complex.I_sq
    _ = _ := by simpa using hc

theorem bell_reduced_purity (r : ℝ) (hr : r*r=(1:ℝ)/2) :
    let R := reduced (bellCoefficients r)
    ∑ i : Fin 2, ∑ j : Fin 2, R i j * R j i = (1:ℂ)/2 := by
  dsimp only
  rw [bell_reduced_density r hr]
  norm_num [Fin.sum_univ_two]

theorem bell_reduced_trace (r : ℝ) (hr : r*r=(1:ℝ)/2) :
    ∑ i : Fin 2, reduced (bellCoefficients r) i i = 1 := by
  rw [bell_reduced_density r hr]
  norm_num [Fin.sum_univ_two]

theorem inverse_sqrt_two_normalizes :
    (1 / Real.sqrt 2) * (1 / Real.sqrt 2) = (1:ℝ)/2 := by
  have hs : (Real.sqrt 2)^2 = (2:ℝ) := Real.sq_sqrt (by norm_num)
  have hne : Real.sqrt 2 ≠ 0 := by positivity
  field_simp
  nlinarith

end

end RelativisticTensor
