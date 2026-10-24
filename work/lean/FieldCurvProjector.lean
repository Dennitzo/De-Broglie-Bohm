import Mathlib

/-! Only the finite real two-dimensional projection algebra of a smooth mode insertion. -/
namespace FieldCurvProjector
def occupied (A B : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![B^2,-A*B; -A*B,A^2]
def inserted (A B : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![A^2,A*B; A*B,B^2]

theorem occupied_idempotent (A B : ℝ) (h : A^2+B^2=1) :
    occupied A B * occupied A B=occupied A B := by
  have hA := congrArg (fun x : ℝ => A^2*x) h
  have hB := congrArg (fun x : ℝ => B^2*x) h
  have hAB := congrArg (fun x : ℝ => A*B*x) h
  dsimp at hA hB hAB
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [occupied,Matrix.mul_apply,Fin.sum_univ_two] <;> nlinarith

theorem inserted_idempotent (A B : ℝ) (h : A^2+B^2=1) :
    inserted A B * inserted A B=inserted A B := by
  have hA := congrArg (fun x : ℝ => A^2*x) h
  have hB := congrArg (fun x : ℝ => B^2*x) h
  have hAB := congrArg (fun x : ℝ => A*B*x) h
  dsimp at hA hB hAB
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [inserted,Matrix.mul_apply,Fin.sum_univ_two] <;> nlinarith

theorem orthogonal_complement (A B : ℝ) (h : A^2+B^2=1) :
    occupied A B * inserted A B=0 ∧ occupied A B+inserted A B=1 := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [occupied,inserted,Matrix.mul_apply,Fin.sum_univ_two] <;> ring
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [occupied,inserted,Matrix.one_apply] <;> nlinarith

end FieldCurvProjector
