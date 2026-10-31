import Mathlib

/-!
Generischer Kontinuitaetskern fuer endliche Feldzahl.

Es sei I eine endliche Indexmenge (n+1 Koordinaten), T^{AB} eine symmetrische
kontragierte Metrik-Masse, u und phi die Real- und Imaginaerteile des Zustands.

Hauptidentitaet:  sum_{A,B} T^{AB} (u_A phi_B - phi_A u_B) = 0.
Beweis:因 T symmetrisch und die Klammer antisymmetrisch ist, verschwindet die Doppelsumme
durch Achsenvertauschung. Dies ist dimensionsfrei für jede endliche Menge.
-/

open BigOperators Finset

namespace GenericCurrent

/-- Die Doppelsumme einer symmetrischen gegen eine antisymmetrische Größe ist Null. -/
theorem sym_contract_antisym {I : Type} [Fintype I] (T : I → I → ℝ) 
    (hu : ∀ i j, T i j = T j i) (A : I → I → ℝ) 
    (hA : ∀ i j, A i j = - A j i) : 
    ∑ i : I, ∑ j : I, T i j * A i j = 0 := by
  have h1 : ∑ i : I, ∑ j : I, T i j * A i j = ∑ j : I, ∑ i : I, T j i * A j i := by
    rw [sum_comm]
  have h2 : ∑ j : I, ∑ i : I, T j i * A j i = ∑ j : I, ∑ i : I, T i j * (-A i j) := by
    apply sum_congr rfl
    intro j _
    apply sum_congr rfl
    intro i _
    rw [hu, hA]
  have h3 : ∑ j : I, ∑ i : I, T i j * (-A i j) = -∑ j : I, ∑ i : I, T i j * A i j := by
    have h4 : ∑ j : I, ∑ i : I, T i j * (-A i j) = ∑ j : I, ∑ i : I, -(T i j * A i j) := by
      apply sum_congr rfl
      intro j _
      apply sum_congr rfl
      intro i _
      ring
    have h5 : ∑ j : I, ∑ i : I, -(T i j * A i j) = -∑ j : I, ∑ i : I, T i j * A i j := by
      simp
    rw [h4, h5]
  have h6 : ∑ i : I, ∑ j : I, T i j * A i j = -∑ i : I, ∑ j : I, T i j * A i j := by
    calc ∑ i : I, ∑ j : I, T i j * A i j
      _ = ∑ j : I, ∑ i : I, T j i * A j i := h1
      _ = ∑ j : I, ∑ i : I, T i j * (-A i j) := h2
      _ = -∑ j : I, ∑ i : I, T i j * A i j := h3
      _ = -∑ i : I, ∑ j : I, T i j * A i j := by rw [← sum_comm]
  have h7 : ∑ i : I, ∑ j : I, T i j * A i j + ∑ i : I, ∑ j : I, T i j * A i j = 0 := by
    linarith
  nlinarith

/-- Anwendung: Die Kontraktion einer symmetrischen Metrik-Masse $T^{AB}$ mit der 
antisymmetrischen Größe $(u_A \phi_B - \phi_A u_B)$ ergibt Null. -/
theorem sym_contract_antisym_apply {I : Type} [Fintype I] (T : I → I → ℝ) (hT : ∀ i j, T i j = T j i)
    (u v : I → ℝ) :
    ∑ i : I, ∑ j : I, T i j * (u i * v j - v i * u j) = 0 := by
  have hA : ∀ i j, (fun i j => u i * v j - v i * u j) i j = -(fun i j => u i * v j - v i * u j) j i := by
    intro i j
    ring
  simpa using sym_contract_antisym T hT (fun i j => u i * v j - v i * u j) hA

end GenericCurrent

#print axioms GenericCurrent.sym_contract_antisym
#print axioms GenericCurrent.sym_contract_antisym_apply
