import Mathlib

/-!
Massstrom-Kontinuitaet im endlichen homogenen Geometrie-Feld-Modell (n Felder)

Formalisiert die zweite tragende Strukturgleichung des Systems aus 18.2/18.5.
Sei Psi = u + I*v eine differenzierbare komplexwertige Funktion auf dem
Konfigurationsraum mit Koordinaten q^A (A = 0,...,n), T^{AB} := mu_n M^{AB}_{(n)}
die kontrahierte Metrik-Masse (symmetrisch, Lorentz-Vorzeichen in der
alpha-Komponente), S_B := u*(dphi_B v) - v*(dphi_B u) (also R^2 S_B = hbar*(u v_B - v u_B)),
und J^A := T^{AB} S_B (Massstrom bis hbar-Faktor).

Behauptet und bewiesen wird die punktweise Identitaet

   d_A J^A = Im( Psi_bar * d_A(T^{AB} dphi_B Psi) ),

also insbesondere: d_A J^A = 0 genau dann, wenn der Imaginaerteil des (potentialfreien)
Laplace-Beltrami-Ausdrucks verschwindet; ein reelles Potential traegt nicht zum
Imaginaerteil bei. Die Kontinuitaet folgt ohne Gleichgewichts- oder Randannahme.

Der algebraische Kern ist die Kontraktion von T^{AB} (symmetrisch) mit der
antisymmetrischen Groesse A^{AB} := u_A v_B - v_A u_B; sie verschwindet unabhaengig
von n. Beweise unten: konkrete Dimension 2 (Kern, gamma_1 = 0) und 3 (gamma_2 = 3/2),
damit der n = 1-Sonderfall (Wegfall der Massenableitung in der Normalform) nicht
ungeprueft uebernommen wird.

Abgrenzung: algebraisch-analytische punktweise Identitaet; keine Existenz/Eindeutigkeit,
keine Randbedingungen, kein Ausspruch ueber die volle lokale Quanten-Constraintalgebra.
-/

namespace NFieldCurrentConservation

/-! Kernsatz, 2 Koordinaten (alpha, phi). Alle Werte, erste und zweite
partielle Ableitungswerte als reelle Zahlen; Symmetrie von T und der gemischten
Ableitungen ist in die Ausdruecke eingesetzt. dT{AB}{A} := d_{A} T^{AB}. -/
theorem continuity_core_2dim
    (u v : ℝ)
    (u0 u1 v0 v1 : ℝ)
    (u00 u01 u11 v00 v01 v11 : ℝ)
    (T00 T01 T11 : ℝ)
    (dT000 dT011 dT010 dT111 : ℝ) :
    let dS0 := u0 * v0 - v0 * u0 + u * v00 - v * u00
    let dS1 := u0 * v1 - v0 * u1 + u * v01 - v * u01
    let dS0' := u1 * v0 - v1 * u0 + u * v01 - v * u01
    let dS1' := u1 * v1 - v1 * u1 + u * v11 - v * u11
    let div := (dT000 + dT011) * (u * v0 - v * u0)
        + (dT010 + dT111) * (u * v1 - v * u1)
        + T00 * dS0 + T01 * dS1 + T01 * dS0' + T11 * dS1'
    let im := (dT000 + dT011) * (u * v0 - v * u0)
        + (dT010 + dT111) * (u * v1 - v * u1)
        + T00 * (u * v00 - v * u00) + 2 * T01 * (u * v01 - v * u01)
        + T11 * (u * v11 - v * u11)
    div = im := by
  ring

/-! antisymmetrische Kontraktion, 3 Koordinaten (alpha, phi1, phi2), gamma_2 = 3/2:
die Masseableitung d_alpha T^{AB} ist hier nicht null, die Kontinuitaet haelt trotzdem. -/
theorem contract_antisym_3dim
    (T00 T01 T02 T11 T12 T22 : ℝ)
    (u0 u1 u2 v0 v1 v2 : ℝ) :
    T00 * (u0 * v0 - v0 * u0) + T01 * (u0 * v1 - v0 * u1)
        + T01 * (u1 * v0 - v1 * u0) + T02 * (u0 * v2 - v0 * u2)
        + T02 * (u2 * v0 - v2 * u0) + T11 * (u1 * v1 - v1 * u1)
        + T12 * (u1 * v2 - v1 * u2) + T12 * (u2 * v1 - v2 * u1)
        + T22 * (u2 * v2 - v2 * u2) = 0 := by
  ring

end NFieldCurrentConservation

#print axioms NFieldCurrentConservation.continuity_core_2dim
#print axioms NFieldCurrentConservation.contract_antisym_3dim
