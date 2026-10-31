import Mathlib

/-!
Kern der Polardarstellung (Phase 7).

Formalisiert werden vier Umrechnungsidentitäten zwischen den kartesischen
Komponenten `u, v` eines Zustands und seiner Polardarstellung `R, theta`:

  * die Ableitung des Betrags `sqrt (u^2 + v^2)`,
  * die Polarkernidentitaet `u * v' - v * u' = R^2 * theta'`,
  * die Stromidentitaet `mu * R^2 * (hbar * theta)' = mu * hbar * (u * v' - v * u')`,
  * und die Positivaet des Polkernnenners unter einer gleichmaessigen
    Amplitude-Untergrenze.

Alle vier Saetze sind identitaetsgleiche Aussagen der Differentialrechnung
entlang einer Parameterlinie. Sie sagen nichts ueber die Existenz, die
Eindeutigkeit oder die Regularitaet einer zugrunde liegenden Loesung aus;
die Existenz bleibt analytischer Inhalt (Cauchy-Kowalevski).
-/

namespace PolarKernel

/-- Ableitung des Betrags aus den kartesischen Komponenten: Ist
`u t0 ^ 2 + v t0 ^ 2 != 0`, so ist `sqrt (u^2 + v^2)` an `t0` differenzierbar
mit Ableitung `(u u' + v v') / sqrt (u^2 + v^2)`. -/
theorem hasDerivAt_sqrt_sum_squares {u v : ℝ → ℝ} {t₀ r₀ v₀ : ℝ}
    (hu : HasDerivAt u r₀ t₀) (hv : HasDerivAt v v₀ t₀)
    (hne : u t₀ ^ 2 + v t₀ ^ 2 ≠ 0) :
    HasDerivAt (fun t => Real.sqrt (u t ^ 2 + v t ^ 2))
      ((u t₀ * r₀ + v t₀ * v₀) / Real.sqrt (u t₀ ^ 2 + v t₀ ^ 2)) t₀ := by
  have h1 : HasDerivAt (fun t => u t ^ 2) (2 * u t₀ * r₀) t₀ := by
    convert hu.pow 2 using 1; ring
  have h2 : HasDerivAt (fun t => v t ^ 2) (2 * v t₀ * v₀) t₀ := by
    convert hv.pow 2 using 1; ring
  have hsum : HasDerivAt (fun t => u t ^ 2 + v t ^ 2) (2 * u t₀ * r₀ + 2 * v t₀ * v₀) t₀ :=
    h1.add h2
  have hpos : 0 < u t₀ ^ 2 + v t₀ ^ 2 := lt_of_le_of_ne (by positivity) (Ne.symm hne)
  have hsqrt_ne : Real.sqrt (u t₀ ^ 2 + v t₀ ^ 2) ≠ 0 := (Real.sqrt_pos.mpr hpos).ne'
  convert (Real.hasDerivAt_sqrt hne).comp t₀ hsum using 1
  field_simp [hne, hsqrt_ne]

/-- Polarkernidentitaet: Fuer `u = R cos theta` und `v = R sin theta` mit
Ableitungen `r0` bzw. `theta0` an `t0` gilt `u * v' - v * u' = R^2 * theta'`.
Damit ist die Phasenableitung gerade die antisymmetrische Kombination der
Komponenten, und die Division durch `R^2` ist die einzig moegliche Umrechnung. -/
theorem polar_wronskian {R θ : ℝ → ℝ} {t₀ r₀ θ₀ : ℝ}
    (hR : HasDerivAt R r₀ t₀) (hθ : HasDerivAt θ θ₀ t₀) :
    let u := fun t => R t * Real.cos (θ t)
    let v := fun t => R t * Real.sin (θ t)
    u t₀ * deriv v t₀ - v t₀ * deriv u t₀ = R t₀ ^ 2 * θ₀ := by
  intro u v
  have hu : HasDerivAt u (r₀ * Real.cos (θ t₀) - R t₀ * Real.sin (θ t₀) * θ₀) t₀ := by
    convert hR.mul (hθ.cos) using 1
    ring
  have hv : HasDerivAt v (r₀ * Real.sin (θ t₀) + R t₀ * Real.cos (θ t₀) * θ₀) t₀ := by
    convert hR.mul (hθ.sin) using 1
    ring
  rw [hu.deriv, hv.deriv]
  linear_combination (Real.cos_sq_add_sin_sq (θ t₀)) * (R t₀ * R t₀ * θ₀)

/-- Stromidentitaet: Mit der Phase `S = hbar * theta` ist `mu * R^2 * S'` gleich
`mu * hbar * (u v' - v u')`. Das ist die im Manuskript benutzte Umrechnung des
Massstroms von der Polar- in die Kartesischform; der Vorfaktor `hbar` sitzt genau
einmal. -/
theorem polar_measure_current (μ ħ : ℝ) {R θ : ℝ → ℝ} {t₀ r₀ θ₀ : ℝ}
    (hR : HasDerivAt R r₀ t₀) (hθ : HasDerivAt θ θ₀ t₀) :
    let u := fun t => R t * Real.cos (θ t)
    let v := fun t => R t * Real.sin (θ t)
    μ * R t₀ ^ 2 * (ħ * θ₀) = μ * ħ * (u t₀ * deriv v t₀ - v t₀ * deriv u t₀) := by
  intro u v
  have hkey := polar_wronskian (R := R) (θ := θ) (t₀ := t₀) (r₀ := r₀) (θ₀ := θ₀) hR hθ
  rw [hkey]
  ring

/-- Positivaet des Polkernnenners: Wo die Amplitude nach unten durch `R0`
abgeschaetzt ist, ist auch `u^2 + v^2` gleichmaessig positiv, und die
Umrechnungsformel der Phasenableitung ist damit ueberall auf diesem Gebiet
zulaessig. -/
theorem amplit_bound_transfer (R₀ : ℝ) {u v : ℝ → ℝ} {t₀ : ℝ}
    (hR₀ : 0 < R₀) (hbound : R₀ ^ 2 ≤ u t₀ ^ 2 + v t₀ ^ 2) :
    0 < u t₀ ^ 2 + v t₀ ^ 2 := by
  nlinarith

end PolarKernel
