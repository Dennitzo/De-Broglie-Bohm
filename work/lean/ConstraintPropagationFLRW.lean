import Mathlib

/-!
# Ausbreitung des FLRW-Zwangs: `Ċ = −3HC`, `(a³C)' = 0` und Konstanz

Formalisierung der Abschnitte 29.6.2 und 29.6.2a der Publikation. Der Friedmann-Rest ist

```
C = H² − Λ/3 + k/a² − κρ/3 ,     ρ = ½Φ̇² + V(Φ) ,     H = ȧ/a            (29.6.2f)
```

Behauptet und bewiesen wird **nicht** `C = 0`, sondern die Identität

```
dC/dt = λ·C ,   λ = −3H = −K                                            (29.6.2g)
```

entlang jeder Lösung der beiden Evolutionsgleichungen (Raychaudhuri und Klein–Gordon).
Der Beweis differenziert die definierte Restgröße `C` tatsächlich: `Ḣ` kommt aus der
Raychaudhuri-Gleichung, `Φ̈` aus der Klein–Gordon-Gleichung, sonst nichts. Es wird weder
`C = 0` vorausgesetzt noch der Ausbreitungsfaktor angenommen.

Damit folgen der Erhaltungssatz `(a³C)' = 0` (29.6.2h), die Konstanz von `a³C` und als
Korollar: ist `C(t₀) = 0` für ein `t₀`, so `C ≡ 0` (29.6.2i).

## Abgrenzung

Der Satz betrifft das **endlich-dimensionale** FLRW-Minisuperspace mit den zwei
Evolutionsgleichungen der reduzierten Wirkung. Er ist kein Beweis der vollständigen
ADM-Zwangsalgebra (dazu braucht es die verschmierten Funktionale `H[N]`, `D[⃗N]` und ihre
Poisson-Klammern, Abschnitt 29.6.4), sichert keine Existenz von Lösungen, nichts über
Kollapsstellen `a → 0` und nichts über die Quantentheorie. Keine numerische Berechnung,
kein `sorry`.
-/

namespace ConstraintPropagationFLRW

/-! ### Die Größen von Gl. (29.6.2f) -/

/-- `H = ȧ/a`. -/
noncomputable def hubble (a da : ℝ → ℝ) (t : ℝ) : ℝ := da t / a t

/-- `ρ_Φ = ½Φ̇² + V(Φ)`. -/
noncomputable def energyDensity (Φ dΦ : ℝ → ℝ) (V : ℝ → ℝ) (t : ℝ) : ℝ :=
  (1 / 2 : ℝ) * (dΦ t) ^ 2 + V (Φ t)

/-- Die Zwangsgröße `C` von Gl. (29.6.2f). -/
noncomputable def constraint (k Λ κ : ℝ) (a da Φ dΦ : ℝ → ℝ) (V : ℝ → ℝ) (t : ℝ) : ℝ :=
  (hubble a da t) ^ 2 - Λ / 3 + k / (a t) ^ 2 - κ / 3 * energyDensity Φ dΦ V t

/-! ### Die Voraussetzungen: Regularität und die beiden Evolutionsgleichungen -/

/-- Regularität: `a` nirgends null, und `a`, `ȧ`, `Φ`, `Φ̇`, `V` sind differenzierbar mit den
angegebenen Ableitungsfunktionen. Das ist die präzise Fassung dessen, was im Manuskript als
„hinreichend glatt“ bezeichnet wird - hier steht sie ausdrücklich, damit der Satz keine
versteckte Annahme enthält. -/
structure Regular (a da dd Φ dΦ ddΦ V V' : ℝ → ℝ) where
  ha : ∀ t, a t ≠ 0
  had : ∀ t, HasDerivAt a (da t) t
  had' : ∀ t, HasDerivAt da (dd t) t
  hphi : ∀ t, HasDerivAt Φ (dΦ t) t
  hphi' : ∀ t, HasDerivAt dΦ (ddΦ t) t
  hV : ∀ x, HasDerivAt V (V' x) x

/-- Raychaudhuri-Gleichung in der Form von Gl. (29.6.2d):
`ä/a = −H²/2 − k/(2a²) + Λ/2 − (κ/4)(Φ̇² − 2V)`. -/
def Raychaudhuri (k Λ κ : ℝ) (a da dd Φ dΦ V : ℝ → ℝ) : Prop :=
  ∀ t, dd t / a t = -(da t / a t) ^ 2 / 2 - k / (2 * (a t) ^ 2) + Λ / 2
    - κ / 4 * ((dΦ t) ^ 2 - 2 * V (Φ t))

/-- Klein–Gordon-Gleichung von Gl. (29.6.2d): `Φ̈ = −3HΦ̇ − V'(Φ)`. -/
def KleinGordon (a da Φ dΦ ddΦ _V V' : ℝ → ℝ) : Prop :=
  ∀ t, ddΦ t = -3 * (da t / a t) * dΦ t - V' (Φ t)

/-! ### Der algebraische Kern von Gl. (29.6.2g) -/

/-- **Algebraischer Kern.** Haben `H, dH, ρ, dρ` die Werte und Richtungsableitungen der
Skalen mit `dH = −3H²/2 − k/(2a²) + Λ/2 − (κ/4)(Φ̇² − 2V)` (Raychaudhuri) und
`dρ = −3HΦ̇²` (Klein–Gordon), so gilt
`2H·dH − 2Hk/a² − κ·dρ/3 = −3H·(H² − Λ/3 + k/a² − κρ/3)`.
Das ist der gesamte Inhalt von (29.6.2g); der Analysisteil unten liefert nur die
Ableitungsregeln, die zu diesen Voraussetzungen führen. -/
theorem propagation_core (a H dH k Λ κ ρ dρ V dΦ : ℝ) (_ha : a ≠ 0)
    (hdH : dH = -(3 * H ^ 2) / 2 - k / (2 * a ^ 2) + Λ / 2 - κ / 4 * (dΦ ^ 2 - 2 * V))
    (hrho : ρ = (1 / 2 : ℝ) * dΦ ^ 2 + V) (hdρ : dρ = -3 * H * dΦ ^ 2) :
    2 * H * dH - 2 * H * k / a ^ 2 - κ * dρ / 3
      = -3 * H * (H ^ 2 - Λ / 3 + k / a ^ 2 - κ * ρ / 3) := by
  rw [hdρ, hrho, hdH]
  ring

/-! ### Der Analysisteil: tatsächliche Differentiation von `C` -/

section Sys

variable {k Λ κ : ℝ} {a da dd Φ dΦ ddΦ V V' : ℝ → ℝ}

/-- `H` ist differenzierbar mit `Ḣ = ä/a − H²`: Definition von `H` plus Quotientenregel,
noch ohne jede Feldgleichung. -/
theorem hasDerivAt_hubble (R : Regular a da dd Φ dΦ ddΦ V V') (t : ℝ) :
    HasDerivAt (hubble a da) (dd t / a t - (da t / a t) ^ 2) t := by
  unfold hubble
  convert (R.had' t).div (R.had t) (R.ha t) using 1
  field_simp [R.ha t]

/-- Die Energiedichte `ρ = ½Φ̇² + V(Φ)` mit Ableitung `Φ̈Φ̇ + V'(Φ)Φ̇`: reine Kettenregel,
ohne Feldgleichung. -/
theorem hasDerivAt_energyDensity (R : Regular a da dd Φ dΦ ddΦ V V') (t : ℝ) :
    HasDerivAt (energyDensity Φ dΦ V) (dΦ t * ddΦ t + V' (Φ t) * dΦ t) t := by
  unfold energyDensity
  convert ((((R.hphi' t).pow 2).const_mul (1 / 2 : ℝ)).add ((R.hV (Φ t)).comp t (R.hphi t)))
    using 1
  ring

/-- **Zwangsausbreitung**, Gl. (29.6.2g): die Ableitung von `C` entlang einer Lösung der
beiden Evolutionsgleichungen ist `-3H·C`. Bewiesen durch tatsächliche Differentiation von
`constraint`: Produkt- und Quotientenregel für `H²` und `k/a²`, Kettenregel für `ρ`, dann
je einmal Raychaudhuri und Klein–Gordon eingesetzt. -/
theorem hasDerivAt_constraint (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (t : ℝ) :
    HasDerivAt (constraint k Λ κ a da Φ dΦ V)
      (-3 * hubble a da t * constraint k Λ κ a da Φ dΦ V t) t := by
  have hH : HasDerivAt (hubble a da) (dd t / a t - (da t / a t) ^ 2) t :=
    hasDerivAt_hubble R t
  have hH2 : HasDerivAt (fun s => (hubble a da s) ^ 2)
      (2 * hubble a da t ^ 1 * (dd t / a t - (da t / a t) ^ 2)) t := hH.pow 2
  have hk : HasDerivAt (fun s => k / (a s) ^ 2) (-2 * (da t / a t) * k / (a t) ^ 2) t := by
    -- Quotientenregel auf die konstante Zählerfunktion und `(a s)²`; so bleibt die
    -- Rechnung in `/`-Schreibweise und es taucht kein punktweises Inverses auf.
    have h2 : HasDerivAt (fun s => (a s) ^ 2) (2 * (a t) ^ 1 * da t) t := (R.had t).pow 2
    convert (hasDerivAt_const t k).div h2 (pow_ne_zero 2 (R.ha t)) using 1
    field_simp [R.ha t, pow_ne_zero 2 (R.ha t)]
    ring
  have hrho : HasDerivAt (energyDensity Φ dΦ V) (dΦ t * ddΦ t + V' (Φ t) * dΦ t) t :=
    hasDerivAt_energyDensity R t
  -- Zusammensetzen der Beiträge von `H²`, `−Λ/3`, `k/a²` und `−κρ/3`.
  have hsum : HasDerivAt (constraint k Λ κ a da Φ dΦ V)
      (2 * hubble a da t ^ 1 * (dd t / a t - (da t / a t) ^ 2)
        - 2 * (da t / a t) * k / (a t) ^ 2
        - κ / 3 * (dΦ t * ddΦ t + V' (Φ t) * dΦ t)) t := by
    convert hH2.sub (hasDerivAt_const t (Λ / 3)) |>.add hk |>.sub (hrho.const_mul (κ / 3))
      using 1
    ring
  -- Jetzt die beiden Feldgleichungen einsetzen und den algebraischen Kern anrufen.
  have hdH : dd t / a t - (da t / a t) ^ 2 =
      -(3 * (da t / a t) ^ 2) / 2 - k / (2 * (a t) ^ 2) + Λ / 2
        - κ / 4 * ((dΦ t) ^ 2 - 2 * V (Φ t)) := by
    rw [hRay t]
    ring
  have hdρ : dΦ t * ddΦ t + V' (Φ t) * dΦ t = -3 * (da t / a t) * (dΦ t) ^ 2 := by
    rw [hKG t]
    ring
  convert hsum using 1
  rw [hdH, hdρ]
  unfold constraint energyDensity hubble
  field_simp [R.ha t]
  ring

/-- Gl. (29.6.2g) in `deriv`-Schreibweise. -/
theorem deriv_constraint (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (t : ℝ) :
    deriv (constraint k Λ κ a da Φ dΦ V) t
      = -3 * hubble a da t * constraint k Λ κ a da Φ dΦ V t :=
  (hasDerivAt_constraint R hRay hKG t).deriv

/-! ### Erhaltungssatz (29.6.2h) und Konstanz (29.6.2i) -/

/-- **Erhaltungssatz**, Gl. (29.6.2h): `(a³C)' = 0`. Die Produktregel liefert
`3a²ȧ·C + a³·(−3HC)`, und mit `H = ȧ/a` hebt sich das exakt - das ist der Grund, warum
gerade die Dichte `a³` und keine andere Potenz erscheint. -/
theorem hasDerivAt_a3_constraint (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (t : ℝ) :
    HasDerivAt (fun s => (a s) ^ 3 * constraint k Λ κ a da Φ dΦ V s) 0 t := by
  convert (((R.had t).pow 3).mul
    (hasDerivAt_constraint R hRay hKG t)) using 1
  simp only [Pi.pow_apply]
  unfold constraint energyDensity hubble
  field_simp [R.ha t]
  ring

/-- Erster Teil von Gl. (29.6.2i): `a³C` ist konstant. -/
theorem a3_constraint_const (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (s u : ℝ) :
    (a s) ^ 3 * constraint k Λ κ a da Φ dΦ V s
      = (a u) ^ 3 * constraint k Λ κ a da Φ dΦ V u := by
  have hf : Differentiable ℝ (fun s => (a s) ^ 3 * constraint k Λ κ a da Φ dΦ V s) :=
    fun t => (hasDerivAt_a3_constraint R hRay hKG t).differentiableAt
  have hd : ∀ t, deriv (fun s => (a s) ^ 3 * constraint k Λ κ a da Φ dΦ V s) t = 0 :=
    fun t => (hasDerivAt_a3_constraint R hRay hKG t).deriv
  exact is_const_of_deriv_eq_zero hf hd s u

/-- Zweiter Teil von Gl. (29.6.2i): gilt `C(t₀) = 0` für ein `t₀`, so ist `C(t) = 0` für
alle `t`. Benutzt wird allein die Konstanz von `a³C` und `a(t) ≠ 0`. -/
theorem constraint_eq_zero (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (t₀ t : ℝ) (h₀ : constraint k Λ κ a da Φ dΦ V t₀ = 0) :
    constraint k Λ κ a da Φ dΦ V t = 0 := by
  have h := a3_constraint_const R hRay hKG t t₀
  rw [h₀, mul_zero] at h
  rcases eq_zero_or_eq_zero_of_mul_eq_zero h with hc | hc
  · exfalso
    exact pow_ne_zero 3 (R.ha t) hc
  · exact hc

/-- Die explizite Zeitentwicklung des Zwangsrests, die zweite Form von Gl. (29.6.2i):
`C(t) = C(t₀)·(a(t₀)/a(t))³`. Der Rest klingt also mit `a⁻³` ab, ohne dass er je null
werden müsste. -/
theorem constraint_evolution (R : Regular a da dd Φ dΦ ddΦ V V')
    (hRay : Raychaudhuri k Λ κ a da dd Φ dΦ V) (hKG : KleinGordon a da Φ dΦ ddΦ V V')
    (t₀ t : ℝ) :
    constraint k Λ κ a da Φ dΦ V t
      = constraint k Λ κ a da Φ dΦ V t₀ * ((a t₀) / (a t)) ^ 3 := by
  have ht : (a t) ^ 3 ≠ 0 := pow_ne_zero 3 (R.ha t)
  have h := a3_constraint_const R hRay hKG t t₀
  have hr : (a t) ^ 3 * ((a t₀) / (a t)) ^ 3 = (a t₀) ^ 3 := by field_simp [R.ha t]
  apply mul_left_cancel₀ ht
  calc (a t) ^ 3 * constraint k Λ κ a da Φ dΦ V t
      = (a t₀) ^ 3 * constraint k Λ κ a da Φ dΦ V t₀ := h
    _ = ((a t) ^ 3 * ((a t₀) / (a t)) ^ 3) * constraint k Λ κ a da Φ dΦ V t₀ := by rw [hr]
    _ = (a t) ^ 3 * (constraint k Λ κ a da Φ dΦ V t₀ * ((a t₀) / (a t)) ^ 3) := by ring

end Sys

end ConstraintPropagationFLRW

#print axioms ConstraintPropagationFLRW.propagation_core
#print axioms ConstraintPropagationFLRW.hasDerivAt_constraint
#print axioms ConstraintPropagationFLRW.hasDerivAt_a3_constraint
#print axioms ConstraintPropagationFLRW.a3_constraint_const
#print axioms ConstraintPropagationFLRW.constraint_eq_zero
#print axioms ConstraintPropagationFLRW.constraint_evolution
