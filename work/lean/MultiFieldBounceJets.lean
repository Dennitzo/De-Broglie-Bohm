import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

/-!
Konstruktive Zweifeld-Instanz der Bounce-Kette (§18.12)

Dieses Modul sichert die algebraisch-differentialrechnungsmäßigen Bausteine
der konstruktiven Instanz für zwei gekoppelte Felder:

* Das Hauptsymbol des reellen Zweikomponentensystems der Zustandsgleichung
  (beide Komponenten tragen $M^{\alpha\alpha}_{(2)}=-\beta e^{-3\alpha}$) hat
  auf der Cauchy-Ebene $\alpha=x$ die Determinante $(\beta e^{-3x})^2\ne0$;
  die Cauchy-Ebene ist nichtcharakteristisch.
* Die analytischen Phasendaten $S(\alpha,\phi^1)=p\,\phi^1-\eta\,\phi^1(\alpha-x)$
  (Fortsetzung konstant in $\phi^2$) erfüllen die geforderten Jets:
  $\partial_\alpha S(x,0)=0$, $\partial_{\phi^1}S(x,0)=p$,
  $\partial_{\alpha}\partial_{\phi^1}S(x,0)=-\eta$,
  $\partial_{\alpha}\partial_{\phi^2}S=0$.
* Mit $V^{\alpha}=M^{\alpha\alpha}\partial_\alpha S=\beta e^{-3\alpha}\eta\phi^1$
  (hier `Vgeo`) und $V^{\phi^1}=M^{\phi^1\phi^1}\partial_{\phi^1}S$
  (hier `Vphi`) gilt die Punktidentität
  $\dot H(q_b)=\beta^2 e^{-6x}\,\eta\,p$, sie ist positiv unter $\eta p>0$
  und verschwindet für $\eta=0$ (Entartungsprobe).

Das Modul formalisiert ausschließlich diese Rechnungen. Es sagt nichts über
die Existenz der zugrunde gelegten analytischen Lösung (Cauchy-Kowalevski),
über die Positivität der Amplitude und nichts über das physikalische Modell.
-/

namespace MultiFieldBounceJets

/-- Hilfslemma: die Ableitung der linearen Funktion $x\mapsto c x$ ist konstant
$c$; die von `HasDerivAt.const_mul` erzeugte Form $c\cdot 1$ wird mit `ring`
aufgelöst. -/
theorem deriv_smul_id (c a : ℝ) : deriv (fun x => c * x) a = c := by
  convert ((hasDerivAt_id a).const_mul c).deriv using 1
  · ring

/-- **Nichtcharakteristik der Cauchy-Ebene.** Das Hauptsymbol des reellen
Zweikomponentensystems ist $M^{\alpha\alpha}_{(2)}(x)$ mal der Identität.
Seine Determinante ist in der Zweierformel $a\cdot a-0\cdot 0$ geschrieben,
also $(\beta e^{-3x})^2$, und für $\beta>0$ ungleich null. -/
theorem symbol_det_ne_zero (β x : ℝ) (hβ : 0 < β) :
    let a := -(β * Real.exp (-3 * x))
    a * a - (0 : ℝ) * 0 ≠ 0 := by
  have hc : -(β * Real.exp (-3 * x)) ≠ 0 := by
    have h : 0 < β * Real.exp (-3 * x) := by positivity
    linarith
  simp only [mul_zero, sub_zero]
  exact mul_ne_zero hc hc

/-- Die analytischen Cauchy-Phasendaten der Instanz: $S$ ist linear in
$\phi^1$ und in $\alpha$ affine mit Neigung $-\eta\phi^1$. -/
def phaseData (x η p : ℝ) (α φ₁ : ℝ) : ℝ := p * φ₁ - η * φ₁ * (α - x)

theorem phaseData_alpha_deriv (x η p α φ₁ : ℝ) :
    HasDerivAt (fun a => phaseData x η p a φ₁) (-(η * φ₁)) α := by
  have h₁ : HasDerivAt (fun a : ℝ => p * φ₁) 0 α := hasDerivAt_const α _
  have h₂ : HasDerivAt (fun a : ℝ => η * φ₁ * (a - x)) (η * φ₁) α := by
    convert ((hasDerivAt_id α).sub (hasDerivAt_const α x)).const_mul (η * φ₁) using 1
    · ring
  show HasDerivAt (fun a => p * φ₁ - η * φ₁ * (a - x)) (-(η * φ₁)) α
  convert h₁.sub h₂ using 1
  · ring

theorem phaseData_alpha_deriv_at (x η p α φ₁ : ℝ) :
    deriv (fun a => phaseData x η p a φ₁) α = -(η * φ₁) :=
  (phaseData_alpha_deriv x η p α φ₁).deriv

theorem phaseData_phi_deriv (x η p α f₀ : ℝ) :
    deriv (fun f => phaseData x η p α f) f₀ = p - η * (α - x) := by
  have : (fun f => phaseData x η p α f) = (fun f => (p - η * (α - x)) * f) := by
    funext f
    simp only [phaseData]
    ring
  rw [this]
  exact deriv_smul_id (p - η * (α - x)) f₀

/-- $\partial_\alpha S$ verschwindet auf der Symmetrieachse $\phi^1=0$,
insbesondere am Bahnpunkt $(x,0)$. -/
theorem phaseData_alpha_at_bounce (x η p : ℝ) :
    deriv (fun a => phaseData x η p a 0) x = 0 := by
  simpa using (phaseData_alpha_deriv_at x η p x 0)

/-- Gemischte Ableitung: $\partial_{\phi^1}\partial_\alpha S(x,0)=-\eta$. -/
theorem phaseData_alpha_phi_mixed (x η p : ℝ) :
    deriv (fun f => deriv (fun a => phaseData x η p a f) x) 0 = -η := by
  have h₁ : (fun f => deriv (fun a => phaseData x η p a f) x) = (fun f => -η * f) := by
    funext f
    rw [phaseData_alpha_deriv_at]
    simp
  rw [h₁]
  exact deriv_smul_id (-η) 0

/-- Uhrengradient am Bahnpunkt: $\partial_{\phi^1}S(x,0)=p$. -/
theorem phaseData_phi_at_bounce (x η p : ℝ) :
    deriv (fun f => phaseData x η p x f) 0 = p := by
  rw [phaseData_phi_deriv x η p x 0]
  simp

/-- Geführte geometrische Geschwindigkeit
$V^{\alpha}=-\beta e^{-3\alpha}\partial_\alpha S=\beta e^{-3\alpha}\eta\phi^1$. -/
noncomputable def Vgeo (β η : ℝ) (α φ₁ : ℝ) : ℝ := β * Real.exp (-3 * α) * η * φ₁

/-- Geführte Uhrfeldgeschwindigkeit
$V^{\phi^1}=\beta e^{-3\alpha}\partial_{\phi^1}S$. -/
noncomputable def Vphi (β x η p : ℝ) (α φ₁ : ℝ) : ℝ :=
  β * Real.exp (-3 * α) * (p - η * (α - x))

/-- Die Führung ist genau die Metrik mal den Daten-Jets:
$V^{\alpha}=M^{\alpha\alpha}\partial_\alpha S$. -/
theorem Vgeo_eq_symbol (β x η p α φ₁ : ℝ) :
    Vgeo β η α φ₁ = -(β * Real.exp (-3 * α)) * deriv (fun a => phaseData x η p a φ₁) α := by
  rw [phaseData_alpha_deriv_at]
  simp only [Vgeo]
  ring

/-- Und $V^{\phi^1}=M^{\phi^1\phi^1}\partial_{\phi^1}S$ am Bahnpunkt. -/
theorem Vphi_eq_symbol (β x η p α : ℝ) :
    Vphi β x η p α 0 = (β * Real.exp (-3 * α)) * deriv (fun f => phaseData x η p α f) 0 := by
  rw [phaseData_phi_deriv]
  simp only [Vphi]

theorem Vgeo_at_bounce (β η x : ℝ) : Vgeo β η x 0 = 0 := by simp [Vgeo]

theorem Vphi_at_bounce (β x η p : ℝ) :
    Vphi β x η p x 0 = β * Real.exp (-3 * x) * p := by
  simp [Vphi]

/-- **Punktidentität der Bounce-Ableitung.** Weil $V^{\alpha}(q_b)=0$ den
$\partial_\alpha V^{\alpha}$-Anteil tilgt, bleibt
$\dot H(q_b)=(\partial_{\phi^1}V^{\alpha})(q_b)\,V^{\phi^1}(q_b)
=\beta^2 e^{-6x}\eta p$. -/
theorem hdot_identity (β x η p : ℝ) :
    (deriv (fun a => Vgeo β η a 0) x) * Vgeo β η x 0
      + (deriv (fun f => Vgeo β η x f) 0) * Vphi β x η p x 0
      = β ^ 2 * Real.exp (-6 * x) * η * p := by
  have h₁ : deriv (fun a => Vgeo β η a 0) x = 0 := by
    have : (fun a => Vgeo β η a 0) = fun _ => (0 : ℝ) := by
      funext a
      simp only [Vgeo, mul_zero]
    rw [this]
    simp
  have h₂ : deriv (fun f => Vgeo β η x f) 0 = β * Real.exp (-3 * x) * η := by
    have : (fun f => Vgeo β η x f) = (fun f => (β * Real.exp (-3 * x) * η) * f) := by
      funext f
      simp only [Vgeo]
    rw [this]
    exact deriv_smul_id (β * Real.exp (-3 * x) * η) 0
  rw [h₁, h₂, Vphi_at_bounce]
  simp only [zero_mul, zero_add]
  have e2 : Real.exp (-3 * x) * Real.exp (-3 * x) = Real.exp (-6 * x) := by
    rw [← Real.exp_add]; ring
  have h3 : (β * Real.exp (-3 * x) * η) * (β * Real.exp (-3 * x) * p)
      = β ^ 2 * (Real.exp (-3 * x) * Real.exp (-3 * x)) * η * p := by ring
  rw [h3, e2]

/-- **Nichtentartung unter $\eta p>0$.** -/
theorem hdot_pos (β x η p : ℝ) (hβ : 0 < β) (hηp : 0 < η * p) :
    0 < β ^ 2 * Real.exp (-6 * x) * η * p := by
  have he : 0 < Real.exp (-6 * x) := Real.exp_pos _
  have hb : 0 < β ^ 2 := pow_pos hβ 2
  have hprod : 0 < (β ^ 2 * Real.exp (-6 * x)) * (η * p) :=
    mul_pos (mul_pos hb he) hηp
  have hrew : β ^ 2 * Real.exp (-6 * x) * η * p = (β ^ 2 * Real.exp (-6 * x)) * (η * p) :=
    by ring
  rw [hrew]
  exact hprod

/-- **Entartungsprobe.** Für $\eta=0$ verschwindet $\dot H(q_b)$; die
Nichtentartung hängt genau an der gemischten Datenableitung. -/
theorem hdot_vanishes_when_eta_zero (β x p : ℝ) :
    (deriv (fun a => Vgeo β 0 a 0) x) * Vgeo β 0 x 0
      + (deriv (fun f => Vgeo β 0 x f) 0) * Vphi β x 0 p x 0 = 0 := by
  simp [Vgeo, Vphi]

end MultiFieldBounceJets
