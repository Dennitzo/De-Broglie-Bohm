import Mathlib

/-!
Constraintfortpflanzung im endlichen homogenen Geometrie-Feld-Modell

Formalisierung des reduzierten homogenen Modells der Kapitel 14–18 der Publikation:
eine Geometriefreiheitsgrad `α` (mit `H = α̇`), `n` homogene Skalarfelder, gemeinsame
Laplace–Beltrami-Quantisierung, Trajektorienwirkung `L_Q = L_cl − N·Q_M`. Der
Hamilton-Zwang ist

```
E = H² − Σ_I v_I² − 2β·𝒰 − 2β·(e^{−3α}·Q_M)
```

Behauptet und bewiesen wird **nicht** `E = 0`, sondern die Identität

```
Ė = −3H·E ,     (e^{3α}·E)′ = 0
```

entlang jeder Lösung der Euler-Lagrange-Gleichungen (geometrische Gleichung und
`n` Feldgleichungen) für ein fest gehaltenes glattes `Q_M`. Es wird weder `E = 0`
vorausgesetzt noch der Ausbreitungsfaktor angenommen.

Der **algebraische Kern** (`nfield_constraint_core`) ist für **beliebiges `n`**
formuliert: die Felder erscheinen ausschließlich über die aggregierten Summen

```
Kv   = Σ_I v_I²,             dK    = Kv′  = −6H·Kv − 2β·dU − 2β·vsdQW
dU                          = Σ_I 𝒰_{,I}·v_I
vsdQW                       = e^{−3α}·Σ_I (Q_M)_{,I}·v_I
qaW                         = e^{−3α}·(Q_M)_{,α}
```

Die geometrische Gleichung (α-Variation) ist `2·Ḣ = −3H² − 3Kv + 6β𝒰 + 2β·qaW`
(doppelt genommen, um Brüche zu vermeiden; durch zwei teilbar). Der Kern ist eine
reine Polynomidentität (`ring`), gültig für jedes `n`.

Abgrenzung (was der Satz NICHT leistet)

* Er formalisiert die algebraische Fortpflanzung entlang einer hinreichend glatten
  Trajektorie; er beweist **keine Existenz oder Eindeutigkeit** von Lösungen und keine
  Regularität an `a → 0`.
* Die mehrdimensionale Kettenregel, die `𝒰(φ(t))` und `e^{−3α}Q_M(α(t),φ(t))` zu
  Funktionen der Zeit macht und deren Totalableitungen `dU`, `dQc` liefert, ist hier
  als **explizite Voraussetzung** (`hU`, `hQc` mit den angegebenen Totalableitungen)
  ausgewiesen und wird nicht aus einer `ℝ^{n+1}`-Differentiation von `Q_M` reproved.
  Das ist eine Standard-Tatsache der Differentialrechnung, kein Gegenstand dieses Satzes.
* Der Satz gilt für das **endliche homogene Minisuperspace-Modell**, nicht für die
  volle ADM- und Funktionstheorie; die lokale Quanten-Constraintalgebra der Operatoren bleibt
  offen. Ohne unbelegte Beweisschritte und ohne Zusatzaxome außer den üblichen.
-/

namespace NFieldConstraintPropagation

/-! Algebraischer Kern, gültig für jedes n -/

/-- **Algebraischer Kern.** Sind `H, Kv, U, Qc` die Werte und `dH, dK, dU, dQc` ihre
Totalableitungswerte entlang einer Trajektorie, und gelten die (doppelt genommene)
geometrische Gleichung `eH`, die aggregierte Feldgleichung `eK` und die Totalableitung
des gewichteten Quantenterms `eQ`, so ist die Totalableitung des Zwangs
`Ė = −3H·E` mit `E = H² − Kv − 2βU − 2βQc`. Rein algebraisch; `Kv` ist dabei
`Σ_I v_I²`, `vsdQW = e^{−3α}Σ_I(Q_M)_{,I}v_I` und `qaW = e^{−3α}(Q_M)_{,α}`. -/
theorem nfield_constraint_core
    (H Kv U Qc qaW vsdQW dH dK dU dQc β : ℝ)
    (eH : 2 * dH = -3 * H ^ 2 - 3 * Kv + 6 * β * U + 2 * β * qaW)
    (eK : dK = -6 * H * Kv - 2 * β * dU - 2 * β * vsdQW)
    (eQ : dQc = -3 * H * Qc + qaW * H + vsdQW) :
    2 * H * dH - dK - 2 * β * dU - 2 * β * dQc
      = -3 * H * (H ^ 2 - Kv - 2 * β * U - 2 * β * Qc) := by
  rw [show 2 * H * dH = H * (2 * dH) by ring, eH, eK, eQ]
  ring

/-- Der Zwang als Funktion der Zeit (Aggregate `Kv`, gewichteter Potentialwert `Uf`
und gewichteter Quantenterm `Qc`). -/
noncomputable def zwang (β : ℝ) (H Kv Uf Qc : ℝ → ℝ) (t : ℝ) : ℝ :=
  (H t) ^ 2 - Kv t - 2 * β * Uf t - 2 * β * Qc t

/-- **Analyse-Satz (Fortpflanzung).** Unter Differenzierbarkeit der Aggregate und den
Evolutionsgleichungen `eH, eK, eQ` (an `t`) ist die Ableitung von `zwang` gleich
`−3H·zwang`. Der Analysisteil liefert die Produkt- und Kettenregeln für `H²`, `Kv`, `Uf`,
`Qc`; der algebraische Inhalt ist exakt `nfield_constraint_core`. -/
theorem zwang_propagation
    (β : ℝ) (H Kv Uf Qc : ℝ → ℝ) (dH dK dU dQc qaW vsdQW : ℝ) (t : ℝ)
    (hH : HasDerivAt H dH t) (hKv : HasDerivAt Kv dK t)
    (hU : HasDerivAt Uf dU t) (hQc : HasDerivAt Qc dQc t)
    (eH : 2 * dH = -3 * (H t) ^ 2 - 3 * Kv t + 6 * β * Uf t + 2 * β * qaW)
    (eK : dK = -6 * (H t) * Kv t - 2 * β * dU - 2 * β * vsdQW)
    (eQ : dQc = -3 * (H t) * Qc t + qaW * (H t) + vsdQW) :
    HasDerivAt (zwang β H Kv Uf Qc)
      (-3 * H t * zwang β H Kv Uf Qc t) t := by
  have hbase := ((hH.pow 2).sub hKv).sub (hU.const_mul (2 * β)) |>.sub (hQc.const_mul (2 * β))
  refine HasDerivAt.congr_deriv hbase ?_
  have hcore := nfield_constraint_core (H t) (Kv t) (Uf t) (Qc t) qaW vsdQW dH dK dU dQc β eH eK eQ
  show 2 * H t ^ (2 - 1) * dH - dK - 2 * β * dU - 2 * β * dQc
      = -3 * H t * ((H t) ^ 2 - Kv t - 2 * β * Uf t - 2 * β * Qc t)
  rw [show H t ^ (2 - 1) = H t by norm_num]
  exact hcore

/-! Gewichtete Erhaltung und Null-Ausbreitung -/

/-- Der Hypothesensatz als punktweise Annahmen: Regularität (`HasDerivAt` für `α`, `H`, Aggregate)
und die drei Evolutionsgleichungen an jedem `t`. Die Aggregate `Uf, Qc` und ihre
Totalableitungen `dUf, dQc` tragen die (vorgabenmäßig ausgewiesene) Kettenregel. -/
structure Sys (β : ℝ) (α H Kv Uf Qc dHf dKf dUf dQcf qaWf vsdQWf : ℝ → ℝ) where
  hα : ∀ t, HasDerivAt α (H t) t
  hH : ∀ t, HasDerivAt H (dHf t) t
  hKv : ∀ t, HasDerivAt Kv (dKf t) t
  hUf : ∀ t, HasDerivAt Uf (dUf t) t
  hQc : ∀ t, HasDerivAt Qc (dQcf t) t
  eH : ∀ t, 2 * dHf t = -3 * (H t) ^ 2 - 3 * Kv t + 6 * β * Uf t + 2 * β * qaWf t
  eK : ∀ t, dKf t = -6 * (H t) * Kv t - 2 * β * dUf t - 2 * β * vsdQWf t
  eQ : ∀ t, dQcf t = -3 * (H t) * Qc t + qaWf t * (H t) + vsdQWf t

/-- Ableitung des Zwangs an jedem `t`: `Ė = −3H·E` (Gl. §15.4, verallgemeinert auf `n`). -/
theorem hasDerivAt_zwang (S : Sys β α H Kv Uf Qc dHf dKf dUf dQcf qaWf vsdQWf) (t : ℝ) :
    HasDerivAt (zwang β H Kv Uf Qc)
      (-3 * H t * zwang β H Kv Uf Qc t) t :=
  zwang_propagation β H Kv Uf Qc (dHf t) (dKf t) (dUf t) (dQcf t) (qaWf t) (vsdQWf t) t
    (S.hH t) (S.hKv t) (S.hUf t) (S.hQc t) (S.eH t) (S.eK t) (S.eQ t)

/-- Gewichtete Erhaltung: `(e^{3α}·E)′ = 0`. Produktregel, `3H·E + E·(−3H) = 0`; das
Gewicht ist gerade `e^{3α} = a³`, nicht eine andere Potenz. -/
theorem weighted_zwang_conserved
    (S : Sys β α H Kv Uf Qc dHf dKf dUf dQcf qaWf vsdQWf) (t : ℝ) :
    HasDerivAt (fun s => Real.exp (3 * α s) * zwang β H Kv Uf Qc s) 0 t := by
  have hexp : HasDerivAt (fun s => Real.exp (3 * α s))
      (Real.exp (3 * α t) * (3 * H t)) t := ((S.hα t).const_mul 3).exp
  convert hexp.mul (hasDerivAt_zwang S t) using 1
  ring

/-- `e^{3α}·E` ist auf ℝ konstant. -/
theorem weighted_zwang_const
    (S : Sys β α H Kv Uf Qc dHf dKf dUf dQcf qaWf vsdQWf) (s u : ℝ) :
    Real.exp (3 * α s) * zwang β H Kv Uf Qc s
      = Real.exp (3 * α u) * zwang β H Kv Uf Qc u := by
  have hf : Differentiable ℝ (fun s => Real.exp (3 * α s) * zwang β H Kv Uf Qc s) :=
    fun t => (weighted_zwang_conserved S t).differentiableAt
  have hd : ∀ t, deriv (fun s => Real.exp (3 * α s) * zwang β H Kv Uf Qc s) t = 0 :=
    fun t => (weighted_zwang_conserved S t).deriv
  exact is_const_of_deriv_eq_zero hf hd s u

/-- **Null-Ausbreitung (Konsistenz des geschlossenen Modells).** Gilt `E(t₀) = 0` an
einem Zeitpunkt, so `E(t) = 0` für alle `t`: die Constraintfläche ist invariant unter
der geführten Dynamik. Benutzt wird nur die Konstanz von `e^{3α}·E` und `exp ≠ 0`. -/
theorem zwang_eq_zero
    (S : Sys β α H Kv Uf Qc dHf dKf dUf dQcf qaWf vsdQWf)
    (t₀ t : ℝ) (h₀ : zwang β H Kv Uf Qc t₀ = 0) :
    zwang β H Kv Uf Qc t = 0 := by
  have h := weighted_zwang_const S t t₀
  rw [h₀, mul_zero] at h
  rcases eq_zero_or_eq_zero_of_mul_eq_zero h with hc | hc
  · exact absurd hc (Real.exp_ne_zero (3 * α t))
  · exact hc

end NFieldConstraintPropagation

#print axioms NFieldConstraintPropagation.nfield_constraint_core
#print axioms NFieldConstraintPropagation.zwang_propagation
#print axioms NFieldConstraintPropagation.hasDerivAt_zwang
#print axioms NFieldConstraintPropagation.weighted_zwang_conserved
#print axioms NFieldConstraintPropagation.weighted_zwang_const
#print axioms NFieldConstraintPropagation.zwang_eq_zero
