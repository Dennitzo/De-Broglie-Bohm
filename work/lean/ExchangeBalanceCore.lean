import Mathlib

/-!
Verallgemeinerte Austauschbilanz im endlichen homogenen Geometrie-Feld-Modell

Formalisierung des algebraischen Kerns der verallgemeinerten Austauschbilanz
(Kapitel 20 der Publikation): für `n` homogene Skalarfelder mit gemeinsamem glattem
Potential und festem glattem Quantenpotential `Q_M` ist die Divergenz der effektiven
Quelle der geführten Konfiguration identisch null,

```
d/dt(ρ_F + ρ_Q) + 3H (ρ_F + p_F + ρ_Q + p_Q) = 0 ,
```

wobei die Feldgleichungen `v̇_I = −3H v_I − β·𝒰_{,I} − β·e^{−3α}·Q_{M,I}` und die
geometrische Gleichung aus der Trajektorienwirkung `L_Q = L_cl − N·Q_M` gelten. Das
verallgemeinert die Einfeld-Austauschbilanz des §15.3 auf beliebige Feldzahl.

Die Felder erscheinen über die Aggregate

```
Kv      = Σ_I v_I²,              dK  = (Σ_I v_I²)′ = 2 Σ_I v_I v̇_I
dU      = Σ_I 𝒰_{,I} v_I
vsdQW   = e^{−3α} Σ_I (Q_M)_{,I} v_I
qaW     = e^{−3α} (Q_M)_{,α}
Qc      = e^{−3α} Q_M
```

Die aggregierte Feldgleichung `dK = −6H·Kv − 2β·dU − 2β·vsdQW` und die Totalableitung
des gewichteten Quantenterms `dQc = −3H·Qc + qaW·H + vsdQW` sind die beiden
Bausteine. Der Kern ist eine reine Polynomidentität (`ring`), gültig für jedes `n`.

Abgrenzung (was der Satz NICHT leistet)

* Der Satz formalisiert die algebraische Bilanz entlang einer hinreichend glatten
  Trajektorie; er beweist **keine Existenz oder Eindeutigkeit** von Lösungen und keine
  Regularität an `a → 0`.
* Die mehrdimensionale Kettenregel, die `𝒰(φ(t))` und `e^{−3α}Q_M(α(t),φ(t))` zu
  Funktionen der Zeit macht und deren Totalableitungen `dU`, `dQc` liefert, ist als
  explizite Voraussetzung (`eK`, `eQ` mit den angegebenen Totalableitungen) ausgewiesen
  und wird nicht aus einer `ℝ^{n+1}`-Differentiation von `Q_M` reproved.
* Der Satz gilt für das **endliche homogene Minisuperspace-Modell**, nicht für die
  volle ADM- und Funktionstheorie. Ohne unbelegte Beweisschritte und ohne Zusatzaxiome
  außer den üblichen.
-/

namespace ExchangeBalanceCore

/-- **Feldbilanz (algebraischer Kern).** Sind `H, Kv` die Aggregate und `dK` die
Totalableitung von `Kv`, und gilt die aggregierte Feldgleichung `eK`, so ist die
Divergenz des klassischen Feldanteils gleich dem negativen Quantenaustausch
`−β·vsdQW`. Rein algebraisch; die Feldzahl `n` geht nur über die Aggregate ein. -/
theorem field_balance_core
    (H Kv vsdQW dK dU β : ℝ)
    (eK : dK = -6 * H * Kv - 2 * β * dU - 2 * β * vsdQW) :
    (1/2) * dK + β * dU + 3 * H * Kv = -β * vsdQW := by
  rw [eK]
  ring

/-- **Quantenbilanz (algebraischer Kern).** Mit `Qc` und der Totalableitung `dQc`
des gewichteten Quantenterms sowie `qaW = e^{−3α}(Q_M)_{,α}` ist die Divergenz des
Quantenanteils (Druck `p_Q = −(e^{−3α}/3)·Q_{M,α}`) gleich dem positiven
Quantenaustausch `β·vsdQW`. Rein algebraisch. -/
theorem quantum_balance_core
    (H Qc qaW vsdQW dQc β : ℝ)
    (eQ : dQc = -3 * H * Qc + qaW * H + vsdQW) :
    β * dQc + 3 * β * H * Qc - β * H * qaW = β * vsdQW := by
  rw [eQ]
  ring

/-- **Totale Divergenz (Austauschbilanz, verallgemeinert).** Sind beide Bausteine
erfüllt — aggregierte Feldgleichung `eK` und Totalableitung des gewichteten
Quantenterms `eQ` — so heben sich die beiden Austauschbeiträge exakt weg, und die
Divergenz der effektiven Quelle der geführten Konfiguration ist null. Gültig für
jedes `n`, da nur die Aggregate eingehen. -/
theorem exchange_core
    (H Kv Qc qaW vsdQW dK dU dQc β : ℝ)
    (eK : dK = -6 * H * Kv - 2 * β * dU - 2 * β * vsdQW)
    (eQ : dQc = -3 * H * Qc + qaW * H + vsdQW) :
    (1/2) * dK + β * dU + 3 * H * Kv + β * dQc + 3 * β * H * Qc - β * H * qaW = 0 := by
  rw [eK, eQ]
  ring

end ExchangeBalanceCore

#print axioms ExchangeBalanceCore.field_balance_core
#print axioms ExchangeBalanceCore.quantum_balance_core
#print axioms ExchangeBalanceCore.exchange_core
