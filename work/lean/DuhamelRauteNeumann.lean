import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecificLimits.Basic

/-!
Brückenlemmata zur charakteristisch abgeschlossenen Raute und zur Neumann-Reihe

Dieses Modul sichert zwei Bausteine der Zustandslösung des gestörten
Wellenproblems: Erstens ist das Kontrollgebiet $\mathcal D_{T,R}$, definiert
durch $|\alpha-\alpha_0|\le T$ und $|\phi|\le R+T-|\alpha-\alpha_0|$, unter
Charakteristikdreiecken abgeschlossen, es gibt also bei der Iteration keinen
Radienverlust. Zweitens gilt für die Iterierten die Potenzabschätzung
$x_k\le B^k x_0$, und der Reihenrest der geometrischen Majorante ist
$c^2/(1-c)$.

Das Modul formalisiert ausschließlich diese mathematischen Bausteine. Es sagt
nichts über das zugrunde liegende physikalische Modell, über die Glattheit des
Potentials oder über die Existenz der zugrunde gelegten Operatoren.
-/

namespace DuhamelRauteNeumann

/-- **Charakteristische Abgeschlossenheit des Rautegebiets.** Liegt
$(\alpha,\phi)$ in $\mathcal D_{T,R}$ und wählt man $s$ zwischen $\alpha_0$ und
$\alpha$ sowie $y$ mit $|y-\phi|\le|\alpha-s|$, so liegt $(s,y)$ wieder in
$\mathcal D_{T,R}$. Der Beweis verwendet nur die Dreiecksungleichung und die
Identität $|\alpha-s|=|\alpha-\alpha_0|-|s-\alpha_0|$ für $s$ zwischen
$\alpha_0$ und $\alpha$. -/
theorem raute_dreieck_inklusion
    {α₀ T R α φ s y : ℝ}
    (hα : |α - α₀| ≤ T)
    (hφ : |φ| ≤ R + T - |α - α₀|)
    (hsℓ : min α₀ α ≤ s) (hsu : s ≤ max α₀ α)
    (hy : |y - φ| ≤ |α - s|) :
    |s - α₀| ≤ T ∧ |y| ≤ R + T - |s - α₀| := by
  obtain ⟨hφℓ, hφr⟩ := abs_le.mp hφ
  obtain ⟨hyℓ, hyr⟩ := abs_le.mp hy
  rcases le_total α₀ α with h | h
  · -- Fall α₀ ≤ α: die Rautekante steigt, s liegt zwischen α₀ und α.
    have hmin : min α₀ α = α₀ := min_eq_left h
    have hmax : max α₀ α = α := max_eq_right h
    rw [hmin] at hsℓ
    rw [hmax] at hsu
    have hαeq : |α - α₀| = α - α₀ := abs_of_nonneg (by linarith)
    have hs0 : |s - α₀| = s - α₀ := abs_of_nonneg (by linarith)
    have has : |α - s| = α - s := abs_of_nonneg (by linarith)
    rw [hαeq] at hα hφℓ hφr
    refine ⟨by rw [hs0]; linarith, ?_⟩
    rw [hs0]
    exact abs_le.mpr ⟨by linarith, by linarith⟩
  · -- Fall α ≤ α₀: spiegelbildliche Rechnung mit umgekehrten Vorzeichen.
    have hmin : min α₀ α = α := min_eq_right h
    have hmax : max α₀ α = α₀ := max_eq_left h
    rw [hmin] at hsℓ
    rw [hmax] at hsu
    have hαeq : |α - α₀| = -(α - α₀) := abs_of_nonpos (by linarith)
    have hs0 : |s - α₀| = -(s - α₀) := abs_of_nonpos (by linarith)
    have has : |α - s| = -(α - s) := abs_of_nonpos (by linarith)
    rw [hαeq] at hα hφℓ hφr
    refine ⟨by rw [hs0]; linarith, ?_⟩
    rw [hs0]
    exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- Die Raute ist nach beiden Seiten in $\alpha$ abgeschlossen: Die
Vorzeichenwahl in `raute_dreieck_inklusion` ist symmetrisch, denn die Aussage
hängt nur an $|\alpha-\alpha_0|$ und an der Zwischenlage von $s$. -/
theorem raute_rand_schranken
    {α₀ T R α φ : ℝ}
    (hα : |α - α₀| ≤ T) (hφ : |φ| ≤ R + T - |α - α₀|) :
    |α - α₀| ≤ T ∧ α ≤ α₀ + T ∧ |φ| ≤ R + T := by
  refine ⟨hα, ?_, ?_⟩
  · linarith [(abs_le.mp hα).2]
  · linarith [abs_nonneg (α - α₀)]

/-- Die von der Rautenbedingung herkommende Randbreite ist nicht negativ; das
ist der Grund, warum $\mathcal D_{T,R}$ für zulässige Punkte nicht leer ist.
Die in der Konstruktion vorausgesetzten $T,R\ge0$ werden für diese Aussage nicht
gebraucht, sie folgen teilweise aus der Bedingung selbst. -/
theorem raute_randbreite_nichtnegativ
    {α₀ T R α φ : ℝ} (hφ : |φ| ≤ R + T - |α - α₀|) :
    0 ≤ R + T - |α - α₀| := by
  linarith [abs_nonneg φ]

/-- **Potenzabschätzung der Iteration.** Erfüllt eine nichtnegative Folge die
Rekursion $x_{k+1}\le B\,x_k$ mit $B\ge 0$, so gilt $x_k\le B^k x_0$. Angewandt
auf $x_k=\|\mathcal T^k u\|$ und $B=\|\mathcal T\|$ liefert das die Majorante
der Neumann-Reihe. -/
theorem folge_potenz_schranke {B : ℝ} (hB : 0 ≤ B) {x : ℕ → ℝ}
    (hrec : ∀ k, x (k + 1) ≤ B * x k) (k : ℕ) :
    x k ≤ B ^ k * x 0 := by
  induction k with
  | zero => simp
  | succ n ih =>
      calc x (n + 1) ≤ B * x n := hrec n
        _ ≤ B * (B ^ n * x 0) := mul_le_mul_of_nonneg_left ih hB
        _ = B ^ (n + 1) * x 0 := by ring

/-- **Rest der geometrischen Majorante.** Für $0\le c<1$ ist
$\sum_{k} c^{k+2}=c^{2}/(1-c)$. Zusammen mit der Potenzabschätzung ergibt das
die explizite Restschranke der Neumann-Reihe nach Abzug der Glieder
$k=0,1$. -/
theorem neumann_reihen_rest {c : ℝ} (h0 : 0 ≤ c) (hc : c < 1) :
    ∑' k : ℕ, c ^ (k + 2) = c ^ 2 / (1 - c) := by
  have hsum : HasSum (fun k : ℕ => c ^ k) (1 - c)⁻¹ := hasSum_geometric_of_lt_one h0 hc
  have hmul : HasSum (fun k : ℕ => c ^ 2 * c ^ k) (c ^ 2 * (1 - c)⁻¹) := hsum.mul_left (c ^ 2)
  have heq : (fun k : ℕ => c ^ (k + 2)) = (fun k : ℕ => c ^ 2 * c ^ k) :=
    funext fun k => by rw [pow_add, mul_comm]
  calc ∑' k : ℕ, c ^ (k + 2)
      = ∑' k : ℕ, c ^ 2 * c ^ k := by rw [heq]
      _ = c ^ 2 * (1 - c)⁻¹ := hmul.tsum_eq
      _ = c ^ 2 / (1 - c) := by field_simp

/-- Konvergenz der Majorante, wie sie für die Invertierbarkeit von
$I-\lambda\mathcal T$ in der Operatornorm gebraucht wird. -/
theorem neumann_majorante_summierbar {B lam : ℝ} (hB : 0 ≤ B) (hlam : 0 ≤ lam)
    (h : B * lam < 1) : Summable (fun k : ℕ => (B * lam) ^ k) :=
  (hasSum_geometric_of_lt_one (mul_nonneg hB hlam) h).summable

/- Das Axiomenaudit dieser Aussagen erfolgt außerhalb dieser deklarativen Datei
durch das Prüfskript, das eine temporäre Kopie mit angehängtem Audit baut. -/

end DuhamelRauteNeumann
