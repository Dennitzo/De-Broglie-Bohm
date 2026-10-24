import Mathlib

/-! Passive boost x'=g(x-vt), t'=g(t-vx), metric signature (+,-).
The two-index current transforms contravariantly; the normal covectors
transform by the inverse boost. All identities retain those normals. -/

namespace RelativisticLorentz

def boostTime (g v t x : ℝ) : ℝ := g * (t - v*x)
def boostSpace (g v t x : ℝ) : ℝ := g * (x - v*t)
def normalTime (g v n0 n1 : ℝ) : ℝ := g * (n0 + v*n1)
def normalSpace (g v n0 n1 : ℝ) : ℝ := g * (n1 + v*n0)
def tensor00 (g v A B C D : ℝ) : ℝ := g^2*(A-v*B-v*C+v^2*D)
def tensor10 (g v A B C D : ℝ) : ℝ := g^2*(B-v*A-v*D+v^2*C)
def tensor01 (g v A B C D : ℝ) : ℝ := g^2*(C-v*A-v*D+v^2*B)
def tensor11 (g v A B C D : ℝ) : ℝ := g^2*(D-v*B-v*C+v^2*A)
def contraction (A B C D n0 n1 m0 m1 : ℝ) : ℝ :=
  A*n0*m0+B*n1*m0+C*n0*m1+D*n1*m1

def boostWeight (s t : ℝ) (z : ℂ) : ℝ :=
  Complex.normSq ((s : ℂ)*(t : ℂ)*z)

theorem spinor_pair_weight (s t : ℝ) (z : ℂ) :
    boostWeight s t z = s^2*t^2*Complex.normSq z := by
  simp only [boostWeight, Complex.normSq_mul, Complex.normSq_ofReal]
  ring

def spinor00 (a b c d : ℂ) : ℝ :=
  Complex.normSq a + Complex.normSq b + Complex.normSq c + Complex.normSq d
def spinor10 (a b c d : ℂ) : ℝ :=
  Complex.normSq a + Complex.normSq b - Complex.normSq c - Complex.normSq d
def spinor01 (a b c d : ℂ) : ℝ :=
  Complex.normSq a - Complex.normSq b + Complex.normSq c - Complex.normSq d
def spinor11 (a b c d : ℂ) : ℝ :=
  Complex.normSq a - Complex.normSq b - Complex.normSq c + Complex.normSq d

/- S = diag(s,t) acts on each particle. This result derives all four
   tensor-current entries from the transformed complex amplitudes. -/
theorem spinor_boost_tensor_covariance (g v s t : ℝ) (a b c d : ℂ)
    (hs : s^2=g*(1-v)) (ht : t^2=g*(1+v)) :
    boostWeight s s a + boostWeight s t b + boostWeight t s c + boostWeight t t d =
      tensor00 g v (spinor00 a b c d) (spinor10 a b c d) (spinor01 a b c d) (spinor11 a b c d) ∧
    boostWeight s s a + boostWeight s t b - boostWeight t s c - boostWeight t t d =
      tensor10 g v (spinor00 a b c d) (spinor10 a b c d) (spinor01 a b c d) (spinor11 a b c d) ∧
    boostWeight s s a - boostWeight s t b + boostWeight t s c - boostWeight t t d =
      tensor01 g v (spinor00 a b c d) (spinor10 a b c d) (spinor01 a b c d) (spinor11 a b c d) ∧
    boostWeight s s a - boostWeight s t b - boostWeight t s c + boostWeight t t d =
      tensor11 g v (spinor00 a b c d) (spinor10 a b c d) (spinor01 a b c d) (spinor11 a b c d) := by
  simp only [spinor_pair_weight, hs, ht]
  unfold tensor00 tensor10 tensor01 tensor11 spinor00 spinor10 spinor01 spinor11
  constructor
  · ring
  constructor
  · ring
  constructor <;> ring

theorem lorentz_metric_preserved (g v t x : ℝ) (hg : g^2*(1-v^2)=1) :
    (boostTime g v t x)^2 - (boostSpace g v t x)^2 = t^2-x^2 := by
  calc
    _ = (g^2*(1-v^2))*(t^2-x^2) := by unfold boostTime boostSpace; ring
    _ = _ := by rw [hg]; ring

theorem boost_inverse (g v t x : ℝ) (hg : g^2*(1-v^2)=1) :
    boostTime g (-v) (boostTime g v t x) (boostSpace g v t x) = t ∧
    boostSpace g (-v) (boostTime g v t x) (boostSpace g v t x) = x := by
  constructor
  · calc
      _ = (g^2*(1-v^2))*t := by unfold boostTime boostSpace; ring
      _ = _ := by rw [hg]; ring
  · calc
      _ = (g^2*(1-v^2))*x := by unfold boostTime boostSpace; ring
      _ = _ := by rw [hg]; ring

theorem vector_normal_contraction_invariant (g v t x n0 n1 : ℝ)
    (hg : g^2*(1-v^2)=1) :
    boostTime g v t x * normalTime g v n0 n1 +
      boostSpace g v t x * normalSpace g v n0 n1 = t*n0+x*n1 := by
  calc
    _ = (g^2*(1-v^2))*(t*n0+x*n1) := by
      unfold boostTime boostSpace normalTime normalSpace
      ring
    _ = _ := by rw [hg]; ring

theorem tensor_density_invariant (g v A B C D n0 n1 m0 m1 : ℝ)
    (hg : g^2*(1-v^2)=1) :
    contraction (tensor00 g v A B C D) (tensor10 g v A B C D)
      (tensor01 g v A B C D) (tensor11 g v A B C D)
      (normalTime g v n0 n1) (normalSpace g v n0 n1)
      (normalTime g v m0 m1) (normalSpace g v m0 m1) =
    contraction A B C D n0 n1 m0 m1 := by
  calc
    _ = (g^2*(1-v^2))^2 * contraction A B C D n0 n1 m0 m1 := by
      unfold contraction tensor00 tensor10 tensor01 tensor11 normalTime normalSpace
      ring
    _ = _ := by rw [hg]; ring

theorem transformed_rest_density (g v A B C D : ℝ)
    (hg : g^2*(1-v^2)=1) :
    contraction (tensor00 g v A B C D) (tensor10 g v A B C D)
      (tensor01 g v A B C D) (tensor11 g v A B C D) g (g*v) g (g*v) = A := by
  have h := tensor_density_invariant g v A B C D 1 0 1 0 hg
  simpa [normalTime, normalSpace, contraction] using h

theorem contracted_current_time_covariant (g v A B C D m0 m1 : ℝ)
    (hg : g^2*(1-v^2)=1) :
    tensor00 g v A B C D * normalTime g v m0 m1 +
      tensor01 g v A B C D * normalSpace g v m0 m1 =
    boostTime g v (A*m0+C*m1) (B*m0+D*m1) := by
  calc
    _ = (g^2*(1-v^2))*boostTime g v (A*m0+C*m1) (B*m0+D*m1) := by
      unfold tensor00 tensor01 normalTime normalSpace boostTime
      ring
    _ = _ := by rw [hg]; ring

theorem contracted_current_space_covariant (g v A B C D m0 m1 : ℝ)
    (hg : g^2*(1-v^2)=1) :
    tensor10 g v A B C D * normalTime g v m0 m1 +
      tensor11 g v A B C D * normalSpace g v m0 m1 =
    boostSpace g v (A*m0+C*m1) (B*m0+D*m1) := by
  calc
    _ = (g^2*(1-v^2))*boostSpace g v (A*m0+C*m1) (B*m0+D*m1) := by
      unfold tensor10 tensor11 normalTime normalSpace boostSpace
      ring
    _ = _ := by rw [hg]; ring

theorem velocity_denominator_positive (u v : ℝ) (hu : |u| ≤ 1) (hv : |v| < 1) :
    0 < 1-v*u := by
  have hprod : |v*u| < 1 := calc
    |v*u| = |v| * |u| := abs_mul v u
    _ ≤ |v| * 1 := mul_le_mul_of_nonneg_left hu (abs_nonneg v)
    _ < 1 := by simpa using hv
  have h := le_abs_self (v*u)
  linarith

theorem velocity_addition_causal (u v : ℝ) (hu : |u| ≤ 1) (hv : |v| < 1) :
    |(u-v)/(1-v*u)| ≤ 1 := by
  have hden := velocity_denominator_positive u v hu hv
  obtain ⟨hul, huh⟩ := abs_le.mp hu
  obtain ⟨hvl, hvh⟩ := abs_lt.mp hv
  have hp := mul_nonneg (by linarith : 0 ≤ 1-u) (by linarith : 0 ≤ 1+v)
  have hm := mul_nonneg (by linarith : 0 ≤ 1+u) (by linarith : 0 ≤ 1-v)
  have hnum : |u-v| ≤ 1-v*u := by
    apply abs_le.mpr
    constructor <;> nlinarith
  rw [abs_div, abs_of_pos hden]
  exact (div_le_one hden).mpr hnum

theorem velocity_addition_identity (g v t x : ℝ) (hg : g ≠ 0) (ht : t ≠ 0) :
    boostSpace g v t x / boostTime g v t x = (x/t-v)/(1-v*(x/t)) := by
  unfold boostSpace boostTime
  field_simp

theorem rest_normal_future_unit (g v : ℝ) (hg : 0 < g) (hv : |v| < 1)
    (hunit : g^2*(1-v^2)=1) : g^2-(g*v)^2=1 ∧ |g*v|<g := by
  constructor
  · nlinarith [hunit]
  · rw [abs_mul, abs_of_pos hg]
    nlinarith

end RelativisticLorentz
