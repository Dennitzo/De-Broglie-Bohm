import Mathlib

/-! The inverse-metric chart for the subcritical specified source.
Endpoint existence and its geometric interpretation remain analytic inputs. -/
namespace FourDimEndpoint
noncomputable section
open Filter
open scoped Topology

def D (L delta r u : ℝ) : ℝ := r^2*u^2*(1-u)+L*(1-u)^2*
  (delta+(2*delta+1/15)*u+(5*delta+7/15)*u^2)
def P (delta u : ℝ) : ℝ := delta+6*delta*u+(9*delta-2/5)*u^2
def Pu (delta u : ℝ) : ℝ := 6*delta+2*(9*delta-2/5)*u

def N (L delta r u : ℝ) : ℝ := r^2*u^2*(1-u)-L*(1-u)^2*P delta u/3
def Nr (r u : ℝ) : ℝ := 2*r*u^2*(1-u)
def Nu (L delta r u : ℝ) : ℝ := r^2*(2*u-3*u^2)+
  2*L*(1-u)*P delta u/3-L*(1-u)^2*Pu delta u/3

def H (L delta r u : ℝ) : ℝ := D L delta r u/(r^3*u^2)
def G (L delta r u : ℝ) : ℝ := N L delta r u/(r^3*u^3)
def Gdot (L delta r u : ℝ) : ℝ :=
  (r*Nr r u-3*N L delta r u)/(r^4*u^3)+
  (u*Nu L delta r u-3*N L delta r u)*D L delta r u/(r^6*u^6)
def lapseRate (L delta r u : ℝ) : ℝ :=
  (4*r^2*u*(1-u)+L*(1-u)^2*((1/15-4*delta)+(13/15-4*delta)*u))/(3*D L delta r u)
def radialTime (L delta r u : ℝ) : ℝ :=
  u*Gdot L delta r u/2+u*(G L delta r u)^2/4+H L delta r u*G L delta r u/4

def normA (L delta r u : ℝ) : ℝ :=
  (2*r^3*u^3*Nr r u-6*r^2*u^3*N L delta r u+
    2*u*D L delta r u*Nu L delta r u-5*D L delta r u*N L delta r u+
    (N L delta r u)^2)/(4*r^6)
def curvature (L delta r u : ℝ) : ℝ :=
  4*((radialTime L delta r u)^2+2*(u*G L delta r u/(2*r))^2+
    2*(H L delta r u/(2*r))^2+((1-u)/r^2)^2)
def normK (L delta r u : ℝ) : ℝ :=
  4*((normA L delta r u)^2+(N L delta r u)^2*u^6/(2*r^8)+
    (D L delta r u)^2*u^6/(2*r^8)+(1-u)^2*u^10/r^4)

theorem source_chart_rate_identities (L delta r u : ℝ) (hr : r ≠ 0) (hu : u ≠ 0) :
    -(1/u)*(1/u-1)/r-L*(1/u-1)^2*(delta*((1/u)^2+2/u+5)+(1/u+7)/15)/r^3 =
      -D L delta r u/(r^3*u^4) ∧
    (1/u-1)/r+L*(1/u-1)^2*(-delta*(1/u+3)^2+2/5)/(3*(1/u)*r^3)=G L delta r u := by
  unfold D G N P
  constructor <;> field_simp <;> ring

theorem denominator_positive (L delta r u : ℝ)
    (hL : 0 < L) (hd : 0 < delta) (hu : 0 < u) (hu1 : u < 1) :
    0 < D L delta r u := by
  unfold D
  have hv : 0 < 1-u := sub_pos.mpr hu1
  positivity

theorem denominator_lower_bound (L delta r u : ℝ)
    (hL : 0 ≤ L) (hd : 0 ≤ delta) (hu : 0 ≤ u) (hu1 : u ≤ 1) :
    L*delta*(1-u)^2 ≤ D L delta r u := by
  have ha : 0 ≤ r^2*u^2*(1-u) := by positivity
  have hb : 0 ≤ L*(1-u)^2*((2*delta+1/15)*u+(5*delta+7/15)*u^2) := by positivity
  unfold D
  nlinarith

theorem N_hasDerivAt (r u : ℝ → ℝ) (L delta x dr du : ℝ)
    (hr : HasDerivAt r dr x) (hu : HasDerivAt u du x) :
    HasDerivAt (fun y => N L delta (r y) (u y))
      (Nr (r x) (u x)*dr+Nu L delta (r x) (u x)*du) x := by
  have hv := (hasDerivAt_const x (1 : ℝ)).sub hu
  have hp := ((hasDerivAt_const x delta).add (hu.const_mul (6*delta))).add
    ((hu.pow 2).const_mul (9*delta-2/5))
  have hn := (((hr.pow 2).mul (hu.pow 2)).mul hv).sub
    ((((hv.pow 2).mul hp).const_mul L).div_const 3)
  convert hn using 1 <;> dsimp [N,Nr,Nu,P,Pu] <;>
    (try funext y) <;> (try dsimp) <;> ring

theorem logFRate_hasDerivAt (u : ℝ → ℝ) (L delta r : ℝ)
    (hr : r ≠ 0) (hu0 : u r ≠ 0)
    (hu : HasDerivAt u (H L delta r (u r)) r) :
    HasDerivAt (fun x => G L delta x (u x)) (Gdot L delta r (u r)) r := by
  have hn := N_hasDerivAt (fun x : ℝ => x) u L delta r 1
    (H L delta r (u r)) (hasDerivAt_id r) hu
  have hd := ((hasDerivAt_id r).pow 3).mul (hu.pow 3)
  convert hn.div hd (by dsimp; positivity) using 1 <;>
    dsimp [G,Gdot,H] <;> (try funext x) <;> (try dsimp) <;> field_simp <;> ring

theorem local_curvature_of_solution (u : ℝ → ℝ) (L delta r : ℝ)
    (hr : r ≠ 0) (hu0 : u r ≠ 0)
    (hu : HasDerivAt u (H L delta r (u r)) r) :
    4*((u r*deriv (fun x => G L delta x (u x)) r/2+
        u r*(G L delta r (u r))^2/4+deriv u r*G L delta r (u r)/4)^2+
      2*(u r*G L delta r (u r)/(2*r))^2+2*(deriv u r/(2*r))^2+
      ((1-u r)/r^2)^2)=curvature L delta r (u r) := by
  rw [(logFRate_hasDerivAt u L delta r hr hu0 hu).deriv,hu.deriv]
  rfl

theorem inverse_radius_square_hasDerivAt (r : ℝ → ℝ) (L delta u : ℝ)
    (hr : r u ≠ 0) (hD : D L delta (r u) u ≠ 0)
    (hode : HasDerivAt r ((r u)^3*u^2/D L delta (r u) u) u) :
    HasDerivAt (fun x => 1/(r x)^2) (-2*u^2/D L delta (r u) u) u := by
  convert (hasDerivAt_const u (1 : ℝ)).div (hode.pow 2) (by dsimp; exact pow_ne_zero 2 hr) using 1 <;>
    dsimp <;> field_simp <;> ring

theorem inverse_radius_rate_lower_bound (L delta r u : ℝ)
    (hL : 0 < L) (hd : 0 < delta) (hu : 0 < u) (hu1 : u < 1) :
    -2*u^2/(L*delta*(1-u)^2) ≤ -2*u^2/D L delta r u := by
  have hp := denominator_positive L delta r u hL hd hu hu1
  have hl := denominator_lower_bound L delta r u (le_of_lt hL) (le_of_lt hd)
    (le_of_lt hu) (le_of_lt hu1)
  have hv : 0 < 1-u := sub_pos.mpr hu1
  have hlp : 0 < L*delta*(1-u)^2 := by positivity
  apply (div_le_div_iff₀ hlp hp).2
  nlinarith [mul_nonneg (sq_nonneg u) (sub_nonneg.mpr hl)]

theorem regularized_lapse_hasDerivAt (f : ℝ → ℝ) (L delta r u : ℝ)
    (hu : u ≠ 0) (hD : D L delta r u ≠ 0) (hf : f u ≠ 0)
    (hode : HasDerivAt f (f u*N L delta r u/(u*D L delta r u)) u) :
    HasDerivAt (fun x => Real.log (f x)+Real.log x/3) (lapseRate L delta r u) u := by
  have h1 := (Real.hasDerivAt_log hf).comp u hode
  have h2 := (Real.hasDerivAt_log hu).div_const 3
  convert h1.add h2 using 1 <;> dsimp [lapseRate] <;>
    (try funext x) <;> (try dsimp) <;> field_simp <;> unfold D N P <;> ring

theorem curvature_normalization (L delta r u : ℝ) (hr : r ≠ 0) (hu : u ≠ 0) :
    u^10*curvature L delta r u=normK L delta r u := by
  unfold curvature radialTime Gdot G H normK normA
  field_simp
  <;> ring

theorem denominator_tendsto (r : ℝ → ℝ) (L delta rD : ℝ)
    (hr : Tendsto r (𝓝[>] (0 : ℝ)) (𝓝 rD)) :
    Tendsto (fun u => D L delta (r u) u) (𝓝[>] (0 : ℝ)) (𝓝 (L*delta)) := by
  have hc : ContinuousAt (fun p : ℝ × ℝ => D L delta p.1 p.2) (rD,0) := by
    unfold D
    fun_prop
  have hi : Tendsto (fun u : ℝ => u) (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
    tendsto_id.mono_left inf_le_left
  have hh := hc.tendsto.comp (hr.prodMk_nhds hi)
  simpa [D] using hh

theorem radius_cubic_endpoint_limit (r : ℝ → ℝ) (L delta rD : ℝ)
    (hL : 0 < L) (hd : 0 < delta)
    (hr : Tendsto r (𝓝[>] (0 : ℝ)) (𝓝 rD))
    (hode : ∀ᶠ u in 𝓝[>] (0 : ℝ),
      HasDerivAt r ((r u)^3*u^2/D L delta (r u) u) u) :
    Tendsto (fun u => (r u-rD)/u^3) (𝓝[>] (0 : ℝ)) (𝓝 (rD^3/(3*L*delta))) := by
  have hD := denominator_tendsto r L delta rD hr
  have hn : Tendsto (fun u => r u-rD) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hr.sub_const rD
  have hi : Tendsto (fun u : ℝ => u) (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
    tendsto_id.mono_left inf_le_left
  have hz : Tendsto (fun u : ℝ => u^3) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hi.pow 3
  have hdn : ∀ᶠ u in 𝓝[>] (0 : ℝ),
      HasDerivAt (fun x => r x-rD) ((r u)^3*u^2/D L delta (r u) u) u := by
    filter_upwards [hode] with u hu
    exact hu.sub_const rD
  have hdd : ∀ᶠ u in 𝓝[>] (0 : ℝ), HasDerivAt (fun x : ℝ => x^3) (3*u^2) u :=
    Eventually.of_forall (fun u => by simpa using (hasDerivAt_id u).pow 3)
  have hnonzero : ∀ᶠ u in 𝓝[>] (0 : ℝ), 3*u^2 ≠ 0 := by
    filter_upwards [self_mem_nhdsWithin] with u hu
    have hu0 : 0 < u := hu
    positivity
  have hrat : Tendsto (fun u => (r u)^3/(3*D L delta (r u) u))
      (𝓝[>] (0 : ℝ)) (𝓝 (rD^3/(3*L*delta))) := by
    convert (hr.pow 3).div (hD.const_mul 3) (by positivity) using 1 <;> ring
  have hratio : Tendsto (fun u => ((r u)^3*u^2/D L delta (r u) u)/(3*u^2))
      (𝓝[>] (0 : ℝ)) (𝓝 (rD^3/(3*L*delta))) := by
    apply hrat.congr'
    filter_upwards [self_mem_nhdsWithin] with u hu
    have hu0 : u ≠ 0 := ne_of_gt hu
    field_simp
  exact HasDerivAt.lhopital_zero_nhdsGT hdn hdd hnonzero hn hz hratio

theorem normalized_curvature_tendsto (r : ℝ → ℝ) (L delta rD : ℝ)
    (hrD : 0 < rD) (hr : Tendsto r (𝓝[>] (0 : ℝ)) (𝓝 rD)) :
    Tendsto (fun u => normK L delta (r u) u) (𝓝[>] (0 : ℝ))
      (𝓝 (64*(L*delta)^4/(81*rD^12))) := by
  have hc : ContinuousAt (fun p : ℝ × ℝ => normK L delta p.1 p.2) (rD,0) := by
    unfold normK normA N Nr Nu P Pu D
    fun_prop (disch := positivity)
  have hi : Tendsto (fun u : ℝ => u) (𝓝[>] (0 : ℝ)) (𝓝 (0 : ℝ)) :=
    tendsto_id.mono_left inf_le_left
  have hh := hc.tendsto.comp (hr.prodMk_nhds hi)
  have he : normK L delta rD 0=64*(L*delta)^4/(81*rD^12) := by
    unfold normK normA N Nr Nu P Pu D
    have hr0 := ne_of_gt hrD
    field_simp
    <;> ring
  simpa only [he] using hh

theorem conditional_curvature_diverges (r : ℝ → ℝ) (L delta rD : ℝ)
    (hL : 0 < L) (hd : 0 < delta) (hrD : 0 < rD)
    (hr : Tendsto r (𝓝[>] (0 : ℝ)) (𝓝 rD)) :
    Tendsto (fun u => curvature L delta (r u) u) (𝓝[>] (0 : ℝ)) atTop := by
  have hn := normalized_curvature_tendsto r L delta rD hrD hr
  have hp : 0 < 64*(L*delta)^4/(81*rD^12) := by positivity
  have hi : Tendsto (fun u : ℝ => (u⁻¹)^10) (𝓝[>] (0 : ℝ)) atTop :=
    (tendsto_pow_atTop (by norm_num : (10 : ℕ) ≠ 0)).comp tendsto_inv_nhdsGT_zero
  apply (Tendsto.pos_mul_atTop hp hn hi).congr'
  have hpos := hr.eventually (eventually_gt_nhds (by linarith : (0 : ℝ) < rD))
  filter_upwards [hpos, self_mem_nhdsWithin] with u hru hu
  have hu0 : u ≠ 0 := ne_of_gt hu
  have hh := curvature_normalization L delta (r u) u (ne_of_gt hru) hu0
  rw [← hh]
  field_simp

end
end FourDimEndpoint
