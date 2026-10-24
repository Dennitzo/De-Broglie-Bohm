import Mathlib

/-! Local stationary Dirac energy and momentum with explicit field derivatives. -/
namespace FieldCurvEnergy
noncomputable section
abbrev Spinor := (ℝ × ℝ) × (ℝ × ℝ)
theorem d11 {u : ℝ → Spinor} {du : Spinor} {x : ℝ} (h : HasDerivAt u du x) :
    HasDerivAt (fun t => (u t).1.1) du.1.1 x := by
  exact hasFDerivAt_fst.comp_hasDerivAt x (hasFDerivAt_fst.comp_hasDerivAt x h)

theorem d12 {u : ℝ → Spinor} {du : Spinor} {x : ℝ} (h : HasDerivAt u du x) :
    HasDerivAt (fun t => (u t).1.2) du.1.2 x := by
  exact hasFDerivAt_snd.comp_hasDerivAt x (hasFDerivAt_fst.comp_hasDerivAt x h)

theorem d21 {u : ℝ → Spinor} {du : Spinor} {x : ℝ} (h : HasDerivAt u du x) :
    HasDerivAt (fun t => (u t).2.1) du.2.1 x := by
  exact hasFDerivAt_fst.comp_hasDerivAt x (hasFDerivAt_snd.comp_hasDerivAt x h)

theorem d22 {u : ℝ → Spinor} {du : Spinor} {x : ℝ} (h : HasDerivAt u du x) :
    HasDerivAt (fun t => (u t).2.2) du.2.2 x := by
  exact hasFDerivAt_snd.comp_hasDerivAt x (hasFDerivAt_snd.comp_hasDerivAt x h)

def evolve (v vp mu w : ℝ) (u ux : Spinor) : Spinor :=
  (((v-1)*ux.1.1+vp/2*u.1.1+mu*u.2.2+w*u.2.1,
    (v-1)*ux.1.2+vp/2*u.1.2-mu*u.2.1+w*u.2.2),
   ((v+1)*ux.2.1+vp/2*u.2.1+mu*u.1.2-w*u.1.1,
    (v+1)*ux.2.2+vp/2*u.2.2-mu*u.1.1-w*u.1.2))
def spatialEvolve (v vp vpp mu w wp : ℝ) (u ux uxx : Spinor) : Spinor :=
  (((v-1)*uxx.1.1+3*vp/2*ux.1.1+vpp/2*u.1.1+mu*ux.2.2+w*ux.2.1+wp*u.2.1,
    (v-1)*uxx.1.2+3*vp/2*ux.1.2+vpp/2*u.1.2-mu*ux.2.1+w*ux.2.2+wp*u.2.2),
   ((v+1)*uxx.2.1+3*vp/2*ux.2.1+vpp/2*u.2.1+mu*ux.1.2-w*ux.1.1-wp*u.1.1,
    (v+1)*uxx.2.2+3*vp/2*ux.2.2+vpp/2*u.2.2-mu*ux.1.1-w*ux.1.2-wp*u.1.2))
def dot (u z : Spinor) : ℝ := u.1.1*z.1.1+u.1.2*z.1.2+u.2.1*z.2.1+u.2.2*z.2.2
def cDot (v : ℝ) (u z : Spinor) : ℝ :=
  (1-v)*(u.1.1*z.1.1+u.1.2*z.1.2)+(-1-v)*(u.2.1*z.2.1+u.2.2*z.2.2)
def iSpin (u : Spinor) : Spinor := ((-u.1.2,u.1.1),(-u.2.2,u.2.1))
def momentum (u ux : Spinor) : ℝ := u.1.1*ux.1.2-u.1.2*ux.1.1+u.2.1*ux.2.2-u.2.2*ux.2.1
def radialKinetic (u ux : Spinor) : ℝ := u.1.1*ux.1.2-u.1.2*ux.1.1-u.2.1*ux.2.2+u.2.2*ux.2.1
def interaction (mu w : ℝ) (u : Spinor) : ℝ :=
  2*mu*(u.1.1*u.2.1+u.1.2*u.2.2)-2*w*(u.1.1*u.2.2-u.1.2*u.2.1)
def energy (v mu w : ℝ) (u ux : Spinor) : ℝ := radialKinetic u ux+interaction mu w u-v*momentum u ux
def energyFlux (v mu w : ℝ) (u ux : Spinor) : ℝ :=
  (1+v*v)*momentum u ux-v*(2*radialKinetic u ux+interaction mu w u)
def momentumFlux (v : ℝ) (u ux : Spinor) : ℝ := radialKinetic u ux-v*momentum u ux

theorem energy_time_identity (v vp mu w : ℝ) (u ux : Spinor) :
    dot u (iSpin (evolve v vp mu w u ux))=energy v mu w u ux := by
  dsimp [dot,iSpin,evolve,energy,radialKinetic,interaction,momentum]
  ring

theorem energy_flux_time_identity (v vp mu w : ℝ) (u ux : Spinor) :
    cDot v u (iSpin (evolve v vp mu w u ux))=energyFlux v mu w u ux := by
  dsimp [cDot,iSpin,evolve,energyFlux,radialKinetic,interaction,momentum]
  ring

theorem evolve_commutes_i (v vp mu w : ℝ) (u ux : Spinor) :
    evolve v vp mu w (iSpin u) (iSpin ux)=iSpin (evolve v vp mu w u ux) := by
  ext <;> dsimp [evolve,iSpin] <;> ring

theorem polarized_continuity (u z : ℝ → ℝ → Spinor) (v : ℝ → ℝ)
    (t r vp mu w : ℝ) (ur zr : Spinor)
    (hv : HasDerivAt v vp r)
    (huR : HasDerivAt (u t) ur r) (hzR : HasDerivAt (z t) zr r)
    (huT : HasDerivAt (fun s => u s r) (evolve (v r) vp mu w (u t r) ur) t)
    (hzT : HasDerivAt (fun s => z s r) (evolve (v r) vp mu w (z t r) zr) t) :
    deriv (fun s => dot (u s r) (z s r)) t+
      deriv (fun x => cDot (v x) (u t x) (z t x)) r=0 := by
  have hT := (((((d11 huT)).mul ((d11 hzT))).add (((d12 huT)).mul ((d12 hzT)))).add (((d21 huT)).mul ((d21 hzT)))).add (((d22 huT)).mul ((d22 hzT)))
  have hR := ((((hasDerivAt_const r (1 : ℝ))).sub (hv)).mul ((((d11 huR)).mul ((d11 hzR))).add (((d12 huR)).mul ((d12 hzR))))).add ((((hasDerivAt_const r (-1 : ℝ))).sub (hv)).mul ((((d21 huR)).mul ((d21 hzR))).add (((d22 huR)).mul ((d22 hzR)))))
  change HasDerivAt (fun s => dot (u s r) (z s r)) _ t at hT
  change HasDerivAt (fun x => cDot (v x) (u t x) (z t x)) _ r at hR
  rw [hT.deriv,hR.deriv]
  dsimp [evolve]
  ring

theorem stationary_energy_conservation (u ux : ℝ → ℝ → Spinor) (v w : ℝ → ℝ)
    (t r vp vpp mu wp : ℝ) (uxx : Spinor)
    (hv : HasDerivAt v vp r) (hw : HasDerivAt w wp r)
    (huR : HasDerivAt (u t) (ux t r) r) (hxR : HasDerivAt (ux t) uxx r)
    (huT : HasDerivAt (fun s => u s r) (evolve (v r) vp mu (w r) (u t r) (ux t r)) t)
    (hxT : HasDerivAt (fun s => ux s r) (spatialEvolve (v r) vp vpp mu (w r) wp (u t r) (ux t r) uxx) t) :
    deriv (fun s => energy (v r) mu (w r) (u s r) (ux s r)) t+
      deriv (fun x => energyFlux (v x) mu (w x) (u t x) (ux t x)) r=0 := by
  have hT := (((((((d11 huT)).mul ((d12 hxT))).sub (((d12 huT)).mul ((d11 hxT)))).sub (((d21 huT)).mul ((d22 hxT)))).add (((d22 huT)).mul ((d21 hxT)))).add (((((hasDerivAt_const t (2 : ℝ))).mul ((hasDerivAt_const t mu))).mul ((((d11 huT)).mul ((d21 huT))).add (((d12 huT)).mul ((d22 huT))))).sub ((((hasDerivAt_const t (2 : ℝ))).mul ((hasDerivAt_const t (w r)))).mul ((((d11 huT)).mul ((d22 huT))).sub (((d12 huT)).mul ((d21 huT))))))).sub (((hasDerivAt_const t (v r))).mul ((((((d11 huT)).mul ((d12 hxT))).sub (((d12 huT)).mul ((d11 hxT)))).add (((d21 huT)).mul ((d22 hxT)))).sub (((d22 huT)).mul ((d21 hxT)))))
  have hR := ((((hasDerivAt_const r (1 : ℝ))).add ((hv).mul (hv))).mul ((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).add (((d21 huR)).mul ((d22 hxR)))).sub (((d22 huR)).mul ((d21 hxR))))).sub ((hv).mul ((((hasDerivAt_const r (2 : ℝ))).mul ((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).sub (((d21 huR)).mul ((d22 hxR)))).add (((d22 huR)).mul ((d21 hxR))))).add (((((hasDerivAt_const r (2 : ℝ))).mul ((hasDerivAt_const r mu))).mul ((((d11 huR)).mul ((d21 huR))).add (((d12 huR)).mul ((d22 huR))))).sub ((((hasDerivAt_const r (2 : ℝ))).mul (hw)).mul ((((d11 huR)).mul ((d22 huR))).sub (((d12 huR)).mul ((d21 huR))))))))
  change HasDerivAt (fun s => energy (v r) mu (w r) (u s r) (ux s r)) _ t at hT
  change HasDerivAt (fun x => energyFlux (v x) mu (w x) (u t x) (ux t x)) _ r at hR
  rw [hT.deriv,hR.deriv]
  dsimp [evolve,spatialEvolve,momentum]
  ring

theorem stationary_momentum_exchange (u ux : ℝ → ℝ → Spinor) (v w : ℝ → ℝ)
    (t r vp vpp mu wp : ℝ) (uxx : Spinor)
    (hv : HasDerivAt v vp r) (hw : HasDerivAt w wp r)
    (huR : HasDerivAt (u t) (ux t r) r) (hxR : HasDerivAt (ux t) uxx r)
    (huT : HasDerivAt (fun s => u s r) (evolve (v r) vp mu (w r) (u t r) (ux t r)) t)
    (hxT : HasDerivAt (fun s => ux s r) (spatialEvolve (v r) vp vpp mu (w r) wp (u t r) (ux t r) uxx) t) :
    deriv (fun s => momentum (u s r) (ux s r)) t+
      deriv (fun x => momentumFlux (v x) (u t x) (ux t x)) r=vp*momentum (u t r) (ux t r)+2*wp*((u t r).1.1*(u t r).2.2-(u t r).1.2*(u t r).2.1) := by
  have hT := (((((d11 huT)).mul ((d12 hxT))).sub (((d12 huT)).mul ((d11 hxT)))).add (((d21 huT)).mul ((d22 hxT)))).sub (((d22 huT)).mul ((d21 hxT)))
  have hR := ((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).sub (((d21 huR)).mul ((d22 hxR)))).add (((d22 huR)).mul ((d21 hxR)))).sub ((hv).mul ((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).add (((d21 huR)).mul ((d22 hxR)))).sub (((d22 huR)).mul ((d21 hxR)))))
  change HasDerivAt (fun s => momentum (u s r) (ux s r)) _ t at hT
  change HasDerivAt (fun x => momentumFlux (v x) (u t x) (ux t x)) _ r at hR
  rw [hT.deriv,hR.deriv]
  dsimp [evolve,spatialEvolve,momentum]
  ring

theorem normal_killing_energy_relation (v mu w : ℝ) (u ux : Spinor) :
    energy v mu w u ux+v*momentum u ux=radialKinetic u ux+interaction mu w u := by
  unfold energy
  ring

theorem energy_flux_stress_relation (v mu w : ℝ) (u ux : Spinor) :
    energyFlux v mu w u ux=(1+v^2)*momentum u ux-
      v*((radialKinetic u ux+interaction mu w u)+radialKinetic u ux) := by
  unfold energyFlux
  ring

theorem spatial_evolve_hasDerivAt (u ux : ℝ → Spinor) (v vp w : ℝ → ℝ)
    (r vpp mu wp : ℝ) (uxx : Spinor)
    (hv : HasDerivAt v (vp r) r) (hvp : HasDerivAt vp vpp r)
    (hw : HasDerivAt w wp r) (hu : HasDerivAt u (ux r) r)
    (hx : HasDerivAt ux uxx r) :
    HasDerivAt (fun x => evolve (v x) (vp x) mu (w x) (u x) (ux x))
      (spatialEvolve (v r) (vp r) vpp mu (w r) wp (u r) (ux r) uxx) r := by
  have h0 := (((((hv).sub ((hasDerivAt_const r (1 : ℝ)))).mul ((d11 hx))).add ((((hasDerivAt_const r (1/2 : ℝ))).mul (hvp)).mul ((d11 hu)))).add (((hasDerivAt_const r mu)).mul ((d22 hu)))).add ((hw).mul ((d21 hu)))
  have h1 := (((((hv).sub ((hasDerivAt_const r (1 : ℝ)))).mul ((d12 hx))).add ((((hasDerivAt_const r (1/2 : ℝ))).mul (hvp)).mul ((d12 hu)))).sub (((hasDerivAt_const r mu)).mul ((d21 hu)))).add ((hw).mul ((d22 hu)))
  have h2 := (((((hv).add ((hasDerivAt_const r (1 : ℝ)))).mul ((d21 hx))).add ((((hasDerivAt_const r (1/2 : ℝ))).mul (hvp)).mul ((d21 hu)))).add (((hasDerivAt_const r mu)).mul ((d12 hu)))).sub ((hw).mul ((d11 hu)))
  have h3 := (((((hv).add ((hasDerivAt_const r (1 : ℝ)))).mul ((d22 hx))).add ((((hasDerivAt_const r (1/2 : ℝ))).mul (hvp)).mul ((d22 hu)))).sub (((hasDerivAt_const r mu)).mul ((d11 hu)))).sub ((hw).mul ((d12 hu)))
  convert (h0.prodMk h1).prodMk (h2.prodMk h3) using 1 <;>
    ext <;> dsimp [evolve,spatialEvolve] <;> ring

theorem normal_observer_work_balance (u ux : ℝ → ℝ → Spinor) (v w : ℝ → ℝ)
    (t r vp vpp mu wp : ℝ) (uxx : Spinor)
    (hv : HasDerivAt v vp r) (hw : HasDerivAt w wp r)
    (huR : HasDerivAt (u t) (ux t r) r) (hxR : HasDerivAt (ux t) uxx r)
    (huT : HasDerivAt (fun s => u s r) (evolve (v r) vp mu (w r) (u t r) (ux t r)) t)
    (hxT : HasDerivAt (fun s => ux s r) (spatialEvolve (v r) vp vpp mu (w r) wp (u t r) (ux t r) uxx) t) :
    deriv (fun s => radialKinetic (u s r) (ux s r)+interaction mu (w r) (u s r)) t+
      deriv (fun x => momentum (u t x) (ux t x)-v x*(radialKinetic (u t x) (ux t x)+interaction mu (w x) (u t x))) r=
      vp*radialKinetic (u t r) (ux t r)+2*v r*wp*((u t r).1.1*(u t r).2.2-(u t r).1.2*(u t r).2.1) := by
  have hT := ((((((d11 huT)).mul ((d12 hxT))).sub (((d12 huT)).mul ((d11 hxT)))).sub (((d21 huT)).mul ((d22 hxT)))).add (((d22 huT)).mul ((d21 hxT)))).add (((((hasDerivAt_const t (2 : ℝ))).mul ((hasDerivAt_const t mu))).mul ((((d11 huT)).mul ((d21 huT))).add (((d12 huT)).mul ((d22 huT))))).sub ((((hasDerivAt_const t (2 : ℝ))).mul ((hasDerivAt_const t (w r)))).mul ((((d11 huT)).mul ((d22 huT))).sub (((d12 huT)).mul ((d21 huT))))))
  have hR := ((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).add (((d21 huR)).mul ((d22 hxR)))).sub (((d22 huR)).mul ((d21 hxR)))).sub ((hv).mul (((((((d11 huR)).mul ((d12 hxR))).sub (((d12 huR)).mul ((d11 hxR)))).sub (((d21 huR)).mul ((d22 hxR)))).add (((d22 huR)).mul ((d21 hxR)))).add (((((hasDerivAt_const r (2 : ℝ))).mul ((hasDerivAt_const r mu))).mul ((((d11 huR)).mul ((d21 huR))).add (((d12 huR)).mul ((d22 huR))))).sub ((((hasDerivAt_const r (2 : ℝ))).mul (hw)).mul ((((d11 huR)).mul ((d22 huR))).sub (((d12 huR)).mul ((d21 huR))))))))
  change HasDerivAt (fun s => radialKinetic (u s r) (ux s r)+interaction mu (w r) (u s r)) _ t at hT
  change HasDerivAt (fun x => momentum (u t x) (ux t x)-v x*(radialKinetic (u t x) (ux t x)+interaction mu (w x) (u t x))) _ r at hR
  rw [hT.deriv,hR.deriv]
  dsimp [evolve,spatialEvolve,radialKinetic]
  ring

end
end FieldCurvEnergy
