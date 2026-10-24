"""Massive Bohm-Dirac packet on fixed Schwarzschild PG geometry.

Exact diagonal characteristic transport and local unitary mass rotations,
Strang splitting, fifth-degree Lagrange interpolation. kappa=0 is the 1+1
benchmark; kappa=1 is the genuine j=1/2, m_j=1/2 spherical Dirac sector.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp, cumulative_simpson, simpson
from scipy.special import beta
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline, PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
M=1.;MU=1.;P0=-2.;RC=12.;WIDTH=3.
LEFT=.5;RIGHT=40.;END=14.

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def initial(r,mu=MU,p0=P0):
    z=(np.asarray(r)-RC)/WIDTH
    amp=np.where(abs(z)<1,np.maximum(1-z*z,0)**4,0)/np.sqrt(WIDTH*beta(.5,9))
    E=np.sqrt(p0*p0+mu*mu)
    spin=np.array([np.sqrt((E+p0)/(2*E)),np.sqrt((E-p0)/(2*E))])
    return amp[:,None]*np.exp(1j*p0*(r-RC))[:,None]*spin[None,:]

def shift(r,M=M):return np.sqrt(2*M/np.asarray(r))

def current(psi,r,M=M):
    rho=abs(psi)**2;v=shift(r,M)
    q=rho.sum(axis=1);j=(1-v)*rho[:,0]-(1+v)*rho[:,1]
    return q,j

def characteristic_maps(r,dt,M=M):
    """Backward flow and its actual spatial Jacobian, for each chirality."""
    n=len(r);outputs=[]
    for sign in [1.,-1.]:
        def rhs(t,y):
            rr=y[:n];v=shift(rr,M)
            return np.r_[v-sign,-v/(2*rr)]
        sol=solve_ivp(rhs,(0,dt),np.r_[r,np.zeros(n)],method='DOP853',rtol=2e-13,atol=2e-14)
        assert sol.success
        feet=sol.y[:n,-1];factor=np.exp(sol.y[n:,-1]/2)
        # Fixed degree-five local interpolation; one-sided near boundaries.
        zz=(feet-r[0])/(r[1]-r[0]);start=np.floor(zz).astype(int)-2
        start=np.clip(start,0,n-6);idx=start[:,None]+np.arange(6)[None,:]
        u=zz-start;weights=np.ones((n,6))
        for k in range(6):
            for l in range(6):
                if k!=l:weights[:,k]*=(u-l)/(k-l)
        valid=(feet>=r[0]-1e-10)&(feet<=r[-1])
        # Inner boundary is an outflow boundary for both components.
        if M>0:assert np.min(feet-r[0])>=-1e-10
        outputs.append((idx,weights,factor,valid))
    return outputs

def rotation(r,dt,mu,kappa):
    """exp(-i dt B/2), B=mu sigma1-kappa/r sigma2 (Hermitian)."""
    w=kappa/r;omega=np.sqrt(mu*mu+w*w)
    c=np.cos(omega*dt/2);s=np.sinc(omega*dt/(2*np.pi))*dt/2
    return c,s*(w-1j*mu),s*(-w-1j*mu)


def step(psi,maps,rot):
    c,plus,minus=rot
    rotated=np.column_stack([c*psi[:,0]+plus*psi[:,1],c*psi[:,1]+minus*psi[:,0]])
    streamed=np.empty_like(psi)
    for k,(idx,weights,factor,valid) in enumerate(maps):
        streamed[:,k]=np.sum(rotated[idx,k]*weights,axis=1)*factor
        streamed[~valid,k]=0
    return np.column_stack([c*streamed[:,0]+plus*streamed[:,1],c*streamed[:,1]+minus*streamed[:,0]])

def evolve(n=3951,dt=.005,mu=MU,M=M,kappa=0,end=END,save_every=.1,left=LEFT,right=RIGHT):
    nstep=round(end/dt);dt=end/nstep
    stride=round(save_every/dt);assert stride>=1 and abs(stride*dt-save_every)<1e-12
    r=np.linspace(left,right,n);psi=initial(r,mu);maps=characteristic_maps(r,dt,M);rot=rotation(r,dt,mu,kappa)
    steps=[0];states=[psi.copy()];boundary=[];mass=[]
    q,j=current(psi,r,M);boundary.append([j[0],j[-1]]);mass.append(simpson(q,x=r))
    for i in range(1,nstep+1):
        psi=step(psi,maps,rot)
        q,j=current(psi,r,M);boundary.append([j[0],j[-1]])
        if i%stride==0 or i==nstep:
            steps.append(i);states.append(psi.copy());mass.append(simpson(q,x=r))
    steps=np.array(steps);states=np.array(states);boundary=np.array(boundary)
    times=np.arange(nstep+1)*dt
    # Integrate the physical boundary currents, independently of volume norm.
    absorbed=cumulative_simpson(-boundary[:,0],x=times,initial=0)
    outgoing=cumulative_simpson(boundary[:,1],x=times,initial=0)
    return {'r':r,'times':steps*dt,'psi':states,'absorbed':absorbed[steps],
      'outer_loss':outgoing[steps],'norm':np.array(mass),'step_times':times,
      'boundary_current':boundary,'absorbed_fine':absorbed,'outer_loss_fine':outgoing,
      'n':n,'dt':dt,'mu':mu,'M':M,'kappa':kappa,'norm_balance':np.array(mass)+absorbed[steps]+outgoing[steps]}

def quantile_curves(run,labels):
    r=run['r'];q=abs(run['psi'])**2;q=q.sum(axis=2)
    cdf=cumulative_simpson(q,x=r,axis=1,initial=0)+run['absorbed'][:,None]
    curves=np.full((len(labels),len(run['times'])),np.nan)
    for it in range(len(run['times'])):
        keep=np.r_[True,np.diff(cdf[it])>1e-16]
        inverse=PchipInterpolator(cdf[it,keep],r[keep],extrapolate=False)
        for iq,lab in enumerate(labels):
            if run['absorbed'][it]<lab<cdf[it,-1]:curves[iq,it]=inverse(lab)
    return curves,cdf

def compare(a,b):
    assert np.max(abs(a['times']-b['times']))<1e-12
    errors=[];density=[]
    for pa,pb in zip(a['psi'],b['psi']):
        interp=CubicSpline(a['r'],pa,axis=0)(b['r']) if len(a['r'])!=len(b['r']) else pa
        errors.append(np.sqrt(simpson(np.sum(abs(interp-pb)**2,axis=1),x=b['r'])))
        density.append(simpson(abs(np.sum(abs(interp)**2,axis=1)-np.sum(abs(pb)**2,axis=1)),x=b['r']))
    return {'max_L2_spinor':float(max(errors)),'max_L1_density':float(max(density)),
      'max_absorbed':float(max(abs(a['absorbed']-b['absorbed'])))}

def primitive(r,M=M):
    r=np.asarray(r);x=np.sqrt(r/(2*M))
    return 4*M*(x*x/2-x+np.log1p(x))

def massless_reference(r,t,M=M):
    # Unique ingoing map g(r_initial)=g(r)+t, with input packet as above.
    old=np.array([brentq(lambda rr:primitive(rr,M)-primitive(x,M)-t,x,x+t*(1+np.sqrt(2*M/x))+1,xtol=2e-13) if t>0 else x for x in r])
    out=initial(old,mu=0)
    out[:,1]*=np.sqrt((1+shift(old,M))/(1+shift(r,M)))
    return out

def geodesics(starts,times,mu=MU):
    rows=[]
    for start in starts:
        def fun(t,y):
            r,p=y;v=shift(r);E=np.sqrt(p*p+mu*mu)
            return [p/E-v,-v*p/(2*r)]
        def stop(t,y):return y[0]-LEFT
        stop.terminal=True;stop.direction=-1
        sol=solve_ivp(fun,(0,times[-1]),[start,P0],method='DOP853',events=stop,dense_output=True,rtol=2e-12,atol=2e-13,max_step=.05)
        valid=times<=sol.t[-1];data=np.full((2,len(times)),np.nan);data[:,valid]=sol.sol(times[valid])
        en=np.sqrt(sol.y[1]**2+mu*mu)-shift(sol.y[0])*sol.y[1]
        rows.append({'r':data[0],'p':data[1],'energy_drift':float(max(abs(en-en[0]))),
          'cutoff_time':float(sol.t_events[0][0]) if len(sol.t_events[0]) else None})
    return rows

def crossing(times,radius,level):
    valid=np.isfinite(radius);t=times[valid];r=radius[valid]
    hits=np.flatnonzero((r[:-1]>=level)&(r[1:]<level))
    if not len(hits):return None
    k=hits[0];return float(t[k]+(t[k+1]-t[k])*(r[k]-level)/(r[k]-r[k+1]))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tag',required=True);ap.add_argument('--pilot',action='store_true')
    args=ap.parse_args();assert args.tag.replace('_','').replace('-','').isalnum()
    out=ROOT/f'work/massive_bohm_{args.tag}.json';assert not out.exists()
    begin=time.perf_counter()
    if args.pilot:
        run=evolve(n=1976,dt=.01,save_every=.1)
        curves,cdf=quantile_curves(run,[.1,.25,.5,.75,.9])
        print(json.dumps({'seconds':time.perf_counter()-begin,'initial_norm':run['norm'][0],
          'norm_balance_error':float(max(abs(run['norm_balance']-1))),
          'absorbed':run['absorbed'][-1],'surviving':run['norm'][-1],
          'horizon_times':[crossing(run['times'],row,2.) for row in curves]},indent=2))
        np.savez_compressed(out.with_suffix('.npz'),**{k:v for k,v in run.items() if isinstance(v,np.ndarray)},curves=curves)
        return
    from massive_suite import run_suite
    run_suite(args.tag,out,begin)

if __name__=='__main__':main()
