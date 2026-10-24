"""Massive Dirac guidance on the regular Simpson--Visser PG background.

This is a test-field calculation on a specified regular geometry with an
additional NEC-violating source, not a derived Dirac quantum backreaction.
No previous research artefacts or PDF files are modified.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from scipy.integrate import solve_ivp,simpson,cumulative_simpson
from scipy.interpolate import CubicSpline,PchipInterpolator,RectBivariateSpline
from scipy.special import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
M=1.;A=1.;MU=1.;KAPPA=1.;P0=-2.;RC=12.;WIDTH=3.
LEFT=-35.;RIGHT=34.;END=16.
LABELS=np.array([.1,.25,.5,.75,.9])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def radius(r,a=A):return np.sqrt(np.asarray(r)**2+a*a)
def shift(r,M=M,a=A):return np.sqrt(2*M/radius(r,a))
def shift_prime(r,M=M,a=A):return -shift(r,M,a)*np.asarray(r)/(2*radius(r,a)**2)
def initial(r,mu=MU,p0=P0):
    u=(np.asarray(r)-RC)/WIDTH
    amp=np.maximum(1-u*u,0)**4/np.sqrt(WIDTH*beta(.5,9))
    E=np.sqrt(p0*p0+mu*mu)
    pol=np.array([np.sqrt((E+p0)/(2*E)),np.sqrt((E-p0)/(2*E))])
    return amp[:,None]*np.exp(1j*p0*(r-RC))[:,None]*pol[None,:]
def density_current(p,r,M=M,a=A):
    q=np.sum(abs(p)**2,axis=-1);v=shift(r,M,a)
    j=(1-v)*abs(p[...,0])**2-(1+v)*abs(p[...,1])**2
    return q,j
def positive_cdf(r,p):
    """Positive cell integrals: exact Gauss integration of |cubic spinor|²."""
    nodes,weights=np.polynomial.legendre.leggauss(4)
    mid=(r[1:]+r[:-1])/2;half=np.diff(r)/2
    values=CubicSpline(r,p,axis=0)(mid[:,None]+half[:,None]*nodes)
    cells=half*np.sum(weights*np.sum(abs(values)**2,axis=2),axis=1)
    return np.r_[0.,np.cumsum(cells)]
def maps(r,dt,M=M,a=A):
    n=len(r);out=[]
    for sign in [1.,-1.]:
        def rhs(t,y):return np.r_[shift(y[:n],M,a)-sign,shift_prime(y[:n],M,a)]
        sol=solve_ivp(rhs,(0,dt),np.r_[r,np.zeros(n)],method='DOP853',rtol=2e-13,atol=2e-14)
        assert sol.success
        feet=sol.y[:n,-1];factor=np.exp(sol.y[n:,-1]/2)
        u=(feet-r[0])/(r[1]-r[0]);start=np.clip(np.floor(u).astype(int)-2,0,n-6)
        ix=start[:,None]+np.arange(6);local=u-start;weights=np.ones((n,6))
        for k in range(6):
            for l in range(6):
                if k!=l:weights[:,k]*=(local-l)/(k-l)
        out.append((ix,weights,factor,(feet>=r[0])&(feet<=r[-1])))
    return out
def evolve(dr=.02,dt=.0025,end=END,save_every=.05,mu=MU,kappa=KAPPA,M=M,a=A):
    n=round((RIGHT-LEFT)/dr)+1;r=np.linspace(LEFT,RIGHT,n)
    nstep=round(end/dt);dt=end/nstep;stride=round(save_every/dt)
    assert abs(stride*dt-save_every)<1e-12
    assert RC-WIDTH-(1+np.sqrt(2*M/a))*end>LEFT
    assert RC+WIDTH+end<RIGHT
    mp=maps(r,dt,M,a);w=kappa/radius(r,a);om=np.sqrt(mu*mu+w*w)
    c=np.cos(om*dt/2);s=np.sinc(om*dt/(2*np.pi))*dt/2
    plus=s*(w-1j*mu);minus=s*(-w-1j*mu)
    p=initial(r,mu);states=[p.copy()];times=[0.];norms=[];norms.append(simpson(abs(p[:,0])**2+abs(p[:,1])**2,x=r))
    rh=np.sqrt(4*M*M-a*a);probe=np.array([-rh,0,rh])
    probe_u=(probe-r[0])/(r[1]-r[0]);start=np.floor(probe_u).astype(int)-2
    probe_ix=start[:,None]+np.arange(6);u=probe_u-start;probe_weights=np.ones((3,6))
    for k in range(6):
        for l in range(6):
            if k!=l:probe_weights[:,k]*=(u-l)/(k-l)
    # Physical currents at exact probe radii, with degree-five interpolation.
    flux=[];boundary=[]
    def record_flux(p):
        pp=np.sum(p[probe_ix]*probe_weights[:,:,None],axis=1)
        flux.append(density_current(pp,probe,M,a)[1]);boundary.append(density_current(p[[0,-1]],r[[0,-1]],M,a)[1])
    record_flux(p)
    for step in range(1,nstep+1):
        pp=np.column_stack([c*p[:,0]+plus*p[:,1],c*p[:,1]+minus*p[:,0]])
        streamed=np.empty_like(p)
        for k,(ix,weights,factor,valid) in enumerate(mp):
            streamed[:,k]=np.sum(pp[ix,k]*weights,axis=1)*factor;streamed[~valid,k]=0
        p=np.column_stack([c*streamed[:,0]+plus*streamed[:,1],c*streamed[:,1]+minus*streamed[:,0]])
        record_flux(p)
        if step%stride==0 or step==nstep:
            states.append(p.copy());times.append(step*dt);norms.append(simpson(np.sum(abs(p)**2,axis=1),x=r))
    times=np.array(times);states=np.array(states);q,j=density_current(states,r,M,a)
    fine_times=np.arange(nstep+1)*dt
    flux=np.array(flux);boundary=np.array(boundary)
    crossed=-cumulative_simpson(flux,x=fine_times,axis=0,initial=0)
    boundary_balance=cumulative_simpson(boundary[:,1]-boundary[:,0],x=fine_times,initial=0)
    sampled=np.rint(times/dt).astype(int)
    cdf=np.array([positive_cdf(r,p) for p in states])
    # Vanishing subnormal CDF slopes have the exact PCHIP limit zero.
    with np.errstate(over='ignore'):
        probabilities=np.array([PchipInterpolator(r,cc)(probe) for cc in cdf])
    curves=np.empty((len(LABELS),len(times)))
    for it,cc in enumerate(cdf):
        keep=np.r_[True,np.diff(cc)>1e-16];inv=PchipInterpolator(cc[keep],r[keep],extrapolate=False)
        curves[:,it]=inv(LABELS)
    return dict(r=r,times=times,psi=states,q=q,j=j,norm=np.array(norms),curves=curves,
      probe_r=probe,probabilities=probabilities,flux_integral=crossed[sampled],
      fine_times=fine_times,probe_flux=flux,boundary_flux=boundary,
      boundary_balance=boundary_balance[sampled],dt=dt,dr=(RIGHT-LEFT)/(n-1),M=M,a=a,mu=mu,kappa=kappa)

def direct_guidance(run):
    r=run['r'];times=run['times'];psi=run['psi']
    interp=[(RectBivariateSpline(times,r,psi[:,:,k].real),RectBivariateSpline(times,r,psi[:,:,k].imag)) for k in range(2)]
    def field(t,r):return np.array([p.ev(t,r)+1j*q.ev(t,r) for p,q in interp])
    trajectories=[];events=[];minimum_density=[]
    for r0 in run['curves'][:,0]:
        def rhs(t,y):
            r,phi,tau=y;p=field(t,r);q=np.sum(abs(p)**2);z=p[0].conj()*p[1]
            return [-shift(r)+(abs(p[0])**2-abs(p[1])**2)/q,-2*z.imag/(radius(r)*q),2*abs(z.real)/q]
        ev=[]
        for level in [np.sqrt(3),0,-np.sqrt(3)]:
            def hit(t,y,level=level):return y[0]-level
            hit.direction=-1;hit.terminal=False;ev.append(hit)
        sol=solve_ivp(rhs,(times[0],times[-1]),[r0,0,0],method='DOP853',events=ev,dense_output=True,rtol=2e-10,atol=2e-11,max_step=.025)
        assert sol.success and all(len(v)==1 for v in sol.t_events)
        path=sol.sol(times);trajectories.append(path)
        events.append([{'T':float(t[0]),'r':float(y[0,0]),'phi':float(y[0,1]),'tau':float(y[0,2])} for t,y in zip(sol.t_events,sol.y_events)])
        minimum_density.append(float(min(np.sum(abs(field(t,x))**2) for t,x in zip(times,path[0]))))
    return np.array(trajectories),events,minimum_density

def compare(a,b):
    assert np.max(abs(a['times']-b['times']))<1e-12
    errors=[];dens=[]
    for pa,pb in zip(a['psi'],b['psi']):
        p=CubicSpline(a['r'],pa,axis=0)(b['r']) if len(a['r'])!=len(b['r']) else pa
        errors.append(np.sqrt(simpson(np.sum(abs(p-pb)**2,axis=1),x=b['r'])))
        dens.append(simpson(abs(np.sum(abs(p)**2,axis=1)-np.sum(abs(pb)**2,axis=1)),x=b['r']))
    return {'max_spinor_L2':float(max(errors)),'max_density_L1':float(max(dens)),
      'max_quantile_difference':float(np.max(abs(a['curves']-b['curves'])))}

def plots(run,paths,tag):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
    colors=plt.cm.viridis(np.linspace(.05,.9,len(LABELS)));files={}
    def save(fig,name):
        for ext in ['png','svg']:
            p=ROOT/f'figures/nonsing_{tag}_{name}.{ext}';assert not p.exists();fig.savefig(p,bbox_inches='tight');files[p.relative_to(ROOT).as_posix()]=sha(p)
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    keep=(run['r']>-32)&(run['r']<17)
    im=ax[0].pcolormesh(run['times'],run['r'][keep],run['q'][:,keep].T,shading='auto',cmap='magma',vmin=0,vmax=.6,rasterized=True)
    fig.colorbar(im,ax=ax[0],label='Radiale Wahrscheinlichkeit q')
    for y in [-np.sqrt(3),np.sqrt(3)]:ax[0].axhline(y,color='cyan',ls='--',lw=.8)
    ax[0].axhline(0,color='white',ls=':',lw=.8)
    for label,c,path in zip(LABELS,colors,paths):
        ax[0].plot(run['times'],path[0],color=c,lw=1)
        ax[1].plot(run['times'],radius(path[0]),color=c,label=f'ζ = {label:g}')
    ax[0].set(xlabel='PG-Zeit T / M',ylabel='Signierte Koordinate r / M',title='Diracpaket durch beide Horizonte')
    ax[1].axhline(1,color='black',ls=':',label='Minimalradius a / M = 1')
    ax[1].set(xlabel='PG-Zeit T / M',ylabel='Flächenradius R / M',title='R wächst nach der Halsdurchquerung');ax[1].legend(fontsize=8)
    fig.suptitle('Massive Bohmführung auf regulärem Black-bounce-Hintergrund: μM = a/M = κ = 1')
    save(fig,'transport')
    fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    prob=run['probabilities'];t=run['times']
    axs[0].plot(t,1-prob[:,2],label='Rechtes Außengebiet');axs[0].plot(t,prob[:,2]-prob[:,0],label='Innenbereich');axs[0].plot(t,prob[:,0],label='Linkes Außengebiet')
    axs[0].plot(t,run['norm'],'k:',label='Gesamtnorm');axs[0].legend(fontsize=8);axs[0].set(xlabel='T / M',ylabel='Wahrscheinlichkeit',title='Kein Verlust am Minimalradius')
    for label,c,path in zip(LABELS,colors,paths):
        axs[1].plot(t,path[1],color=c,label=f'ζ = {label:g}')
        axs[2].plot(t,path[2],color=c)
    axs[1].set(xlabel='T / M',ylabel='Azimut φ / rad',title='Winkelstrom läuft glatt weiter');axs[1].legend(fontsize=8)
    axs[2].plot(t,t,'k:',label='τ = T');axs[2].set(xlabel='T / M',ylabel='Bahn-Eigenzeit τ / M',title='Vollständige kausale Führung');axs[2].legend(fontsize=8)
    save(fig,'fate')
    fig=plt.figure(figsize=(12,5.8),layout='constrained')
    for panel,limit in [(1,16),(2,4)]:
        ax=fig.add_subplot(1,2,panel,projection='3d')
        rr,ph=np.meshgrid(np.linspace(-limit,limit,110),np.linspace(0,2*np.pi,70))
        R=radius(rr);Z=np.arcsinh(rr)
        ax.plot_surface(R*np.cos(ph),R*np.sin(ph),Z,color='#86afc5',alpha=.12,linewidth=0,rstride=3,cstride=3)
        angle=np.linspace(0,2*np.pi,150)
        for rv,col,sty in [(0,'black','-'),(np.sqrt(3),'#5274a2','--'),(-np.sqrt(3),'#5274a2','--')]:
            ax.plot(radius(rv)*np.cos(angle),radius(rv)*np.sin(angle),np.full_like(angle,np.arcsinh(rv)),color=col,ls=sty,lw=.8)
        for i,(path,c) in enumerate(zip(paths,colors)):
            r,phi=path[:2];ph=phi+.9*i;sel=abs(r)<=limit
            ax.plot(radius(r[sel])*np.cos(ph[sel]),radius(r[sel])*np.sin(ph[sel]),np.arcsinh(r[sel]),color=c,lw=1.8,label=f'ζ = {LABELS[i]:g}')
        ax.set(xlabel='X / M',ylabel='Y / M',zlabel='Z / M',title='Äquatoriale Raumschnitt-Einbettung' if panel==1 else 'Glatte Passage durch R = a')
        ax.set_box_aspect((1,1,.65));ax.view_init(elev=18,azim=-65)
        if panel==1:ax.legend(fontsize=8,loc='upper left')
    fig.suptitle('X = R cos φ, Y = R sin φ, Z = a arsinh(r/a); keine Raumzeit-Einbettung')
    save(fig,'trajectories3d')
    return files

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tag',required=True);ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
    assert args.tag.replace('_','').replace('-','').isalnum();out=ROOT/f'work/nonsing_bohm_{args.tag}.json';assert not out.exists()
    begin=time.perf_counter()
    if args.pilot:
        run=evolve(dr=.05,dt=.005,save_every=.1);paths,events,density=direct_guidance(run)
        print(json.dumps({'seconds':time.perf_counter()-begin,'norm_error':float(max(abs(run['norm']-1))),
          'events':events,'end_probabilities':run['probabilities'][-1].tolist(),
          'minimum_path_density':density,'guidance_quantile_error':float(np.max(abs(paths[:,0]-run['curves'])))},indent=2),flush=True)
        out.write_text(json.dumps({'pilot':True,'source_sha256':sha(__file__),'events':events},indent=2),encoding='utf8')
        return
    print('Time refinements',flush=True)
    coarse=evolve(dr=.005,dt=.0025);medium=evolve(dr=.005,dt=.00125);run=evolve(dr=.005,dt=.000625)
    conv={'time_coarse_medium':compare(coarse,medium),'time_medium_fine':compare(medium,run)}
    del coarse,medium
    print('Space refinement',flush=True)
    fine=evolve(dr=.0025,dt=.000625);conv['space']=compare(run,fine);del fine
    print('Independent guidance on stored spinors',flush=True)
    paths,events,pathdensity=direct_guidance(run)
    norms=run['norm'];balance=norms+run['boundary_balance']
    fluxerr=float(np.max(abs(run['probabilities']-run['probabilities'][0]-run['flux_integral'])))
    radialerr=float(np.max(abs(paths[:,0]-run['curves'])))
    checks={'unit_norm':max(abs(norms-1))<2e-6,'boundary_balance':max(abs(balance-1))<2e-6,
      'internal_flux_balance':fluxerr<5e-6,'all_paths_cross_both_horizons':all(len(v)==3 for v in events),
      'quantile_guidance_match':radialerr<2e-4,'positive_path_densities':min(pathdensity)>1e-5,
      'time_refinement':conv['time_medium_fine']['max_spinor_L2']<.35*conv['time_coarse_medium']['max_spinor_L2'],
      'space_refinement':conv['space']['max_spinor_L2']<1e-5,
      'positive_radius':np.min(radius(paths[:,0]))>=A,
      'causal_proper_time':np.max(paths[:,2]-run['times'])<1e-9,
      'all_arrays_finite':all(np.isfinite(v).all() for v in run.values() if isinstance(v,np.ndarray)) and np.isfinite(paths).all()}
    arrays={k:v for k,v in run.items() if isinstance(v,np.ndarray)};arrays.update(direct_paths=paths,quantile_labels=LABELS)
    np.savez_compressed(out.with_suffix('.npz'),**arrays)
    figures=plots(run,paths,args.tag)
    result={'schema':'regular-black-bounce-massive-bohm-v1','source_sha256':sha(__file__),
      'scope':'Massive 3+1 test Dirac on specified Simpson--Visser metric, not a self-consistent renormalized Dirac source',
      'parameters':dict(M=M,a=A,mu=MU,kappa=KAPPA,p0=P0,center=RC,halfwidth=WIDTH,left=LEFT,right=RIGHT,end=END,theta=float(np.pi/2)),
      'grid':dict(dr=run['dr'],dt=run['dt'],points=len(run['r']),profiles=len(run['times'])),
      'convergence':conv,'events':events,'minimum_path_density':pathdensity,
      'max_norm_error':float(max(abs(norms-1))),'max_flux_balance_error':fluxerr,
      'max_guidance_quantile_error':radialerr,'final_region_probabilities':dict(right_exterior=float(norms[-1]-run['probabilities'][-1,2]),interior=float(run['probabilities'][-1,2]-run['probabilities'][-1,0]),left_exterior=float(run['probabilities'][-1,0])),
      'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':bool(all(checks.values())),
      'arrays_file':out.with_suffix('.npz').relative_to(ROOT).as_posix(),'arrays_sha256':sha(out.with_suffix('.npz')),
      'figures':figures,'seconds':time.perf_counter()-begin}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps(result,indent=2),flush=True)
    assert result['all_checks_pass']
if __name__=='__main__':main()
