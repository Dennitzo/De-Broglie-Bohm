"""UV-finite relative four-dimensional mode sum on a minimal OR geometry.

Compares inner Neumann and Dirichlet self-adjoint domains, with the same
outer regulator. This is a complete relative Green function prescription,
not the absolute Hadamard-renormalized RSET or a backreacted solution.
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre
from scipy.sparse import diags, bmat, csr_matrix
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
L=.75; DELTA=11/60

def geometry(R=200.):
    d=DELTA;m=1.
    wc={1:-2*m,4:-32*L*m*m*(15*d+1)/45,5:-2*L*m**3*(60*d+1)/15,6:-16*L*m**4*(75*d+1)/75}
    fc={1:-2*m,2:-2*m*m,3:-8*m**3/3,4:2*m*m*(60*L*d-7*L-90*m*m)/45,
        5:2*m**3*(300*L*d-91*L-720*m*m)/225,6:-4*m**4*(107*L+600*m*m)/225}
    wmax=1+sum(v/R**k for k,v in wc.items())
    f0=sum(v/R**k for k,v in fc.items())
    def inverse(w,y):
        r,g=y;aa=d+(2*d+1/15)*w+(5*d+7/15)*w*w
        D=r*r*w*w*(1-w)+L*(1-w)**2*aa
        return [r**3*w*w/D,(4*r*r*w*(1-w)+L*(1-w)**2*((1/15-4*d)+(13/15-4*d)*w))/(3*D)]
    back=solve_ivp(inverse,(wmax,0),[R,f0+np.log(wmax)/3],rtol=2e-12,atol=2e-14,dense_output=True,method='DOP853')
    assert back.success
    rD,gD=back.y[:,-1]
    def forward(z,y):
        r,g,x=y;w=z**3
        aa=d+(2*d+1/15)*w+(5*d+7/15)*w*w
        D=r*r*w*w*(1-w)+L*(1-w)**2*aa
        rr=3*r**3*z**8/D
        gg=3*z*z*(4*r*r*w*(1-w)+L*(1-w)**2*((1/15-4*d)+(13/15-4*d)*w))/(3*D)
        xx=3*r**3*np.exp(-g/2)*z**7/D
        return [rr,gg,xx]
    sol=solve_ivp(forward,(0,wmax**(1/3)),[rD,gD,0],rtol=2e-12,atol=2e-14,dense_output=True,method='DOP853')
    assert sol.success and abs(sol.y[0,-1]/R-1)<2e-8
    return sol

def coefficients(z,sol):
    r,g,x=sol.sol(z);w=z**3
    aa=DELTA+(2*DELTA+1/15)*w+(5*DELTA+7/15)*w*w
    D=r*r*w*w*(1-w)+L*(1-w)**2*aa
    # Keep x_z*V regular at z=0, where V itself is divergent.
    xz=3*r**3*np.exp(-g/2)*z**7/D
    v0xz=3*np.exp(g/2)/D*(r*z**6*(1-w)+
             L*(1-w)**2/(2*r)*(2*DELTA/3+w/15+(2*DELTA+3/5)*w*w))
    v1xz=3*r*np.exp(g/2)*z**6/D
    return r,g,x,xz,v0xz,v1xz

def compute(sol,probes,ellmax=24,nfreq=40,kmax=8.,rtol=2e-8):
    zmax=sol.t[-1]
    zp=np.array([brentq(lambda z:sol.sol(z)[0]-r,0,zmax,xtol=1e-14) for r in probes])
    nodes,weights=roots_legendre(nfreq)
    # Quadratic frequency coordinate resolves the low cavity eigenfrequencies.
    u=(nodes+1)/2
    ks=kmax*u*u; kw=kmax*u*weights
    ls=np.arange(ellmax+1)
    kappa,ell=np.meshgrid(ks,ls,indexing='ij')
    kappa=kappa.ravel();ell=ell.ravel();n=len(ell)
    lam=ell*(ell+1)
    rB,gB,xB,xz,v0xz,v1xz=coefficients(zmax,sol)
    fB=np.exp(gB)/zmax;wB=zmax**3
    # Common self-adjoint Neumann outer wall. Its radius is varied separately;
    # a finite box must not silently be called the asymptotically flat vacuum.
    zetab=np.zeros(n)
    def fun(z,y):
        r,g,x,xz,v0xz,v1xz=coefficients(z,sol)
        zeta=y[:n]
        return np.r_[xz*(zeta*zeta-kappa*kappa)-v0xz-v1xz*lam,-xz*zeta]
    def jac(z,y):
        xz=coefficients(z,sol)[3]
        return bmat([[diags(2*xz*y[:n]),None],[diags(np.full(n,-xz)),csr_matrix((n,n))]],format='csc')
    states=solve_ivp(fun,(zmax,0),np.r_[zetab,np.zeros(n)],method='Radau',
      rtol=rtol,atol=rtol*.01,jac=jac,dense_output=True)
    assert states.success
    zeta0=states.y[:n,-1];logq0=states.y[n:,-1]
    assert np.min(zeta0)>0
    reports=[];raw={'kappa':ks,'frequency_weights':kw,'ell':ls,'zeta0':zeta0.reshape(nfreq,-1)}
    for i,z in enumerate(zp):
        yy=states.sol(z);zeta=yy[:n];logrel=yy[n:]-logq0
        r,g,x=sol.sol(z);f=np.exp(g)/z;w=z**3
        q2over=np.exp(2*logrel)/zeta0
        measure=np.tile(2*ls+1,nfreq)*np.repeat(kw,ellmax+1)/(4*np.pi*np.pi*r*r)
        rxr=np.sqrt(f*w)/r
        variance=np.sum(measure*q2over)
        et=np.sum(measure*(-kappa*kappa/f)*q2over)
        er=np.sum(measure*(zeta+rxr)**2/f*q2over)
        ea=np.sum(measure*lam/r**2*q2over)
        epsilon=(et+er+ea)/2;pr=(et+er-ea)/2;pt=(et-er)/2
        reports.append({'r':float(r),'x':float(x),'w':float(w),'f':float(f),'gD':float(sol.y[1,0]),'delta_phi2':float(variance),
          'delta_epsilon':float(epsilon),'delta_pr':float(pr),'delta_pt':float(pt),
          'temporal_kinetic':float(et),'radial_kinetic':float(er),'angular_kinetic':float(ea)})
        raw[f'p{i}_relative_kernel']=q2over.reshape(nfreq,-1)
        raw[f'p{i}_zeta']=zeta.reshape(nfreq,-1)
    return {'models':reports,'ellmax':ellmax,'nfreq':nfreq,'kmax':kmax,
       'rtol':rtol,'R_outer':float(sol.y[0,-1]),'solver_steps':len(states.t)},raw

def tail_bounds(row,ellmax,kmax):
    """Joint exponential comparison bound, conditional on exact geometry data.

    g_w >= -40/33 is algebraic for the minimal OR equation. The 1e-6
    downward margin covers observed geometry errors; it is NOT an interval
    arithmetic certification of the floating-point geometry.
    """
    r,x,w,f,gD=[row[k] for k in ['r','x','w','f','gD']]
    lowerf=np.exp(gD-1e-6-40*w/33)*w**(-1/3)
    c1=lowerf/r**2;c0=c1*(1-w)
    theta=np.tanh(x*np.sqrt(c0));a=np.sqrt(2)*x
    q=np.exp(-a*np.sqrt(c1))
    pref=4/(np.sqrt(c0)*theta)/(4*np.pi*np.pi*r*r)
    def sums(n):
        m0=1/(1-q);m1=q/(1-q)**2;m2=q*(1+q)/(1-q)**3;m3=q*(1+4*q+q*q)/(1-q)**4
        s0=q**n*((2*n+1)*m0+2*m1)
        s1=q**n*((2*n**3+3*n*n+n)*m0+(6*n*n+6*n+1)*m1+(6*n+3)*m2+2*m3)
        return s0,s1
    all0,all1=sums(0);tail0,tail1=sums(ellmax+1)
    i0=1/a;i2=2/a**3
    j0=np.exp(-a*kmax)/a
    j2=np.exp(-a*kmax)*(kmax*kmax/a+2*kmax/a**2+2/a**3)
    field=pref*(j0*all0+i0*tail0)
    Fmax=f*np.exp((1-w)/w)
    B1=Fmax/r**2
    B0=Fmax*((1-w)/r**2+L*(1-w)**2*(8*DELTA/3+2/3)/(2*r**4*w*w))
    rxr2=f*w/(r*r)
    cK=3/(2*f);cA=B1/f+1/(2*r*r);cZ=(B0+rxr2)/f
    stress=pref*(cK*(j2*all0+i2*tail0)+cA*(j0*all1+i0*tail1)+cZ*(j0*all0+i0*tail0))
    return {'c0':float(c0),'c1':float(c1),'B0':float(B0),'B1':float(B1),
      'phi2_tail_bound':float(field),'each_stress_component_tail_bound':float(stress),
      'scope':'analytic comparison bound evaluated on numerical geometry with conservative margin; not an interval certificate'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def vector(row):return np.array([row[k] for k in ['delta_phi2','delta_epsilon','delta_pr','delta_pt']])

def compare(a,b):
    return float(max(np.max(abs(vector(x)-vector(y))/np.maximum(abs(vector(y)),1e-14)) for x,y in zip(a['models'],b['models'])))

def main():
    p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--tag',required=True);a=p.parse_args()
    assert a.tag.replace('_','').replace('-','').isalnum()
    out=ROOT/f'work/fullrset_spectral_{a.tag}.json'
    assert not out.exists() and not out.with_suffix('.npz').exists()
    begin=time.perf_counter();sol=geometry()
    probes=np.array([3.,4.,6.])
    if a.pilot:
        rep,raw=compute(sol,probes,ellmax=24,nfreq=32,kmax=16,rtol=1e-7)
        print(json.dumps(rep,indent=2));return
    # The lower resolution runs vary frequency and angular cutoffs separately.
    levels=[]
    for ellmax,nfreq in [(24,80),(64,80),(128,80),(128,128)]:
        rep,_=compute(sol,probes,ellmax=ellmax,nfreq=nfreq,kmax=16,rtol=2e-9)
        levels.append(rep)
        print(json.dumps({'level':(ellmax,nfreq),'seconds':time.perf_counter()-begin}),flush=True)
    # Main source profile and true finite-difference tests of conservation.
    grid=np.geomspace(3,8,45)
    fdsteps=(.003,.0015)
    requested=list(grid)+list(probes)
    for hh in fdsteps:
        for rr in probes:
            requested.extend(rr+hh*np.array([-2,-1,1,2]))
    allprobes=np.unique(requested)
    # The independent second spatial derivative needs tighter dense-output
    # accuracy than values and first derivatives alone.
    mainrep,raw=compute(sol,allprobes,ellmax=128,nfreq=128,kmax=16,rtol=2e-12)
    raw['mode_probe_r']=allprobes
    raw['geometry_z']=np.linspace(0,sol.t[-1],12001)
    for key,vals in zip(['r','g','x'],sol.sol(raw['geometry_z'])):
        raw['geometry_'+key]=vals
    def nearest(rr):return mainrep['models'][int(np.argmin(abs(allprobes-rr)))]
    main={'models':[nearest(rr) for rr in probes]}
    fd=[]
    for hh in fdsteps:
        for rr in probes:
            center=nearest(rr);r=center['r'];w=center['w']
            vals=[nearest(rr+k*hh) for k in [-2,-1,0,1,2]]
            dr=lambda key:(vals[0][key]-8*vals[1][key]+8*vals[3][key]-vals[4][key])/(12*hh)
            d2=lambda key:(-vals[0][key]+16*vals[1][key]-30*vals[2][key]+16*vals[3][key]-vals[4][key])/(12*hh*hh)
            aa=DELTA+(2*DELTA+1/15)*w+(5*DELTA+7/15)*w*w
            bb=-DELTA-6*DELTA*w+(2/5-9*DELTA)*w*w
            wr=(1-w)/r+L*(1-w)**2*aa/(r**3*w*w)
            J=(1-w)/(r*w)+L*(1-w)**2*bb/(3*r**3*w**3)
            ward=dr('delta_pr')+J*(center['delta_epsilon']+center['delta_pr'])/2+2*(center['delta_pr']-center['delta_pt'])/r
            box=w*d2('delta_phi2')+(2*w/r+(wr+w*J)/2)*dr('delta_phi2')
            tr=-center['delta_epsilon']+center['delta_pr']+2*center['delta_pt']
            scale=sum(abs(center[k]) for k in ['delta_epsilon','delta_pr','delta_pt'])
            fd.append({'r':rr,'step':hh,'Ward_scaled':float(abs(ward)/scale),
              'minimal_trace_scaled':float(abs(tr+box/2)/scale)})
    outer=[]
    for radius in [400.,800.]:
        ext=geometry(radius)
        rep,_=compute(ext,probes,ellmax=48,nfreq=128,kmax=16,rtol=2e-9)
        outer.append(rep)
        print(json.dumps({'outer_R':radius,'seconds':time.perf_counter()-begin}),flush=True)
    tails=[tail_bounds(row,128,16) for row in main['models']]
    convergence={'angular_24_to_64':compare(levels[0],levels[1]),
      'angular_64_to_128':compare(levels[1],levels[2]),
      'frequency_80_to_128':compare(levels[2],levels[3]),
      'solver_tolerance':compare(levels[3],main),
      'outer_200_to_400':compare(main,outer[0]),'outer_400_to_800':compare(outer[0],outer[1])}
    checks={'nonzero_domain_dependence':all(row['delta_phi2']>0 for row in main['models']),
      'correct_Wick_sign':all(row['temporal_kinetic']<0 for row in main['models']),
      'angular_convergence':convergence['angular_64_to_128']<2e-6,
      'frequency_convergence':convergence['frequency_80_to_128']<2e-6,
      'solver_convergence':convergence['solver_tolerance']<2e-6,
      'finite_outer_limit_comparison':convergence['outer_400_to_800']<.002,
      'Ward_finite_difference':max(z['Ward_scaled'] for z in fd)<2e-5,
      'minimal_trace_finite_difference':max(z['minimal_trace_scaled'] for z in fd)<2e-5,
      'angular_frequency_remainder':all(b['each_stress_component_tail_bound']<2e-8 for b in tails),
      'raw_arrays_finite':all(np.isfinite(v).all() for v in raw.values())}
    for key in ['r','x','w','f','delta_phi2','delta_epsilon','delta_pr','delta_pt']:
        raw['profile_'+key]=np.array([nearest(rr)[key] for rr in grid])
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'svg.fonttype':'none'})
    fig,ax=plt.subplots(1,3,figsize=(13.8,4.1),layout='constrained')
    rr=raw['profile_r']
    ax[0].loglog(rr,raw['profile_delta_phi2'],color='#0072B2')
    ax[0].set(xlabel='r / √α',ylabel='α Δ⟨Φ²⟩ / (Nℏ)',title='Innere Neumann- minus Dirichlet-Domäne')
    for key,label,col in [('delta_epsilon','Δε','#0072B2'),('delta_pr','Δpᵣ','#D55E00'),('delta_pt','−Δpₜ','#009E73')]:
        vals=raw['profile_'+key]*( -1 if key=='delta_pt' else 1)
        ax[1].loglog(rr,vals,label=label,color=col)
    ax[1].set(xlabel='r / √α',ylabel='α² × Quellendifferenz / (Nℏ)',title='Endliche lokale Stressdifferenz');ax[1].legend()
    ii=np.arange(3)
    ax[2].semilogy(ii,[max(z['Ward_scaled'] for z in fd if z['r']==r) for r in probes],'o-',label='Ward-Rest')
    ax[2].semilogy(ii,[max(z['minimal_trace_scaled'] for z in fd if z['r']==r) for r in probes],'s-',label='Spur-Rest')
    ax[2].set(xticks=ii,xticklabels=['3','4','6'],xlabel='Prüfradius r / √α',ylabel='Skalierter Differenzenrest',title='Unabhängige räumliche Ableitungen');ax[2].legend()
    fig.suptitle('Relative vierdimensionale Modensumme · minimales Feld · feste OR-Geometrie · äußere Wand R=200√α')
    for ext in ('png','svg'):
        dst=ROOT/f'figures/fullrset_{a.tag}_spectral_difference.{ext}';assert not dst.exists()
        fig.savefig(dst,dpi=190,bbox_inches='tight')
    plt.close(fig)
    np.savez_compressed(out.with_suffix('.npz'),**raw)
    report={'scope':'Full UV-finite D/N relative mode prescription; numerical quadrature with analytic joint tails on fixed minimal OR geometry and finite common Neumann outer wall. Not absolute Hadamard RSET or self-consistent backreaction.',
      'normalization':'alpha=1, N*hbar=1 in listed field and stress values; M=1 background',
      'inner_boundary':'Dirichlet psi(0)=0 versus Neumann psi_x(0)=0; difference Neumann minus Dirichlet',
      'outer_boundary':'Common Neumann at R=200; comparisons R=400,800, no claim of exact infinite-radius limit',
      'models':main['models'],'tail_bounds':tails,'convergence':convergence,'fd_checks':fd,
      'checks':checks,'all_checks_pass':all(checks.values()),'refinement_runs':levels,'outer_comparisons':outer,
      'main_run':{k:v for k,v in mainrep.items() if k!='models'},
      'source_sha256':sha(__file__),'arrays_sha256':sha(out.with_suffix('.npz')),
      'background_receipt_sha256':sha(ROOT/'work/fourdim_numerics_20261018_final.json'),
      'elapsed_seconds':time.perf_counter()-begin,'raw_arrays':len(raw)}
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'checks':checks,'convergence':convergence,'tails':tails,'fd':fd,'seconds':report['elapsed_seconds']},indent=2))
    assert report['all_checks_pass']

if __name__=='__main__':main()
