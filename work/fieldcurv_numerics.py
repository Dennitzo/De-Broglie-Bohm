"""Smooth Dirac excitation: stress, geometric work and Killing energy.

Two fixed geometries; no absolute vacuum RSET or nonlinear backreaction.
Uses the frozen, separately validated DG and characteristic spatial cores.
All new outputs use fresh tags; previous research is never overwritten.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad,simpson,cumulative_simpson
from scipy.interpolate import CubicSpline
from nonsing_independent import BounceDG
from nonsing_bohm import maps

ROOT=Path(__file__).resolve().parents[1]
END=12.;DT_OUT=.025
TIMES=np.linspace(0,END,round(END/DT_OUT)+1)
NORMALIZATION=3*quad(lambda x:np.exp(-2/(1-x*x)) if abs(x)<1 else 0.,-1,1,epsabs=2e-14,epsrel=2e-14)[0]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def initial(r):
    z=(np.asarray(r)-12)/3;env=np.zeros_like(z);inside=abs(z)<1
    env[inside]=np.exp(-1/(1-z[inside]**2))/np.sqrt(NORMALIZATION)
    p=-2.;E=np.hypot(p,1);pol=np.sqrt([(E+p)/(2*E),(E-p)/(2*E)])
    return env[...,None]*np.exp(1j*p*(np.asarray(r)-12))[...,None]*pol

def geometry(r,a):
    R=np.sqrt(r*r+a*a);v=np.sqrt(2/R);rp=r/R;vp=-v*r/(2*R**2)
    ca=(2*r*r-a*a)/R**5;cb=-r*r/R**5
    cc=(2*a*a-r*r-a*a*R)/R**5;cd=(2*r*r+a*a*R)/R**5
    K=4*(ca*ca+2*cb*cb+2*cc*cc+cd*cd)
    return R,v,rp,vp,K

def stress(p,pr,r,a):
    """Last axis is spinor; sphere-integrated projections have no 4piR²."""
    R,v,rp,vp,K=geometry(r,a);z=np.conj(p[...,0])*p[...,1]
    q=np.sum(abs(p)**2,axis=-1)
    Pi=np.imag(np.sum(np.conj(p)*pr,axis=-1))
    P=np.imag(np.conj(p[...,0])*pr[...,0]-np.conj(p[...,1])*pr[...,1])
    s1=2*z.real;s2=2*z.imag;W=s1-s2/R
    e=P+W;pt=-s2/(2*R)
    ek=e-v*Pi;fk=(1+v*v)*Pi-v*(2*P+W)
    fe=Pi-v*e;work=vp*P+2*v*rp/R*pt
    number=(1-v)*abs(p[...,0])**2-(1+v)*abs(p[...,1])**2
    factor=4*np.pi*R*R
    # Positive frame norm, explicitly not a Lorentz scalar.
    stressnorm=np.sqrt(e*e+P*P+2*pt*pt+2*Pi*Pi)/factor
    return dict(q=q,number=number,e=e,P=P,Pi=Pi,pt=pt,ek=ek,fk=fk,fe=fe,work=work,
        rho=e/factor,pressure=P/factor,pt_phys=pt/factor,momentum=Pi/factor,
        trace=(-e+P+2*pt)/factor,scalar=s1/factor,
        nec_minus=(e+P+2*Pi)/factor,
        source_scale=8*np.pi*stressnorm/np.sqrt(K))

def mesh(a,h):
    pieces=[(.5,3,h/2),(3,30,h)] if a==0 else [(-30,-2.3,h),(-2.3,-1.2,h/4),(-1.2,30,h)]
    out=[]
    for lo,hi,step in pieces:
        x=np.linspace(lo,hi,int(np.ceil((hi-lo)/step))+1)
        out.extend(x if not out else x[1:])
    return np.array(out)

def dg_sampling(solver,points,derivative=False):
    cell=np.clip(np.searchsorted(solver.edges,points,side='right')-1,0,solver.cells-1)
    local=2*(points-(solver.edges[cell]+solver.edges[cell+1])/2)/solver.h[cell]
    weights=np.empty((len(points),solver.n))
    for i in range(solver.n):
        poly=np.poly1d([1.])
        for j in range(solver.n):
            if i!=j:poly*=np.poly1d([1.,-solver.xi[j]])/(solver.xi[i]-solver.xi[j])
        weights[:,i]=(np.polyder(poly)(local)*2/solver.h[cell] if derivative else poly(local))
    def evaluate(state):return np.sum(state[..., :,cell,:]*weights,axis=-1).swapaxes(-2,-1)
    return evaluate

def char_run(a,dr=.005,dt=.000625):
    left=.5 if a==0 else -30.;right=30.
    r=np.linspace(left,right,round((right-left)/dr)+1);mp=maps(r,dt,a=a)
    R=np.sqrt(r*r+a*a);om=np.sqrt(1+1/R**2)
    c=np.cos(om*dt/2);s=np.sinc(om*dt/(2*np.pi))*dt/2
    upper=s*(1/R-1j);lower=s*(-1/R-1j)
    p=initial(r);states=[p.copy()];stride=round(DT_OUT/dt)
    for step in range(1,round(END/dt)+1):
        mid=np.column_stack([c*p[:,0]+upper*p[:,1],c*p[:,1]+lower*p[:,0]])
        streamed=np.empty_like(mid)
        for k,(ix,w,jac,valid) in enumerate(mp):
            streamed[:,k]=np.sum(mid[ix,k]*w,axis=1)*jac
            streamed[~valid,k]=0
        p=np.column_stack([c*streamed[:,0]+upper*streamed[:,1],c*streamed[:,1]+lower*streamed[:,0]])
        if step%stride==0:states.append(p.copy())
    assert len(states)==len(TIMES)
    return r,np.array(states)

def dg_run(a,h):
    solver=BounceDG(1,edges=mesh(a,h),a=a)
    p=initial(solver.r).transpose(2,0,1)
    state=solver.evolve(p,TIMES)
    return solver,state

def quadrature(solver):
    x,w=leggauss(9)
    r=(solver.edges[:-1,None]+solver.edges[1:,None])/2+solver.h[:,None]*x/2
    return r.ravel(),(solver.h[:,None]*w/2).ravel()

def diagnose(solver,state,a):
    r,w=quadrature(solver);sample=dg_sampling(solver,r);derivative=dg_sampling(solver,r,True)
    p=sample(state);pr=derivative(state);d=stress(p,pr,r,a)
    rb=solver.edges[[0,-1]]
    pb=dg_sampling(solver,rb)(state);prb=dg_sampling(solver,rb,True)(state)
    bd=stress(pb,prb,rb,a)
    result={key:np.sum(d[key]*w,axis=-1) for key in ['q','e','ek','work']}
    for key in ['number','fk','fe']:
        result[key+'_boundary']=bd[key][:,1]-bd[key][:,0]
        result[key+'_cumulative']=cumulative_simpson(result[key+'_boundary'],x=TIMES,initial=0)
    result['work_cumulative']=cumulative_simpson(result['work'],x=TIMES,initial=0)
    out,jump,penalty=solver.rates(state)
    result['jump_cumulative']=cumulative_simpson(jump+penalty,x=TIMES,initial=0)
    result['norm_balance']=result['q']-result['q'][0]+result['number_cumulative']
    result['killing_balance']=result['ek']-result['ek'][0]+result['fk_cumulative']
    result['euler_balance']=result['e']-result['e'][0]+result['fe_cumulative']-result['work_cumulative']
    result['source_scale_max']=d['source_scale'].max(axis=1)
    result['trace_identity_max']=np.max(abs(d['trace']+d['scalar']))
    result['null_source_min']=np.min(d['nec_minus'],axis=1)
    return result,(r,w,p,pr,d)

def compare(solver,state,r_other,psi_other,a):
    r,w=quadrature(solver);p=dg_sampling(solver,r)(state);pr=dg_sampling(solver,r,True)(state)
    d=stress(p,pr,r,a)
    diffs={k:[] for k in ['spinor_L2','derivative_L2','energy_density_L1','killing_density_L1','work_L1']}
    for i,t in enumerate(TIMES):
        interp=CubicSpline(r_other,psi_other[i],axis=0);po=interp(r);pro=interp(r,1)
        do=stress(po,pro,r,a)
        diffs['spinor_L2'].append(np.sqrt(np.sum(abs(po-p[i])**2*w[:,None])))
        diffs['derivative_L2'].append(np.sqrt(np.sum(abs(pro-pr[i])**2*w[:,None])))
        for label,key in [('energy_density_L1','e'),('killing_density_L1','ek'),('work_L1','work')]:
            diffs[label].append(np.sum(abs(do[key]-d[key][i])*w))
    return {k:float(max(v)) for k,v in diffs.items()}

def dg_compare(coarse,cs,fine,fs,a):
    r,w=quadrature(fine)
    p=dg_sampling(fine,r)(fs);pr=dg_sampling(fine,r,True)(fs)
    cp=dg_sampling(coarse,r)(cs);cr=dg_sampling(coarse,r,True)(cs)
    d=stress(p,pr,r,a);dc=stress(cp,cr,r,a)
    return dict(spinor_L2=float(np.max(np.sqrt(np.sum(abs(p-cp)**2*w[None,:,None],axis=(1,2))))),
        derivative_L2=float(np.max(np.sqrt(np.sum(abs(pr-cr)**2*w[None,:,None],axis=(1,2))))),
        energy_density_L1=float(np.max(np.sum(abs(d['e']-dc['e'])*w,axis=1))),
        killing_density_L1=float(np.max(np.sum(abs(d['ek']-dc['ek'])*w,axis=1))),
        work_L1=float(np.max(np.sum(abs(d['work']-dc['work'])*w,axis=1))))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tag',required=True)
    ap.add_argument('--model',type=int,choices=[0,1],required=True)
    ap.add_argument('--h',type=float,default=.05);ap.add_argument('--dr',type=float,default=.005)
    ap.add_argument('--dt',type=float,default=.000625);args=ap.parse_args()
    stem=ROOT/f'work/fieldcurv_{args.tag}_a{args.model}';assert not stem.with_suffix('.json').exists() and not stem.with_suffix('.npz').exists()
    started=time.perf_counter();a=args.model
    coarse,cs=dg_run(a,args.h*2);fine,fs=dg_run(a,args.h)
    diag,(r,w,p,pr,d)=diagnose(fine,fs,a)
    refine=dg_compare(coarse,cs,fine,fs,a)
    print(json.dumps(dict(phase='DG_done',model=a,seconds=time.perf_counter()-started,refinement=refine)),flush=True)
    rc,pc=char_run(a,args.dr,args.dt);comparison=compare(fine,fs,rc,pc,a)
    balances={key:float(np.max(abs(diag[key]))) for key in ['norm_balance','killing_balance','euler_balance']}
    # Second independent temporal quadrature resolution, on unchanged DG fields.
    coarse_time=TIMES[::2]
    time_quad={}
    for key in ['number','fk','fe']:
        val=cumulative_simpson(diag[key+'_boundary'][::2],x=coarse_time,initial=0)
        time_quad[key]=float(np.max(abs(val-diag[key+'_cumulative'][::2])))
    val=cumulative_simpson(diag['work'][::2],x=coarse_time,initial=0)
    time_quad['work']=float(np.max(abs(val-diag['work_cumulative'][::2])))
    checks=dict(initial_norm=abs(diag['q'][0]-1)<1e-9,
        norm_balance=balances['norm_balance']<1e-6,
        killing_balance=balances['killing_balance']<2e-5,
        euler_balance=balances['euler_balance']<2e-5,
        trace_identity=diag['trace_identity_max']<1e-12,
        DG_field_refinement=refine['spinor_L2']<2e-5,
        DG_derivative_refinement=refine['derivative_L2']<5e-4,
        DG_stress_refinement=refine['energy_density_L1']<1e-4,
        independent_field=comparison['spinor_L2']<2e-5,
        independent_derivative=comparison['derivative_L2']<5e-4,
        independent_stress=comparison['energy_density_L1']<1e-4,
        independent_work=comparison['work_L1']<1e-4,
        temporal_quadrature=max(time_quad.values())<1e-4,
        negligible_jump_loss=max(diag['jump_cumulative'])<1e-7,
        finite_fields=np.isfinite(fs).all() and np.isfinite(pc).all())
    # Plot profiles at .01 spacing use actual DG polynomials, including derivatives.
    plot_r=np.linspace(.5 if a==0 else -20,16,1551 if a==0 else 3601)
    plot_p=dg_sampling(fine,plot_r)(fs)
    plot_pr=dg_sampling(fine,plot_r,True)(fs)
    plot_d=stress(plot_p,plot_pr,plot_r,a)
    arrays=dict(times=TIMES,edges=fine.edges,nodes=fine.r,states=fs,
        coarse_edges=coarse.edges,coarse_nodes=coarse.r,coarse_states=cs,
        char_r=rc,char_psi=pc,plot_r=plot_r,plot_psi=plot_p,plot_derivative=plot_pr)
    arrays.update({'diag_'+k:v for k,v in diag.items()})
    arrays.update({'plot_'+k:v for k,v in plot_d.items()})
    assert all(np.isfinite(v).all() for v in arrays.values())
    np.savez_compressed(stem.with_suffix('.npz'),**arrays)
    record=dict(schema='field-curvature-smooth-excitation-v1',source_sha256=sha(__file__),
        input_hashes={f'work/{n}':sha(ROOT/'work'/n) for n in ['nonsing_independent.py','nonsing_bohm.py']},
        a=a,M=1,mu=1,kappa=1,p0=-2,center=12,halfwidth=3,normalization=NORMALIZATION,
        end=END,output_step=DT_OUT,interval=[float(fine.left),float(fine.right)],
        grid=dict(coarse_cells=coarse.cells,fine_cells=fine.cells,degree=5,hmin=float(min(fine.h)),hmax=float(max(fine.h)),char_dr=args.dr,char_dt=args.dt),
        initial=dict(norm=float(diag['q'][0]),euler_energy=float(diag['e'][0]),killing_energy=float(diag['ek'][0])),
        final=dict(norm=float(diag['q'][-1]),euler_energy=float(diag['e'][-1]),killing_energy=float(diag['ek'][-1]),integrated_work=float(diag['work_cumulative'][-1])),
        balances=balances,temporal_quadrature=time_quad,refinement=refine,comparison=comparison,
        max_source_scale=float(max(diag['source_scale_max'])),min_null_excitation=float(min(diag['null_source_min'])),
        jump_loss=float(max(diag['jump_cumulative'])),checks={k:bool(v) for k,v in checks.items()},all_checks_pass=bool(all(checks.values())),
        arrays_file=stem.with_suffix('.npz').relative_to(ROOT).as_posix(),arrays_sha256=sha(stem.with_suffix('.npz')),seconds=time.perf_counter()-started)
    stem.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps(record,indent=2),flush=True)
    assert record['all_checks_pass'],{k:v for k,v in checks.items() if not v}

if __name__=='__main__':main()
