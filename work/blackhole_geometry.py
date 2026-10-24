"""Exact Schwarzschild geometry, 2D quantum stress, and a chiral Bohm current.

Three explicitly distinct models: 4D scalar trace on a prescribed metric;
an exact 2D conformal scalar stress tensor; and a 1+1 massless Dirac current.
No claimed self-consistent four-dimensional evaporating black-hole solution.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import scipy
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tensor_and_symbolic_audit():
    t,r,theta,phi=sp.symbols('t r theta phi',real=True)
    m,G=sp.symbols('m G',positive=True)
    coordinates=[t,r,theta,phi]
    f=1-2*m/r
    metric=[-f,1/f,r*r,r*r*sp.sin(theta)**2]
    gamma={}
    def christ(a,b,c):
        if (a,b,c) not in gamma:
            value=0
            if a==c:value+=sp.diff(metric[a],coordinates[b])
            if a==b:value+=sp.diff(metric[a],coordinates[c])
            if b==c:value-=sp.diff(metric[b],coordinates[a])
            gamma[a,b,c]=sp.simplify(value/(2*metric[a]))
        return gamma[a,b,c]
    riemann={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(c+1,4):
                    val=sp.diff(christ(a,b,d),coordinates[c])-sp.diff(christ(a,b,c),coordinates[d])
                    val+=sum(christ(a,e,c)*christ(e,b,d)-christ(a,e,d)*christ(e,b,c) for e in range(4))
                    val=sp.trigsimp(sp.simplify(val))
                    if val!=0:
                        riemann[a,b,c,d]=val
                        riemann[a,b,d,c]=-val
    for b in range(4):
        for d in range(4):
            assert sp.simplify(sum(riemann.get((a,b,a,d),0) for a in range(4)))==0
    K=sp.simplify(sum(metric[a]*v*v/(metric[b]*metric[c]*metric[d])
                      for (a,b,c,d),v in riemann.items()))
    assert sp.simplify(K-48*m*m/r**6)==0
    anomaly=sp.simplify((sp.Rational(1,120)-sp.Rational(1,360))*K/(16*sp.pi**2))
    assert anomaly==m*m/(60*sp.pi**2*r**6)
    leading_R=sp.simplify(-8*sp.pi*G*anomaly)
    assert leading_R==-2*G*m*m/(15*sp.pi*r**6)
    # Exact 2D Polyakov stress in ds^2=-f du dv.  State terms are constants.
    fp=sp.diff(f,r);fpp=sp.diff(fp,r)
    B=(2*f*fpp-fp*fp)/(192*sp.pi)
    uv=f*fpp/(96*sp.pi)
    F=1/(768*sp.pi*m*m)
    outgoing=B+F
    assert sp.simplify(-4*uv/f+fpp/(24*sp.pi))==0
    assert sp.simplify(f*(sp.diff(B,r)-sp.diff(uv,r))+fp*uv)==0
    factored=(r-2*m)**2*(r*r+4*m*r+12*m*m)/(768*sp.pi*m*m*r**4)
    assert sp.simplify(outgoing-factored)==0
    assert sp.limit(outgoing/f**2,r,2*m)==1/(128*sp.pi*m*m)
    assert sp.simplify(B.subs(r,2*m)+F)==0
    w=sp.symbols('w',positive=True)
    # Infall from rest at infinity: du/ds=1/(1-w), dv/ds=1/(1+w).
    uu_ratio=1+2*w*w+3*w**4
    rhoU=(1+w)**2*uu_ratio-8*w**6+w**6*(3*w*w-4)/(1+w)**2
    rhoHH=(1+w)**2*uu_ratio-8*w**6+(1-w)**2*uu_ratio
    rhoB=(3*w**8-4*w**6)*((1-w)**-2+(1+w)**-2)-8*w**6
    assert sp.simplify(rhoU.subs(w,1)-sp.Rational(63,4))==0
    assert sp.simplify(rhoHH.subs(w,1)-16)==0
    assert sp.limit((1-w)**2*rhoB,w,1)==-1
    # PG null coordinate: use r=2m x^2 to avoid sqrt assumptions.
    x=sp.symbols('x',positive=True)
    gg=2*m*x*x-4*m*x+4*m*sp.log(1+x)
    gr=sp.simplify(sp.diff(gg,x)/(4*m*x))
    assert sp.simplify(gr-x/(1+x))==0
    assert sp.limit(gg,x,0)==0
    assert sp.simplify(gg.subs(x,1)-m*(4*sp.log(2)-2))==0
    v=1/x;speed=-1-v
    assert sp.simplify(-(1-v*v)+2*v*speed+speed*speed)==0
    # A C^1 chiral amplitude sqrt(P), compactly supported in Vc-L<V<Vc+L.
    z=sp.symbols('z',real=True)
    prob=sp.Rational(315,256)*(1-z*z)**4
    assert sp.integrate(prob,(z,-1,1))==1
    primitive=sp.Rational(315,256)*(z-4*z**3/3+6*z**5/5-4*z**7/7+z**9/9)
    assert sp.diff(primitive,z)==prob.expand()
    return {'Ricci_tensor_zero_by_Christoffel_contraction':True,
        'nonzero_mixed_Riemann_components':len(riemann),'Kretschmann':str(K),
        'four_dimensional_trace':str(anomaly),'leading_induced_Ricci_scalar':str(leading_R),
        'two_dimensional_trace_and_Ward_identity':True,
        'Unruh_outgoing_double_zero_at_horizon':True,
        'Unruh_regular_ratio_at_horizon':'1/(128*pi*m^2)',
        'freely_falling_horizon_densities_in_units_F':{'Unruh':15.75,'Hartle_Hawking':16.,'Boulware':'diverges negatively'},
        'PG_ingoing_coordinate_derivative_and_null_speed':True,
        'compact_chiral_probability_normalized':True}


def pg_g(r,m=1.):
    x=np.sqrt(np.asarray(r)/(2*m))
    # g=4m(x^2/2-x+log(1+x)); the x^3 series avoids cancellation near r=0.
    small=x<.01
    regular=4*m*(x*x/2-x+np.log1p(x))
    series=4*m*sum((-1)**(n+1)*x**n/n for n in range(3,13))
    return np.where(small,series,regular)


def prob_V(V,center=8.,halfwidth=2.):
    z=(np.asarray(V)-center)/halfwidth
    return np.where(abs(z)<1,315/(256*halfwidth)*(1-z*z)**4,0.)


def cdf_V(V,center=8.,halfwidth=2.):
    z=np.clip((np.asarray(V)-center)/halfwidth,-1.,1.)
    value=.5+315/256*(z-4*z**3/3+6*z**5/5-4*z**7/7+z**9/9)
    return np.clip(value,0.,1.)


def rho2d(r,state):
    w=np.sqrt(2/np.asarray(r));F=1/(768*np.pi)
    common=1+2*w*w+3*w**4
    if state=='Unruh':
        val=(1+w)**2*common-8*w**6+w**6*(3*w*w-4)/(1+w)**2
    elif state=='Hartle_Hawking':
        val=((1+w)**2+(1-w)**2)*common-8*w**6
    elif state=='Boulware':
        val=(3*w**8-4*w**6)*((1-w)**-2+(1+w)**-2)-8*w**6
    else:raise ValueError(state)
    return F*val


def numerical_audit():
    arrays={};checks={};report={}
    r=np.unique(np.r_[np.geomspace(.35,1.99,160),2.,np.geomspace(2.01,25,160)])
    arrays['radius']=r
    arrays['Kretschmann_4d']=48/r**6
    arrays['trace_4d']=1/(60*np.pi**2*r**6)
    arrays['relative_leading_curvature_per_G_over_m2']=1/(30*np.pi*np.sqrt(3)*r**3)
    for state in ('Unruh','Hartle_Hawking'):arrays['rho2d_'+state]=rho2d(r,state)
    exterior=2+np.geomspace(1e-4,18,250)
    arrays['exterior_radius']=exterior
    for state in ('Unruh','Hartle_Hawking','Boulware'):arrays['rho2d_exterior_'+state]=rho2d(exterior,state)
    checks['Unruh_horizon_density']=abs(rho2d(2.,'Unruh')-21/(1024*np.pi))<1e-15
    checks['HH_horizon_density']=abs(rho2d(2.,'Hartle_Hawking')-1/(48*np.pi))<1e-15
    # Compare stable horizon expressions against directly transformed tensors
    # away from their coordinate-singular cancellation region.
    error=0.
    for rr in (.5,1.,1.7,2.3,3.,7.,20.):
        f=1-2/rr;w=np.sqrt(2/rr);B=(3-2*rr)/(48*np.pi*rr**4)
        uv=-f/(24*np.pi*rr**3);F=1/(768*np.pi)
        for state,tu,tv in [('Unruh',F,0.),('Hartle_Hawking',F,F),('Boulware',0.,0.)]:
            direct=(B+tu)/(1-w)**2+2*uv/f+(B+tv)/(1+w)**2
            error=max(error,abs(direct-rho2d(rr,state)))
    checks['direct_tensor_frame_comparison']=error<1e-12
    report['frame_transformation_max_error']=error
    mass_errors=[]
    times=np.linspace(0,10,201);norm=1-cdf_V(times)
    for T in np.linspace(0,10,21):
        lower=max(0.,6-T);upper=max(0.,10-T)
        if upper==0:
            integral=0.
        else:
            rlo=0. if lower==0 else brentq(lambda rr:pg_g(rr)-lower,0,40,xtol=1e-13)
            rhi=brentq(lambda rr:pg_g(rr)-upper,0,40,xtol=1e-13)
            integral=quad(lambda rr:prob_V(T+pg_g(rr))/(1+np.sqrt(2/rr)),rlo,rhi,
                          epsabs=2e-12,epsrel=2e-12)[0]
        mass_errors.append(abs(integral-(1-cdf_V(T))))
    checks['probability_mass_and_outflow']=max(mass_errors)<1e-10
    arrays['probability_time']=times;arrays['surviving_probability']=norm
    arrays['singularity_outflow']=prob_V(times)
    # Fixed labels are quantiles of the exact initial normalized distribution.
    labels=np.array([brentq(lambda V:cdf_V(V)-p,6,10,xtol=1e-13) for p in np.linspace(.05,.95,13)])
    rows=[];curve=np.full((len(labels),len(times)),np.nan)
    for i,V in enumerate(labels):
        r0=brentq(lambda rr:pg_g(rr)-V,0,40,xtol=1e-13)
        rstop=.02;stop_time=V-float(pg_g(rstop))
        evaluation=np.linspace(0,stop_time,151)
        def rhs(T,y):return [-1-np.sqrt(2/max(y[0],1e-30))]
        sol=solve_ivp(rhs,(0,stop_time),[r0],t_eval=evaluation,method='DOP853',rtol=2e-12,atol=2e-13)
        assert sol.success,sol.message
        residual=np.max(abs(evaluation+pg_g(sol.y[0])-V))
        rerr=abs(sol.y[0,-1]-rstop)
        assert residual<2e-9 and rerr<1e-8
        valid=times<V
        curve[i,valid]=[brentq(lambda rr:pg_g(rr)-(V-T),0,40,xtol=1e-12) for T in times[valid]]
        rows.append({'label_V':float(V),'initial_r':r0,'horizon_PG_time':float(V-pg_g(2.)),
                     'singularity_PG_time':float(V),'integrated_to_r':rstop,
                     'null_coordinate_max_drift':float(residual),'stop_radius_error':float(rerr)})
    arrays['trajectory_time']=times;arrays['trajectory_radius']=curve;arrays['trajectory_labels']=labels
    checks['horizon_crossing_and_finite_endpoint']=all(0<x['horizon_PG_time']<x['singularity_PG_time'] for x in rows)
    checks['independent_trajectory_ODE']=all(x['null_coordinate_max_drift']<2e-9 for x in rows)
    # A two-dimensional grid checks true finite-difference current divergence.
    residuals=[]
    for step in (1e-3,5e-4,2.5e-4):
        worst=0.
        for T in (0.,4.,6.,8.):
            for rr in (.1,.5,2.,4.,8.,12.):
                def density(tt,x):return prob_V(tt+pg_g(x))/(1+np.sqrt(2/x))
                def flux(tt,x):return -prob_V(tt+pg_g(x))
                div=(density(T+step,rr)-density(T-step,rr))/(2*step)
                div+=(flux(T,rr+step)-flux(T,rr-step))/(2*step)
                worst=max(worst,abs(div))
        residuals.append({'step':step,'current_divergence':float(worst)})
    checks['current_continuity_finite_difference']=residuals[-1]['current_divergence']<2e-7
    checks['FD_convergence']=residuals[-1]['current_divergence']<residuals[0]['current_divergence']/10
    report.update({'m':1.,'probability_mass_max_error':float(max(mass_errors)),
                   'current_finite_differences':residuals,'trajectories':rows,
                   'checks':{k:bool(v) for k,v in checks.items()}})
    return report,arrays


def plot(root,stem,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22})
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.5),layout='constrained')
    r=arrays['radius'];ext=arrays['exterior_radius']
    for key,label in [('Unruh','Unruh'),('Hartle_Hawking','Hartle–Hawking'),('Boulware','Boulware')]:
        axes[0].plot(ext-2,arrays['rho2d_exterior_'+key],label=label)
    axes[0].set_yscale('symlog',linthresh=1e-3)
    axes[0].set(xscale='log',xlabel='(r − 2m)/m',ylabel='m² ρ beim freien Fall',title='2D-Stress: Zustandswahl am Horizont')
    axes[0].legend()
    axes[1].loglog(r,arrays['Kretschmann_4d'],label='m⁴ K (4D)')
    axes[1].loglog(r,arrays['trace_4d'],label='m⁴ ⟨Tᵃₐ⟩ (4D)')
    axes[1].axvline(2,color='gray',ls=':',label='Horizont r=2m')
    axes[1].set(xlabel='r/m',ylabel='Dimensionslose Invarianten',title='Schwarzschild: Horizont endlich, Zentrum singulär')
    axes[1].legend()
    for suffix in ('.png','.svg'):fig.savefig(root/'figures'/(stem+'_stress'+suffix),dpi=165,metadata={'Date':None} if suffix=='.svg' else None)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.7),layout='constrained')
    for i,rr in enumerate(arrays['trajectory_radius']):axes[0].plot(arrays['trajectory_time'],rr,lw=1.2)
    axes[0].axhline(2,color='black',ls='--',label='Horizont')
    axes[0].axhline(0,color='black',lw=.8)
    axes[0].set(xlabel='Horizontreguläre PG-Zeit T/m',ylabel='r/m',title='Bohmscher einlaufender 1+1-Diracstrom',ylim=(0,22))
    axes[0].legend()
    axes[1].plot(arrays['probability_time'],arrays['surviving_probability'],label='Verbleibende Wahrscheinlichkeit N(T)')
    axes[1].plot(arrays['probability_time'],arrays['singularity_outflow'],label='Ausfluss am Singularitätsrand')
    axes[1].set(xlabel='T/m',ylabel='N bzw. m × Ausflussrate',title='Erhaltung mit offenem Rand bei r=0')
    axes[1].legend()
    for suffix in ('.png','.svg'):fig.savefig(root/'figures'/(stem+'_current'+suffix),dpi=165,metadata={'Date':None} if suffix=='.svg' else None)
    plt.close(fig)
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='20261018')
    parser.add_argument('--pilot',action='store_true')
    args=parser.parse_args()
    if not args.tag.replace('_','').replace('-','').isalnum():parser.error('Invalid tag')
    stem='blackhole_geometry_'+args.tag
    outputs=[ROOT/'work'/(stem+e) for e in ('.json','.npz')]
    outputs += [ROOT/'figures'/(stem+s+e) for s in ('_stress','_current') for e in ('.png','.svg')]
    if not args.pilot and any(p.exists() for p in outputs):parser.error('Use a fresh tag')
    symbolic=tensor_and_symbolic_audit()
    numerical,arrays=numerical_audit()
    assert all(numerical['checks'].values()),numerical['checks']
    result={'schema':'bohm-blackhole-geometry-v1','utc':datetime.now(timezone.utc).isoformat(),
            'symbolic':symbolic,'numerical':numerical,'all_checks_pass':True,
            'scope':'Prescribed Schwarzschild geometry: 4D anomaly; separate exact 2D conformal scalar stress and 1+1 Dirac current. No self-consistent interior solution.'}
    if args.pilot:
        print(json.dumps(result,indent=2));return
    version=plot(ROOT,stem,arrays)
    np.savez_compressed(outputs[1],**arrays)
    result['provenance']={'source_sha256':sha(Path(__file__)),'arrays_sha256':sha(outputs[1]),
                         'python':platform.python_version(),'numpy':np.__version__,
                         'scipy':scipy.__version__,'sympy':sp.__version__,'matplotlib':version}
    outputs[0].write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'success':True,'checks':len(numerical['checks']),'files':[p.relative_to(ROOT).as_posix() for p in outputs]},indent=2))


if __name__=='__main__':main()
