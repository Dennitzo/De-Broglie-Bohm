"""Static dilaton s-wave closure: original Fabbri equations (30)--(32).

This is the explicit normal-ordered t_++=t_--=0
state closure, not a complete spectral renormalized four-dimensional RSET.
alpha=NG*hbar/(12*pi)=1; rho is normalized at the outer boundary and must
not be shifted afterward without transforming the state prescription.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
MASSES=(1.,2.,5.)
R_FACTORS=(50.,100.,200.)
TOLS=(2e-8,2e-10,2e-12)
S_MAX=30.
BOUNDS={'chart_comparison':2e-8,'cross_solver':2e-8,'outer_start_convergence':2e-7,
        'tolerance_convergence':2e-7,'Einstein_scaled':2e-9,'Ward_FD_scaled':2e-5,
        'original_32_scaled':2e-8,'invariant_slack':2e-10,'bounds_slack':2e-9}


def proper_rhs(x,y):
    r,rho,u,v=y
    U=r/3*(1-u*u-2*r*u*v)-r*v*v/3-3*u*v+2*rho*u*u/r
    V=(r+6/r)*u*v-6*rho*u*u/r**2-r*U
    return np.array([u,v,U,V])


def proper_jac(x,y):
    r,p,u,v=y;U=proper_rhs(x,y)[2]
    Ur=(1-u*u-v*v)/3-4*r*u*v/3-2*p*u*u/r**2
    Up=2*u*u/r
    Uu=-2*r*u/3-2*r*r*v/3-3*v+4*p*u/r
    Uv=-2*r*r*u/3-2*r*v/3-3*u
    Vr=(1-6/r**2)*u*v+12*p*u*u/r**3-U-r*Ur
    Vp=-6*u*u/r**2-r*Up
    Vu=(r+6/r)*v-12*p*u/r**2-r*Uu
    Vv=(r+6/r)*u-r*Uv
    return np.array([[0,0,1,0],[0,0,0,1],[Ur,Up,Uu,Uv],[Vr,Vp,Vu,Vv]])


def areal_rhs(r,y):
    rho,w,logh,x=y
    h=np.exp(logh)
    P=r/3*(h-1-w*w-(2*r+9/r)*w+6*rho/r**2)
    W=(r+6/r)*w-6*rho/r**2-(r+w)*P
    return np.array([w,W,-2*P,-np.sqrt(h)])


def redshift_rhs(s,y):
    r,q,logB,x,affine,optical=y
    B=np.exp(logB);rho=-s
    F=B-1-q*q+(2*r+9/r)*q+6*rho*q*q/r**2
    G=-(r+6/r)*q-6*rho*q*q/r**2-r*r*F/3
    return np.array([q,r*F/3+q*G,2*G,-np.exp(logB/2),np.exp(-s+logB/2),np.exp(s+logB/2)])


def outer_data(m,R,order=6):
    if order:
        a={1:-m,2:-m*m,3:-4*m**3/3+13*m/6,4:-2*m**4+19*m*m/12,
           5:-16*m**5/5+43*m**3/60-89*m/40,
           6:-16*m**6/3-3*m**4/2+172*m*m/45}
        b={1:-2*m,2:0,3:-5*m,4:5*m*m/3,5:m**3/2+107*m/4,6:-199*m*m/30}
        rho=sum(a[j]/R**j for j in range(1,order+1))
        f=1+sum(b[j]/R**j for j in range(1,order+1))
        rhop=-sum(j*a[j]/R**(j+1) for j in range(1,order+1))
        return np.array([R,rho,np.sqrt(f),np.sqrt(f)*rhop])
    f=1-2*m/R
    return np.array([R,.5*np.log(f),np.sqrt(f),m/(R*R*np.sqrt(f))])


def integrate_exterior(m,Rfactor=200.,rtol=2e-12):
    R=Rfactor*m
    def throat(x,y):return y[2]
    throat.terminal=True
    sol=solve_ivp(proper_rhs,(0,-2*R),outer_data(m,R),method='Radau',
                  jac=proper_jac,rtol=rtol,atol=rtol*.01,events=throat,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1,sol.message
    return sol


def integrate_inner(throat,smax=S_MAX,rtol=2e-12,method='DOP853'):
    r,rho,u,v=throat
    sol=solve_ivp(redshift_rhs,(-rho,smax),[r,0,-2*np.log(v),0,0,0],
                  method=method,rtol=rtol,atol=rtol*.01,dense_output=True)
    assert sol.success,sol.message
    return sol


def pilot(m,Rfactor=50.,rtol=2e-10):
    R=Rfactor*m
    y0=outer_data(m,R)
    def throat(x,y):return y[2]
    def stop_redshift(x,y):return y[1]+25
    def stop_velocity(x,y):return max(abs(y[2]),abs(y[3]))-1e6
    def stop_radius(x,y):return y[0]-.02
    stop_redshift.terminal=True;stop_velocity.terminal=True;stop_radius.terminal=True
    t=time.perf_counter()
    sol=solve_ivp(proper_rhs,(0,-2*R),y0,method='Radau',rtol=rtol,atol=rtol*.01,
                  jac=proper_jac,events=[throat,stop_redshift,stop_velocity,stop_radius],dense_output=True)
    report={'mass_boundary':m,'outer_radius':R,'success':sol.success,'message':sol.message,
            'nfev':sol.nfev,'elapsed':time.perf_counter()-t,'final_x':float(sol.t[-1]),
            'final_state':sol.y[:,-1].tolist(),'events':[e.tolist() for e in sol.t_events],
            'throat_states':sol.y_events[0].tolist()}
    if len(sol.y_events[0]):
        rb,rhob,ub,vb=sol.y_events[0][0]
        inner=solve_ivp(redshift_rhs,(-rhob,40),[rb,0,-2*np.log(vb),0,0,0],
                        method='Radau',rtol=rtol,atol=rtol*.01,dense_output=True)
        report['redshift_chart']={'success':inner.success,'last_s':float(inner.t[-1]),
            'last_state':inner.y[:,-1].tolist(),
            'samples':{str(s):inner.sol(s).tolist() for s in (5,10,20,40)}}
    return report,sol


def symbolic_audit():
    r,p,u,v,U,V=sp.symbols('r p u v U V',real=True)
    eq30=r*U+V-(r+6/r)*u*v+6*p*u*u/r**2
    eq31=(r+3/r)*U+V-1+u*u+v*v+(r+3/r)*u*v
    solution=sp.solve([eq30,eq31],[U,V])
    Uf=r*(1-u*u-2*r*u*v)/3-r*v*v/3-3*u*v+2*p*u*u/r
    Vf=(r+6/r)*u*v-6*p*u*u/r**2-r*Uf
    assert sp.simplify(solution[U]-Uf)==0 and sp.simplify(solution[V]-Vf)==0
    flow={r:u,p:v,u:Uf,v:Vf}
    derivative=lambda f:sum(sp.diff(f,a)*b for a,b in flow.items())
    W=Uf+u*v
    ratio=sp.factor((derivative(W)+2*v*W)/u)
    assert u not in sp.denom(ratio).free_symbols  # removable throat quotient
    angular_bracket=3*(Vf+v*v)-6*u*v/r+6*p*(-Uf/r+u*u/r**2-u*v/r)+sp.Rational(3,2)*(ratio-W/r)
    angular=-angular_bracket/r**2
    eps=(2*Vf+v*v-3*u*v/r+6*p*u*u/r**2+3*Uf/r)/r**2
    pr=(-v*v-9*u*v/r+6*p*u*u/r**2-3*Uf/r)/r**2
    Einstein=[(1-u*u-2*r*Uf)/r**2,(u*u-1+2*r*u*v)/r**2,Vf+v*v+W/r]
    for a,b in zip(Einstein,[eps,pr,angular]):assert sp.factor(a-b)==0
    ward=sp.factor(derivative(pr)+v*(eps+pr)+2*u/r*(pr-angular))
    assert ward==0
    assert sp.factor(-(W/r+Vf+v*v)-angular_bracket/r**2)==0
    # A finite positive areal radius and regular rho(r) cannot produce h=0.
    w,h=sp.symbols('w h',real=True)
    P=r/3*(h-1-w*w-(2*r+9/r)*w+6*p/r**2)
    assert sp.simplify(Uf.subs({u:1/sp.sqrt(h),v:w/sp.sqrt(h)})*h-P)==0
    # Reconstruct the fixed asymptotic coefficients from the two source ODEs.
    z,m=sp.symbols('z m');aa=sp.symbols('a2:7');bb=sp.symbols('b2:7')
    rho=-m*z+sum(aa[j-2]*z**j for j in range(2,7))
    A=1-2*m*z+sum(bb[j-2]*z**j for j in range(2,7))
    rp=-z*z*sp.diff(rho,z);Ap=-z*z*sp.diff(A,z)
    rpp=2*z**3*sp.diff(rho,z)+z**4*sp.diff(rho,z,2)
    E1=sp.expand(3*z*Ap/2-(1-A-2*A*rp/z)+A*rp**2+9*z*A*rp-6*rho*A*z*z)
    E2=sp.expand(A*rpp+Ap*rp/2-(1/z+6*z)*A*rp+6*rho*A*z*z+Ap/(2*z))
    coefficients={}
    for j in range(2,7):
        equations=[sp.expand(E.subs(coefficients)).coeff(z,j) for E in (E1,E2)]
        coefficients.update(sp.solve(equations,[aa[j-2],bb[j-2]]))
    expected_a=[-m*m,-4*m**3/3+13*m/6,-2*m**4+19*m*m/12,
                -16*m**5/5+43*m**3/60-89*m/40,-16*m**6/3-3*m**4/2+172*m*m/45]
    expected_b=[0,-5*m,5*m*m/3,m**3/2+107*m/4,-199*m*m/30]
    for key,expected in zip((*aa,*bb),(*expected_a,*expected_b)):
        assert sp.simplify(coefficients[key]-expected)==0
    source=sp.lambdify((r,p,u,v),(eps,pr,angular),'numpy',cse=True)
    return {'original_equations_30_31_rederived':True,
            'original_equation_32_identically_satisfied':True,
            'throat_angular_quotient_removable':True,
            'full_Einstein_source_and_Ward_identity':True,
            'areal_h_log_derivative_regular':True,
            'asymptotic_coefficients_derived_to_order_6':True},source


def invariant_arrays(s,Y,entry_s):
    r,q,logb,x,affine,optical=Y
    b=np.exp(logb);rho=-s;z=r*r;k=r*q
    F0=b-1-q*q+(2*r+9/r)*q+6*rho*q*q/r**2
    Ubar=r*F0/3
    G=-(r+6/r)*q-6*rho*q*q/r**2-r*r*F0/3
    Kbar=4*((G+1)**2+2*(q/r)**2+2*(Ubar/r)**2+((b-q*q)/r**2)**2)
    logK=np.log(Kbar)-2*logb
    coefficient=(Ubar-q)/r
    logV=-2*s-logb+np.log(np.maximum(abs(coefficient),np.finfo(float).tiny))
    idx=int(np.argmin(abs(s-entry_s)))
    assert abs(s[idx]-entry_s)<1e-10
    delta=s-entry_s;z0=z[idx];lb0=logb[idx]
    return {'s':s,'r':r,'q':q,'z':z,'k':k,'logb':logb,'proper':x,
            'affine':affine,'optical':optical,'G':G,'Ubar':Ubar,'logK':logK,
            'potential_sign':np.sign(coefficient),'potential_logabs':logV,
            # The areal chart is undefined at the throat q=0.  A masked
            # finite sentinel keeps portable arrays without inventing a value.
            'logh':np.where(q>0,logb-2*np.log(np.maximum(q,np.finfo(float).tiny)),0.),
            'logh_valid':q>0,
            'bound_start_s':np.array([entry_s]),
            'radius_lower_bound':z0+1.5*delta,'radius_upper_bound':z0+2*delta,
            'logb_upper_bound':lb0-(z0/3+2.5)*delta-delta*delta/4,
            'logK_lower_bound':np.log(4.5)-2*np.log(z)-2*logb,
            'proper_tail_log_bound':logb/2-np.log(z/6+1.25),
            'affine_tail_log_bound':-s+logb/2-np.log(z/6+2.25),
            'optical_tail_log_bound':s+logb/2-np.log(z/6+.25)}


def inner_checks(raw,entry_s):
    mask=raw['s']>=entry_s
    slack=BOUNDS['invariant_slack'];d=BOUNDS['bounds_slack']
    return {'strict_entry_in_invariant_region':bool(raw['z'][mask][0]>max(entry_s,16) and .75<raw['k'][mask][0]<1 and raw['logb'][mask][0]<0),
            'invariant_region_samples':bool(np.all(raw['z'][mask]>=np.maximum(raw['s'][mask],16)-slack) and np.all(raw['k'][mask]>=.75-slack) and np.all(raw['k'][mask]<=1+slack) and np.all(raw['logb'][mask]<=slack)),
            'radius_growth_bounds':bool(np.all(raw['z'][mask]>=raw['radius_lower_bound'][mask]-d) and np.all(raw['z'][mask]<=raw['radius_upper_bound'][mask]+d)),
            'Gaussian_logb_bound':bool(np.all(raw['logb'][mask]<=raw['logb_upper_bound'][mask]+d)),
            'curvature_lower_bound':bool(np.all(raw['logK'][mask]>=raw['logK_lower_bound'][mask]-d)),
            'G_barrier_bound':bool(np.all(raw['G'][mask]<=-raw['z'][mask]/6-1.25+d)),
            'logh_finite_on_every_sample':bool(np.all(np.isfinite(raw['logh'][mask])))}


def source_numeric_checks(m,ext,inner,source):
    xb=ext.t[-1];yb=ext.y[:,-1];sb=-yb[1]
    # Fixed overlap by velocity, not by redshift: at larger mass a unit
    # redshift interval already contains exponentially unresolved distances.
    s_overlap=brentq(lambda ss:inner.sol(ss)[2]+12,sb,S_MAX)
    def stop(x,y):return y[3]-np.exp(6.)
    stop.terminal=True
    ins=solve_ivp(proper_rhs,(0,-3),yb,method='Radau',jac=proper_jac,
                  rtol=2e-12,atol=2e-14,events=stop,dense_output=True,max_step=.004)
    assert ins.success and len(ins.t_events[0])==1
    chart=[]
    for s in np.linspace(sb,sb+.95*(s_overlap-sb),31):
        yy=inner.sol(s);rr,q,lb,xx=yy[:4]
        expected=np.array([rr,-s,-q*np.exp(-lb/2),np.exp(-lb/2)])
        chart.append(np.max(abs(ins.sol(xx)-expected)/(1+abs(expected))))
    def point(x):return ins.sol(x) if x<0 else ext.sol(xb+x)
    points=np.r_[np.linspace(ins.t[-1]*.7,ins.t[-1]*.2,5),np.linspace(.1*m,2*m,9)]
    algebra=[]
    for xx in points:
        y=point(xx);r,p,u,v=y;U,V=proper_rhs(xx,y)[2:]
        ein=np.array([(1-u*u-2*r*U)/r**2,(u*u-1+2*r*u*v)/r**2,V+v*v+(U+u*v)/r])
        src=np.asarray(source(r,p,u,v),float)
        algebra.append(np.max(abs(ein-src)/(1+abs(ein)+abs(src))))
    rows=[]
    for h in (.002,.001,.0005):
        worst=0.
        for xx in points:
            y=point(xx);r,p,u,v=y;ep,pr,pt=source(r,p,u,v)
            values=[source(*point(xx+j*h))[1] for j in (-2,-1,1,2)]
            dp=(values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
            terms=[dp,v*(ep+pr),2*u/r*(pr-pt)]
            worst=max(worst,abs(sum(terms))/(1+sum(map(abs,terms))))
        rows.append({'step':h,'Ward_max_scaled_residual':float(worst)})
    return {'proper_vs_redshift_max_scaled':float(max(chart)),
            'Einstein_all_components_scaled':float(max(algebra)),
            'Ward_finite_differences':rows}


def summarize_inner(raw,entry_s):
    i=int(np.argmin(abs(raw['s']-entry_s)))
    return {'entry':{'s':entry_s,'z':float(raw['z'][i]),'k':float(raw['k'][i]),
                'logb':float(raw['logb'][i]),'b':float(np.exp(raw['logb'][i]))},
            'last':{name:float(raw[name][-1]) for name in ('s','r','k','logb','proper','affine','optical','logK','logh','proper_tail_log_bound','affine_tail_log_bound','optical_tail_log_bound')},
            'tail_convention':'E=1 for affine distance; rho fixed by asymptotic state, no later shift',
            'bounds_scope':'Analytic bounds conditional on exact invariant-region entry; exterior entry established numerically, not interval-certified'}


def run_model(m,source):
    ext=integrate_exterior(m);yb=ext.y[:,-1];sb=-yb[1]
    inner=integrate_inner(yb)
    s=np.unique(np.r_[np.linspace(sb,S_MAX,801),10.])
    raw=invariant_arrays(s,inner.sol(s),10.)
    cross=integrate_inner(yb,method='Radau')
    cross_error=float(np.max(abs(cross.sol(s)-inner.sol(s))/(1+abs(inner.sol(s)))))
    fd=source_numeric_checks(m,ext,inner,source)
    radius_rows=[]
    for factor in R_FACTORS:
        other=integrate_exterior(m,factor)
        radius_rows.append({'R_over_m':factor,'throat_r':float(other.y[0,-1]),
            'throat_rho':float(other.y[1,-1]),'throat_v':float(other.y[3,-1]),
            'max_scaled_throat_difference':float(np.max(abs(other.y[:,-1]-yb)/(1+abs(yb))))})
    tolerance_rows=[]
    for tol in TOLS:
        other=integrate_exterior(m,rtol=tol);oi=integrate_inner(other.y[:,-1],rtol=tol)
        ss=np.linspace(10,S_MAX,101)
        error=float(np.max(abs(oi.sol(ss)-inner.sol(ss))/(1+abs(inner.sol(ss)))))
        tolerance_rows.append({'rtol':tol,'max_scaled_inner_difference':error,
            'throat_r':float(other.y[0,-1]),'throat_rho':float(other.y[1,-1])})
    # Exterior plot geometry, no arbitrary inner scattering boundary.
    x=np.linspace(0,8*m,501);Y=ext.sol(ext.t[-1]+x)
    UV=np.stack([proper_rhs(xx,Y[:,j])[2:] for j,xx in enumerate(x)],axis=1)
    U,V=UV;r,p,u,v=Y
    potential=np.exp(2*p)*(U+u*v)/r
    optical=np.zeros_like(x);optical_error=0.
    for j in range(1,len(x)):
        segment,error=quad(lambda xx:np.exp(-ext.sol(ext.t[-1]+xx)[1]),x[j-1],x[j],epsabs=1e-13,epsrel=1e-12)
        optical[j]=optical[j-1]+segment;optical_error+=error
    raw.update({'outer_x_from_throat':x,'outer_r':r,'outer_rho':p,'outer_u':u,'outer_v':v,
                'outer_U':U,'outer_V':V,'outer_potential':potential,'outer_optical':optical,
                'outer_Schwarzschild_potential':(1-2*m/r)*2*m/r**3})
    checks=inner_checks(raw,10.)
    checks.update({'independent_proper_chart':fd['proper_vs_redshift_max_scaled']<BOUNDS['chart_comparison'],
        'DOP853_vs_Radau':cross_error<BOUNDS['cross_solver'],
        'three_outer_start_radii':max(row['max_scaled_throat_difference'] for row in radius_rows)<BOUNDS['outer_start_convergence'],
        'three_tolerance_levels':tolerance_rows[1]['max_scaled_inner_difference']<BOUNDS['tolerance_convergence'],
        'tolerance_refinement':tolerance_rows[1]['max_scaled_inner_difference']<tolerance_rows[0]['max_scaled_inner_difference'],
        'all_Einstein_components':fd['Einstein_all_components_scaled']<BOUNDS['Einstein_scaled'],
        'source_conservation_finite_difference':fd['Ward_finite_differences'][-1]['Ward_max_scaled_residual']<BOUNDS['Ward_FD_scaled'],
        'throat_flare_and_positive_lapse':bool(0<yb[3]<1 and np.exp(2*yb[1])>0),
        'nonzero_radial_backscatter_potential':float(np.max(abs(potential)))>1e-6})
    report={'asymptotic_mass_parameter':m,'outer_radius':200*m,
        'throat':{'r':float(yb[0]),'rho':float(yb[1]),'v':float(yb[3]),'lapse_squared':float(np.exp(2*yb[1]))},
        'inner':summarize_inner(raw,10.),'cross_solver_scaled':cross_error,
        'outer_optical_quadrature_error_estimate':float(optical_error),
        'source_checks':fd,'outer_radius_convergence':radius_rows,
        'tolerance_convergence':tolerance_rows,'checks':{k:bool(v) for k,v in checks.items()}}
    return report,raw


def canonical_model():
    s0=16.;z0=32.;k0=7/8;b0=.5
    y0=[np.sqrt(z0),k0/np.sqrt(z0),np.log(b0),0,0,0]
    sol=solve_ivp(redshift_rhs,(s0,s0+30),y0,method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    alt=solve_ivp(redshift_rhs,(s0,s0+30),y0,method='Radau',rtol=2e-12,atol=2e-14,dense_output=True)
    assert sol.success and alt.success
    s=np.linspace(s0,s0+30,601);raw=invariant_arrays(s,sol.sol(s),s0)
    error=float(np.max(abs(alt.sol(s)-sol.sol(s))/(1+abs(sol.sol(s)))))
    checks=inner_checks(raw,s0);checks['DOP853_vs_Radau']=error<BOUNDS['cross_solver']
    report={'exact_data':{'s0':16,'z0':32,'k0':'7/8','b0':'1/2'},
            'scope':'Exact admissible interior data; no asymptotically flat exterior attachment asserted',
            'inner':summarize_inner(raw,s0),'cross_solver_scaled':error,'checks':checks}
    report['inner']['tail_convention']='Distances start at the exact interior datum s0=16; affine E=1 with rho=-s; no exterior time normalization asserted'
    report['inner']['bounds_scope']='Exact rational initial data lie strictly in the invariant region; mathematical continuation theorem applies without a numerical exterior-entry hypothesis'
    return report,raw


def plot_results(stem,models,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22})
    colors=['#137c8b','#c95632','#6750a4']
    def save(fig,part):
        fig.savefig(ROOT/'figures'/(stem+part+'.png'),dpi=170)
        fig.savefig(ROOT/'figures'/(stem+part+'.svg'),metadata={'Date':None})
        plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.6),layout='constrained')
    for row,col in zip(models,colors):
        m=row['asymptotic_mass_parameter'];a=arrays[f'm{int(m)}']
        near=a['s']<=-row['throat']['rho']+1.5
        x=np.r_[a['proper'][near][::-1],a['outer_x_from_throat'][1:]]/m
        r=np.r_[a['r'][near][::-1],a['outer_r'][1:]]/m
        axs[0].plot(x,r,color=col,label=f'm/√α={m:g}')
        axs[1].plot(a['outer_r']/m,m*m*a['outer_potential'],color=col,label=f'm/√α={m:g}')
    r=np.linspace(2.01,8,301)
    axs[1].plot(r,(1-2/r)*2/r**3,'--',color='#444444',label='Schwarzschild-Vergleich')
    axs[0].axvline(0,color='gray',ls=':');axs[0].set(xlim=(-.7,4),ylim=(1.95,5),xlabel='Eigener Abstand vom Hals ℓ/m',ylabel='Flächenradius r/m',title='Regulärer Hals im dilatonabhängigen Modell')
    axs[1].set(xlim=(2,8),xlabel='Flächenradius r/m',ylabel='Radiales Testmodenpotential m²V₀',title='Rückstreupotential auf der berechneten Geometrie')
    for ax in axs:ax.legend(fontsize=8)
    fig.suptitle('Statischer s-Wellen-Schluss · α=1 · kein Hawking-Nettofluss')
    save(fig,'_geometry')
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.6),layout='constrained')
    for row,col in zip(models,colors):
        m=row['asymptotic_mass_parameter'];a=arrays[f'm{int(m)}'];mask=a['s']>=10
        axs[0].plot(a['s'],a['z'],color=col,label=f'm/√α={m:g}')
        axs[0].fill_between(a['s'][mask],a['radius_lower_bound'][mask],a['radius_upper_bound'][mask],color=col,alpha=.13)
        axs[1].plot(a['s'],a['logK']/np.log(10),color=col,label=f'm/√α={m:g}')
        axs[1].plot(a['s'][mask],a['logK_lower_bound'][mask]/np.log(10),color=col,ls='--',alpha=.6)
    a=arrays['canonical'];axs[0].plot(a['s'],a['z'],':',color='black',label='Exakte innere Anfangsdaten')
    axs[1].plot(a['s'],a['logK']/np.log(10),':',color='black',label='Exakte innere Anfangsdaten')
    axs[0].set(xlabel='Rotverschiebungsvariable s=−ρ',ylabel='z=r²/α',title='Kein endlicher maximaler Flächenradius')
    axs[1].set(xlabel='Rotverschiebungsvariable s=−ρ',ylabel='log₁₀(α²K)',title='Krümmung und analytische Untergrenze')
    for ax in axs:ax.legend(fontsize=8)
    fig.suptitle('Fortsetzung im stabilen Chart und kontrollierter Singularitätsverlauf')
    save(fig,'_interior')
    fig,axs=plt.subplots(1,2,figsize=(11.4,4.6),layout='constrained')
    for row,col in zip(models,colors):
        m=row['asymptotic_mass_parameter'];a=arrays[f'm{int(m)}'];mask=a['s']>=10
        axs[0].plot(a['r'][mask],a['logh'][mask]/np.log(10),color=col,label=f'm/√α={m:g}')
        axs[1].plot(a['s'][mask],a['proper_tail_log_bound'][mask]/np.log(10),color=col,label=f'm/√α={m:g}')
        axs[1].plot(a['s'][mask],a['optical_tail_log_bound'][mask]/np.log(10),color=col,ls='--',alpha=.65)
    axs[0].set(xlabel='Flächenradius r/√α',ylabel='log₁₀ h, h=(dℓ/dr)²',title='Positiver Metrikfaktor bleibt in Logform auflösbar')
    axs[1].set(xlabel='Rotverschiebungsvariable s=−ρ',ylabel='log₁₀ der verbleibenden Längenschranke / √α',title='Endliche Eigen- und optische Distanz')
    for ax in axs:ax.legend(fontsize=8)
    fig.suptitle('Keine Singularitätsdiagnose aus Unterlauf · gestrichelt rechts: optischer Rest')
    save(fig,'_bounds')
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--mass',type=float,default=1.)
    parser.add_argument('--radius-factor',type=float,default=50.)
    parser.add_argument('--tag',default='20261018')
    args=parser.parse_args()
    if args.pilot:
        report,_=pilot(args.mass,args.radius_factor)
        print(json.dumps(report,indent=2),flush=True);return
    if not args.tag.replace('_','').replace('-','').isalnum():parser.error('Invalid tag')
    stem='swave_numerics_'+args.tag;figstem='swave_'+args.tag
    paths=[ROOT/'work'/(stem+ext) for ext in ('.json','.npz','.md')]
    paths += [ROOT/'figures'/(figstem+part+ext) for part in ('_geometry','_interior','_bounds') for ext in ('.png','.svg')]
    if any(p.exists() for p in paths):parser.error('Historical output exists; use a new tag')
    start=time.perf_counter();symbolic,source=symbolic_audit();models=[];arrays={}
    for m in MASSES:
        print(f'Computing asymptotic mass m={m:g}...',flush=True)
        report,raw=run_model(m,source);models.append(report);arrays[f'm{int(m)}']=raw
    canonical,arrays['canonical']=canonical_model()
    mpl=plot_results(figstem,models,arrays)
    flat={f'{name}_{key}':value for name,raw in arrays.items() for key,value in raw.items()}
    assert all(np.all(np.isfinite(a)) for a in flat.values())
    np.savez_compressed(paths[1],**flat)
    checks={f'm{int(row["asymptotic_mass_parameter"])}_{name}':bool(value) for row in models for name,value in row['checks'].items()}
    checks.update({f'canonical_{name}':bool(value) for name,value in canonical['checks'].items()})
    report={'model':{'alpha':1,'asymptotic_masses':MASSES,'outer_radius_over_mass':200,'inner_s_max':S_MAX,
        'source_closure':'Fabbri equations (26)--(32), t_++=t_--=0 with nonzero dilaton angular pressure',
        'scope':'Specified static local state closure; not a complete spectral RSET, no evaporation or thermal emission',
        'no_arbitrary_time_shift':True,'distances':'Exterior-attached models: throat origin, null affine E=1, rho fixed at infinity by series. Canonical interior data: own origin and time convention specified separately.',
        'floating_point':'logb, logh and logK retain deep interior. logh is defined only where logh_valid is true; its throat sentinel 0 is not physical. exp(logb) below float range contributes numerical zero only to exponentially negligible additive terms; positivity follows analytically, not from such underflowed values.'},
        'acceptance_bounds':BOUNDS,'symbolic':symbolic,'models':models,'canonical':canonical,
        'checks':checks,'all_checks_pass':all(checks.values()),'provenance':{
        'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,
        'matplotlib':mpl,'platform':platform.platform(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'primary_source_sha256':hashlib.sha256((ROOT/'inputs/bhback_source_fabbri_2006.html').read_bytes()).hexdigest(),
        'arrays_sha256':hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        'elapsed_seconds':time.perf_counter()-start,'command':f'python work/swave_numerics.py --tag {args.tag}'}}
    paths[0].write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    passed=sum(checks.values());total=len(checks)
    paths[2].write_text(f'''# Dilatonabhängiger statischer s-Wellen-Schluss

Finaler Lauf `{args.tag}`: {passed}/{total} numerische Kontrollen bestanden.
Geometrien mit asymptotischen Massenparametern m/√α=1,2,5, α=1.
Originalgleichungen (30)/(31) werden mit dem nichtverschwindenden Winkeldruck
geschlossen; Gleichung (32) und die volle Erhaltung sind symbolisch geprüft.

Die Rotverschiebungsintegration setzt über den scheinbaren Endpunkt des
Proper-distance-Charts hinaus fort. Alle Außenlösungen liegen bei s=10
numerisch strikt im analytisch invarianten Gebiet. Zusätzlich werden die
exakten inneren Anfangsdaten (s,z,k,b)=(16,32,7/8,1/2) integriert. Für diese
Daten ist der Eintritt keine numerische Annahme. Das bewiesene Gebiet erzwingt
r→∞, K→∞ und endliche Eigen-, affine und optische Randlängen. Ein endlicher
maximaler Flächenradius wird für diese Lösungen nicht reproduziert.

Das radiale Testmodenpotential wird aus derselben Geometrie berechnet. Es
ist keine Berechnung des vollständigen spektralen renormierten Stresses und
keine Hawking-Leistung. Keine unbestimmte innere Streurandbedingung wird
ergänzt. Der Schluss ist statisch und hat keinen Hawking-Nettofluss.

Source SHA256: `{report['provenance']['source_sha256']}`.
NPZ SHA256: `{report['provenance']['arrays_sha256']}`.
Reproduktion mit neuem Tag: `python work/swave_numerics.py --tag NEUER_TAG`.
''',encoding='utf-8')
    print(json.dumps({'passed':passed,'total':total,'failed':[key for key,value in checks.items() if not value],
        'models':models,'canonical':canonical,'provenance':report['provenance']},indent=2),flush=True)
    if not report['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
