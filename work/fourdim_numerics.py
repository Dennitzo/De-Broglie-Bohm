"""Self-consistent vacuum OR-AHS source, Arrechea et al. (2023), eq. (9).

The angular stress is completed by the exact static Ward identity.  This is
the analytic order-reduced approximation, not an exact four-dimensional RSET.
The printed evolution equations (15),(16) have errors; the source (9), rather
than those transcriptions, defines this implementation.  alpha=GN*hbar/(12*pi)
is held at one, matching the preceding Polyakov and dilaton comparisons.

Run: python work/fourdim_numerics.py --tag 20261018
Existing tagged outputs are never overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 1.0
L = 3 * ALPHA / 4
XI_C = 11 / 60
MASSES = (1.0, 2.0, 5.0)
COUPLINGS = (0.0, 1 / 6)
RTOLS = (2e-8, 2e-10, 2e-12)
OUTER_LEVELS = (50.0, 100.0, 200.0)
BOUNDS = {
    'cross_solver_scaled': 2e-8,
    'log_chart_scaled': 2e-8,
    'areal_chart_scaled': 2e-8,
    'outer_start_scaled': 2e-7,
    'tolerance_scaled': 2e-7,
    'Einstein_scaled': 2e-9,
    'Ward_FD_scaled': 2e-5,
    'endpoint_scaling_relative': 2e-4,
    'distance_scaling_relative': 2e-4,
    'curvature_log_slope': 2e-4,
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scaled(a, b):
    return float(np.max(np.abs(np.asarray(a) - b) / (1 + np.abs(b))))


def polynomials(w, delta):
    aa = delta + (2 * delta + 1 / 15) * w + (5 * delta + 7 / 15) * w * w
    bb = -delta - 6 * delta * w + (2 / 5 - 9 * delta) * w * w
    return aa, bb


def exterior_data(mass, xi, outer_multiple):
    """Independent asymptotic solution through r^-6, f(infinity)=1."""
    m, d, R = mass, XI_C - xi, outer_multiple * mass
    wc = {1: -2*m, 4: -32*L*m*m*(15*d+1)/45,
          5: -2*L*m**3*(60*d+1)/15, 6: -16*L*m**4*(75*d+1)/75}
    fc = {1: -2*m, 2: -2*m*m, 3: -8*m**3/3,
          4: 2*m*m*(60*L*d-7*L-90*m*m)/45,
          5: 2*m**3*(300*L*d-91*L-720*m*m)/225,
          6: -4*m**4*(107*L+600*m*m)/225}
    w = 1 + sum(v/R**n for n,v in wc.items())
    logf = sum(v/R**n for n,v in fc.items())
    return R, w, logf


def inverse_rhs(w, state, delta):
    """Regular at w=0 for delta>0; state=(r,g), g=log(f)+log(w)/3."""
    r, g = state
    aa, _ = polynomials(w, delta)
    dd = r*r*w*w*(1-w) + L*(1-w)**2*aa
    rw = r**3*w*w/dd
    gw = (4*r*r*w*(1-w) + L*(1-w)**2 *
          ((1/15-4*delta)+(13/15-4*delta)*w))/(3*dd)
    return np.array([rw, gw])


def inverse_solution(mass, xi, outer_multiple=200.0, rtol=2e-12, method='DOP853'):
    R,w0,F0 = exterior_data(mass,xi,outer_multiple)
    sol = solve_ivp(lambda w,y: inverse_rhs(w,y,XI_C-xi), (w0,0),
                    (R,F0+np.log(w0)/3), method=method, rtol=rtol,
                    atol=rtol*0.01, dense_output=True)
    if not sol.success or sol.t[-1] != 0:
        raise RuntimeError('Regular inverse chart did not reach its specified endpoint: '+sol.message)
    return sol


def build_symbolics():
    r,w,d,l=sp.symbols('r w d l', positive=True)
    aa=d+(2*d+sp.Rational(1,15))*w+(5*d+sp.Rational(7,15))*w**2
    bb=-d-6*d*w+(sp.Rational(2,5)-9*d)*w**2
    wr=(1-w)/r+l*(1-w)**2*aa/(r**3*w**2)
    p=(1-w)/(r*w)+l*(1-w)**2*bb/(3*r**3*w**3)
    deriv=lambda z: sp.diff(z,r)+sp.diff(z,w)*wr
    eps=-l*(1-w)**2*aa/(r**4*w**2)
    pr=l*(1-w)**2*bb/(3*r**4*w**2)
    pt=pr+r*deriv(pr)/2+r*p*(eps+pr)/4
    a=w*(deriv(p)/2+p*p/4)+wr*p/4
    b=w*p/(2*r);c=wr/(2*r);dd=(1-w)/r**2
    tensor=(-(dd-2*c),2*b-dd,a+b+c)
    xi=sp.Rational(11,60)-d
    source9_tt=8*l*(1-w)**2/r**4*((1/w**2+6/w+33)/480-(xi-sp.Rational(1,6))*(1/w**2+2/w+5)/8)
    source9_rr=8*l*(1-w)**2/r**4*(-(1/w**2+6/w-15)/1440+(xi-sp.Rational(1,6))*(1/w+3)**2/24)
    checks={
        'source9_tt_normalization':sp.factor(source9_tt+eps)==0,
        'source9_rr_normalization':sp.factor(source9_rr-pr)==0,
        'Einstein_tt_from_source9':sp.factor(tensor[0]+eps)==0,
        'Einstein_rr_from_source9':sp.factor(tensor[1]-pr)==0,
        'Einstein_angular_Ward_completion':sp.factor(tensor[2]-pt)==0,
        'Ward_identity_off_shell_along_flow':sp.factor(deriv(pr)+p*(pr+eps)/2+2*(pr-pt)/r)==0,
    }
    # Check coefficients of the exterior series independently from the stored ones.
    x,m=sp.symbols('x m');ww=1-2*m*x
    wc={4:-32*l*m*m*(15*d+1)/45,5:-2*l*m**3*(60*d+1)/15,
        6:-16*l*m**4*(75*d+1)/75}
    ww+=sum(v*x**n for n,v in wc.items())
    ff=-2*m*x-2*m*m*x*x-8*m**3*x**3/3+2*m*m*(60*l*d-7*l-90*m*m)*x**4/45
    ff+=2*m**3*(300*l*d-91*l-720*m*m)*x**5/225-4*m**4*(107*l+600*m*m)*x**6/225
    aa2=aa.subs(w,ww);bb2=bb.subs(w,ww)
    wrseries=sp.Poly(sp.expand(-x*sp.diff(ww,x)*ww**2-(1-ww)*ww**2-l*(1-ww)**2*aa2*x*x),x)
    pfseries=sp.Poly(sp.expand(-3*x*sp.diff(ff,x)*ww**3-3*ww**2*(1-ww)-l*(1-ww)**2*bb2*x*x),x)
    checks['outer_series_w_through_order6']=all(sp.simplify(wrseries.nth(n))==0 for n in range(7))
    checks['outer_series_logf_through_order6']=all(sp.simplify(pfseries.nth(n))==0 for n in range(7))
    funcs=sp.lambdify((r,w,d,l),(wr,p,deriv(p),eps,pr,pt,a,b,c,dd), 'numpy',cse=True)
    return checks, funcs


def diagnostic(r,w,delta,funcs):
    wr,p,prr,eps,pr,pt,a,b,c,dd=[np.asarray(v) for v in funcs(r,w,delta,L)]
    K=4*(a*a+2*b*b+2*c*c+dd*dd)
    Ricci=-2*a-4*b-4*c+2*dd
    return {'w_r':wr,'logf_r':p,'logf_rr':prr,'epsilon_8piG':eps,
            'radial_pressure_8piG':pr,'angular_pressure_8piG':pt,
            'Ricci':Ricci,'Kretschmann':K,'logK':np.log(K),
            'G_energy':dd-2*c,'G_radial':2*b-dd,'G_angular':a+b+c}


def integrals_to_endpoint(w,sol,delta):
    """No cancellation of r-rD: quadrature from the regular endpoint chart."""
    def f(z,which):
        if z == 0:return 0.0
        r,g=sol.sol(z);rw=inverse_rhs(z,(r,g),delta)[0]
        if which==0:return rw
        if which==1:return rw/np.sqrt(z)
        if which==2:return rw*np.exp(g/2)*z**(-2/3)
        return rw*np.exp(-g/2)*z**(-1/3)
    vals=[]
    for which in range(4):
        # Rescaling keeps relative quadrature accurate even for tiny gaps.
        val=quad(lambda t:f(w*t,which),0,1,epsabs=1e-200,epsrel=2e-11,limit=100)[0]*w
        vals.append(val)
    return np.array(vals)


def one_model(mass,xi,funcs):
    delta=XI_C-xi
    sol=inverse_solution(mass,xi)
    w0=sol.t[0];rD,gD=sol.y[:,-1]
    w=np.geomspace(w0,1e-10,601)
    r,g=sol.sol(w);logf=g-np.log(w)/3
    ar={'w':w,'r':r,'g':g,'logf':logf,'logh':-np.log(w)}
    ar.update(diagnostic(r,w,delta,funcs))
    # Independently chosen integration coordinate s=log(h), state=(r,log(f)).
    def log_rhs(s,y):
        ww=np.exp(-s);rr,F=y;aa,bb=polynomials(ww,delta)
        dd=rr*rr*ww*ww*(1-ww)+L*(1-ww)**2*aa
        N=rr*rr*ww*ww*(1-ww)+L*(1-ww)**2*bb/3
        return [-rr**3*ww**3/dd,-N/dd]
    ss=-np.log(w)
    other=solve_ivp(log_rhs,(ss[0],ss[-1]),(r[0],logf[0]),rtol=2e-12,
                    atol=2e-14,method='Radau',t_eval=ss)
    if not other.success:raise RuntimeError(other.message)
    log_error=scaled(other.y,np.vstack((r,logf)))
    ar['independent_logchart_r']=other.y[0];ar['independent_logchart_logf']=other.y[1]
    rad=inverse_solution(mass,xi,method='Radau')
    cross=scaled(rad.sol(w),np.vstack((r,g)))
    # Original areal-radius equations integrated up to a non-stiff overlap boundary.
    w_stop=.08
    rend=float(sol.sol(w_stop)[0])
    def radial_rhs(rr,y):
        ww,F=y;aa,bb=polynomials(ww,delta)
        return [(1-ww)/rr+L*(1-ww)**2*aa/(rr**3*ww*ww),
                (1-ww)/(rr*ww)+L*(1-ww)**2*bb/(3*rr**3*ww**3)]
    direct=solve_ivp(radial_rhs,(r[0],rend),(w0,logf[0]),rtol=2e-12,
                     atol=2e-14,method='DOP853',dense_output=True)
    if not direct.success:raise RuntimeError(direct.message)
    ww=np.geomspace(w0,w_stop,101);rr,gg=sol.sol(ww)
    areal_error=scaled(direct.sol(rr),np.vstack((ww,gg-np.log(ww)/3)))
    # Fixed outer-start and tolerance comparisons; endpoint regularity permits direct comparisons.
    radius=[];tolerances=[]
    for Rm in OUTER_LEVELS:
        z=inverse_solution(mass,xi,outer_multiple=Rm)
        radius.append({'R_over_M':Rm,'endpoint_r':float(z.y[0,-1]),
                       'endpoint_g':float(z.y[1,-1]),'scaled_endpoint_difference':scaled(z.y[:,-1],(rD,gD))})
    for tol in RTOLS:
        z=inverse_solution(mass,xi,rtol=tol)
        tolerances.append({'rtol':tol,'scaled_state_difference':scaled(z.sol(w),np.vstack((r,g))),
                           'endpoint_r':float(z.y[0,-1])})
    # Full Einstein components from independent tetrad curvature and source9+Ward.
    ei=max(scaled(ar['G_energy'],ar['epsilon_8piG']),
           scaled(ar['G_radial'],ar['radial_pressure_8piG']),
           scaled(ar['G_angular'],ar['angular_pressure_8piG']))
    # Five-point differentiation of p_r in the smooth inverse chart, compared to Ward.
    wf=np.geomspace(.6,2e-3,61);rf,gf=sol.sol(wf)
    diag=diagnostic(rf,wf,delta,funcs)
    ward=[]
    for step in (.004,.002,.001):
        derivatives=0.0
        for off,coef in ((-2,1),(-1,-8),(1,8),(2,-1)):
            wp=wf*(1+off*step);rp=sol.sol(wp)[0]
            derivatives+=coef*diagnostic(rp,wp,delta,funcs)['radial_pressure_8piG']
        dprdw=derivatives/(12*step*wf)
        dd=dprdw*diag['w_r']
        b=diag['logf_r']*(diag['epsilon_8piG']+diag['radial_pressure_8piG'])/2
        c=2*(diag['radial_pressure_8piG']-diag['angular_pressure_8piG'])/rf
        err=float(np.max(np.abs(dd+b+c)/(1+np.abs(dd)+np.abs(b)+np.abs(c))))
        ward.append({'relative_w_step':step,'max_scaled_residual':err})
    # Endpoint asymptotics using integrated gaps, not subtraction of nearly equal radii.
    tailw=np.geomspace(1e-3,1e-9,61)
    integ=np.array([integrals_to_endpoint(z,sol,delta) for z in tailw])
    tr,tg=sol.sol(tailw);tf=np.exp(tg)*tailw**(-1/3)
    td=diagnostic(tr,tailw,delta,funcs)
    C=rD**3/(3*L*delta);a=L*delta/rD**3
    fD=np.exp(gD)*C**(1/9)
    chi=8*a**(1/3)/3**(11/3)
    ratios={
        'radius_gap_ratio':integ[:,0]/(C*tailw**3),
        'lapse_asymptotic_ratio':tf*integ[:,0]**(1/9)/fD,
        'Ricci_asymptotic_ratio':(-td['Ricci'])*integ[:,0]**(5/3)/chi,
        'proper_tail_ratio':integ[:,1]/((2/5)/a*tailw**(5/2)),
        'affine_tail_ratio':integ[:,2]/((3/7)*np.exp(gD/2)/a*tailw**(7/3)),
        'optical_tail_ratio':integ[:,3]/((3/8)*np.exp(-gD/2)/a*tailw**(8/3)),
    }
    ar.update({'tail_w':tailw,'tail_r':tr,'tail_g':tg,'tail_logf':tg-np.log(tailw)/3,
               'tail_radius_gap':integ[:,0],'tail_proper':integ[:,1],
               'tail_affine_E1':integ[:,2],'tail_optical':integ[:,3],
               'tail_logK':td['logK'],'tail_Ricci':td['Ricci'],**ratios})
    slope=float(np.polyfit(np.log(integ[-15:,0]),td['logK'][-15:],1)[0])
    # Turning point of f is unique here; p changes sign while f remains positive.
    def p_of_w(z):
        rr=sol.sol(z)[0]
        return float(diagnostic(rr,z,delta,funcs)['logf_r'])
    wm=brentq(p_of_w,1e-7,.95,xtol=1e-14)
    rm,gm=sol.sol(wm);fm=float(np.exp(gm)*wm**(-1/3))
    checks={
        'regular_inverse_chart_reaches_w_zero':bool(sol.success and sol.t[-1]==0),
        'singular_radius_outside_2M':bool(rD>2*mass),
        'positive_lapse_minimum':bool(fm>0),
        'DOP853_vs_Radau':cross<BOUNDS['cross_solver_scaled'],
        'independent_log_chart':log_error<BOUNDS['log_chart_scaled'],
        'original_areal_chart':areal_error<BOUNDS['areal_chart_scaled'],
        'three_outer_start_levels':max(v['scaled_endpoint_difference'] for v in radius)<BOUNDS['outer_start_scaled'],
        'three_tolerance_levels':tolerances[1]['scaled_state_difference']<BOUNDS['tolerance_scaled'],
        'tolerance_refinement':tolerances[1]['scaled_state_difference']<tolerances[0]['scaled_state_difference'],
        'all_Einstein_components':ei<BOUNDS['Einstein_scaled'],
        'Ward_finite_difference':ward[-1]['max_scaled_residual']<BOUNDS['Ward_FD_scaled'],
        'Ward_refinement':ward[-1]['max_scaled_residual']<ward[0]['max_scaled_residual'],
        'endpoint_radius_lapse_Ricci_scaling':max(abs(ratios[k][-1]-1) for k in list(ratios)[:3])<BOUNDS['endpoint_scaling_relative'],
        'finite_distance_tail_scalings':max(abs(ratios[k][-1]-1) for k in list(ratios)[3:])<BOUNDS['distance_scaling_relative'],
        'Kretschmann_power_minus_10_over_3':abs(slope+10/3)<BOUNDS['curvature_log_slope'],
    }
    checks={k:bool(v) for k,v in checks.items()}
    rep={'mass':mass,'xi':xi,'delta':delta,'endpoint':{'rD':float(rD),'gD':float(gD),
          'gap_coefficient_C':float(C),'fD':float(fD),'Ricci_chi_corrected':float(chi)},
          'lapse_minimum':{'w':float(wm),'r':float(rm),'f':fm},
          'cross_solver_scaled':cross,'log_chart_scaled':log_error,'areal_chart_scaled':areal_error,
          'Einstein_scaled':ei,'Ward_finite_differences':ward,
          'outer_start_convergence':radius,'tolerance_convergence':tolerances,
          'final_endpoint_ratios':{k:float(v[-1]) for k,v in ratios.items()},
          'Kretschmann_log_gap_slope':slope,'checks':checks}
    return rep,ar


def plot_results(models,arrays,tag,comparison):
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,
                          'svg.fonttype':'none','figure.dpi':120})
    colors=['#0072B2','#D55E00','#009E73']
    def save(fig,name):
        fig.savefig(ROOT/'figures'/f'fourdim_{tag}_{name}.png',dpi=190,bbox_inches='tight')
        fig.savefig(ROOT/'figures'/f'fourdim_{tag}_{name}.svg',bbox_inches='tight')
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11.2,4.2),layout='constrained')
    for model in models:
        a=arrays[model['prefix']];m=model['mass'];color=colors[MASSES.index(m)]
        ls='-' if model['xi']==0 else '--'; mask=a['r']<8*m
        label=f'M/√α={m:g}, ξ='+('0' if model['xi']==0 else '1/6')
        ax[0].plot(a['r'][mask]/m,np.exp(a['logf'][mask]),color=color,ls=ls,label=label)
        ax[1].plot(a['r'][mask]/m,a['w'][mask],color=color,ls=ls)
    rr=np.linspace(2,8,300)
    for aa in ax:aa.plot(rr,1-2/rr,':',color='#555555',label='Schwarzschild');aa.set_xlabel('r / M');aa.set_xlim(1.9,8)
    ax[0].set(yscale='log',ylabel='f = −gₜₜ',ylim=(1e-4,20),title='Positives Lapse-Minimum; anschließend Divergenz')
    ax[1].set(ylabel='w = 1 / h',title='Endlicher Singularitätsradius, kein regulärer Hals')
    ax[0].legend(fontsize=8,ncol=2)
    fig.suptitle('Vierdimensionale OR-AHS-Näherung · α=1 · durchgezogen ξ=0, gestrichelt ξ=1/6')
    save(fig,'geometry')
    fig,ax=plt.subplots(1,2,figsize=(11.2,4.2),layout='constrained')
    for model in models:
        a=arrays[model['prefix']];m=model['mass'];color=colors[MASSES.index(m)];ls='-' if model['xi']==0 else '--'
        gap=a['tail_radius_gap'];h=1/a['tail_w']
        ax[0].plot(np.log10(gap),np.log10(h),color=color,ls=ls,label=f'M={m:g}, ξ='+('0' if model['xi']==0 else '1/6'))
        ax[1].plot(np.log10(gap),a['tail_logK']/np.log(10),color=color,ls=ls)
    ax[0].set(xlabel='log₁₀[(r − rD) / √α]',ylabel='log₁₀ h',title='h ∝ (r − rD)⁻¹ᐟ³')
    ax[1].set(xlabel='log₁₀[(r − rD) / √α]',ylabel='log₁₀(α² K)',title='K ∝ (r − rD)⁻¹⁰ᐟ³')
    ax[0].legend(fontsize=8,ncol=2)
    fig.suptitle('Regulär integrierter Endpunkt · Radiusabstand ohne Auslöschung quadraturberechnet')
    save(fig,'endpoint')
    fig,ax=plt.subplots(1,2,figsize=(11.2,4.2),layout='constrained')
    a=arrays['m1_minimal'];mask=(a['w']<.5)&(a['w']>1e-5)
    for key,label,col in [('epsilon_8piG','−ε','#0072B2'),('radial_pressure_8piG','−pᵣ','#D55E00'),('angular_pressure_8piG','pₜ','#009E73')]:
        sign=-1 if key!='angular_pressure_8piG' else 1
        ax[0].plot(-np.log10(a['w'][mask]),np.log10(sign*a[key][mask]),label=label,color=col)
    ax[0].set(xlabel='log₁₀ h',ylabel='log₁₀[8πG α × Betrag]',title='Anisotrope Quelle, M/√α=1, ξ=0');ax[0].legend()
    ms=np.array(MASSES)
    for xi,label,col in [(0.0,'OR: Singularitätsradius, ξ=0','#0072B2'),(1/6,'OR: Singularitätsradius, ξ=1/6','#D55E00')]:
        radii=np.array([z['endpoint']['rD'] for z in models if z['xi']==xi])
        ax[1].plot(ms,radii-2*ms,'o-',label=label,color=col)
    if comparison:
        ax[1].plot(ms,[v['throat_r']-2*v['mass'] for v in comparison],'^--',color='#777777',label='Dilaton-s-Welle: regulärer Hals, ξ=0')
    ax[1].set(xlabel='ADM-Masse M / √α',ylabel='(rₛₚₑz − 2M) / √α',title='Modellvergleich bei gleichem M und α');ax[1].legend(fontsize=8)
    fig.suptitle('Erwartungswertquelle und Modellvergleich · statischer Boulware-Zustand')
    save(fig,'source_comparison')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tag',required=True);parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    sys.stdout.reconfigure(encoding='utf-8')
    stem=ROOT/'work'/f'fourdim_numerics_{args.tag}'
    paths=[stem.with_suffix(z) for z in ('.json','.npz','.md')]
    if any(p.exists() for p in paths):raise SystemExit('Tagged outputs exist; choose a fresh tag.')
    start=time.perf_counter();symbolic,funcs=build_symbolics();models=[];arrays={}
    xis=(0.0,) if args.pilot else COUPLINGS
    masses=(1.0,) if args.pilot else MASSES
    for xi in xis:
        for mass in masses:
            key=f'm{mass:g}_'+('minimal' if xi==0 else 'conformal')
            rep,ar=one_model(mass,xi,funcs);rep['prefix']=key;models.append(rep);arrays[key]=ar
            print(json.dumps({'model':key,'endpoint':rep['endpoint'],'checks_passed':sum(rep['checks'].values()),'checks':len(rep['checks'])}),flush=True)
    checks={**symbolic,**{m['prefix']+'_'+k:bool(v) for m in models for k,v in m['checks'].items()}}
    raw={key+'_'+name:v for key,a in arrays.items() for name,v in a.items()}
    checks['all_raw_arrays_finite']=all(np.isfinite(v).all() for v in raw.values())
    compare_path=ROOT/'work/swave_numerics_20261018_final.json'
    comparison=[]
    if not args.pilot:
        cj=json.loads(compare_path.read_text(encoding='utf-8'))
        comparison=[{'mass':m['asymptotic_mass_parameter'],'throat_r':m['throat']['r']} for m in cj['models']]
        plot_results(models,arrays,args.tag,comparison)
    np.savez_compressed(paths[1],**raw)
    report={'model':{'alpha':ALPHA,'L_GNhbar_over_16pi':L,'xi_c':XI_C,'masses':list(masses),'couplings':list(xis),
            'state':'Static Boulware, kappa=0, f(infinity)=1','source':'Vacuum OR-AHS source from Arrechea 2023 equation (9), angular component completed through exact Ward identity',
            'scope':'Self-consistent analytic approximation, not the exact full four-dimensional spectral RSET; no evaporation, no experimental fit',
            'normalization':'8*pi*G times stresses are archived; alpha=G*N*hbar/(12*pi)=1 matches earlier Polyakov and dilaton models',
            'endpoint':'w=1/h reaches zero regularly; g=log(f)+log(w)/3 finite; positive coupling gap delta=xi_c-xi',
            'distance_convention':'Distances are positive outward from the singular endpoint; affine normalization Killing E=1 with f(infinity)=1'},
            'acceptance_bounds':BOUNDS,'symbolic':symbolic,'models':models,'checks':checks,'all_checks_pass':all(checks.values()),
            'comparison':comparison,'provenance':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform(),
                'source_sha256':sha(__file__),'primary_source_sha256':sha(ROOT/'inputs/bhback_source_arrechea_2023.html'),
                'arrays_sha256':sha(paths[1]),'comparison_receipt_sha256':sha(compare_path) if comparison else None,
                'elapsed_seconds':time.perf_counter()-start,'command':f'python work/fourdim_numerics.py --tag {args.tag}'+(' --pilot' if args.pilot else '')}}
    paths[0].write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    txt='# Vierdimensionale OR-AHS-Numerik\n\n'
    txt+=f'Gültiger Lauf: `{args.tag}`. Prüfungen: {sum(checks.values())}/{len(checks)}; Roharrays: {len(raw)}.\n\n'
    txt+='Quelle ist die vierdimensionale analytische OR-Näherung im statischen Boulware-Zustand. Sie ist keine exakte spektrale RSET. Die Winkelspannung wird aus der vollständigen Ward-Identität rekonstruiert. α=GNℏ/(12π)=1; L=GNℏ/(16π)=3/4.\n\n'
    txt+='Die reguläre inverse Chart erreicht w=1/h=0 ohne Solver-Abbruch. Die asymptotische Fortsetzung besitzt einen endlichen rD, f∝(r−rD)^−1/9, h∝(r−rD)^−1/3 und K∝(r−rD)^−10/3. Eigen-, affine (E=1) und optische Reststrecken sind endlich. Abstände zum Endpunkt werden durch positive Quadratur bestimmt und nicht durch Subtraktion fast gleicher Radien.\n\n'
    txt+='Die gedruckten Gleichungen (15) und (16) enthalten Fehler; maßgeblich ist Quelle (9). Der Mittelkoeffizient in h′ ist 2(ξ−13/60), die Konstante im f′/f-Polynom 9(ξ−5/36). Auch der Ricci-Vorfaktor aus (27) wird unabhängig nachgerechnet: χ=8(Lδ/rD³)^(1/3)/3^(11/3).\n\n'
    txt+='| M/√α | ξ | rD/√α | minimales f |\n|---:|---:|---:|---:|\n'
    for m in models:txt+=f"| {m['mass']:g} | {m['xi']:.9g} | {m['endpoint']['rD']:.12g} | {m['lapse_minimum']['f']:.12g} |\n"
    txt+='\nDas Vergleichsbild kennzeichnet den regulären s-Wellen-Hals und den OR-Singularitätsradius als verschiedene geometrische Objekte. Gleiche α- und ADM-Massennormierung bedeutet keine Gleichsetzung der Quellenmodelle.\n\n'
    txt+='```json\n'+json.dumps({'checks':checks,'provenance':report['provenance']},indent=2,ensure_ascii=False)+'\n```\n'
    paths[2].write_text(txt,encoding='utf-8')
    print(json.dumps({'all_checks_pass':report['all_checks_pass'],'passed':sum(checks.values()),'checks':len(checks),'seconds':report['provenance']['elapsed_seconds']}),flush=True)
    if not report['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
