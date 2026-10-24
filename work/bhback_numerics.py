"""Self-consistent static spherical Einstein--Polyakov Boulware model.

alpha=GN*hbar/(12*pi)=1.  The conserved 2D Polyakov stress is lifted by
1/(4*pi*r**2), with zero angular pressure.  This is a specified effective
model, not the full 4D scalar RSET and not an evaporating Unruh state.
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
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
RADII=(2.,4.,10.)
S_MAX=18.
TOLS=(2e-8,2e-10,2e-12)
OUTER_FACTORS=(50.,150.,450.,1350.)
# These acceptance bounds are set before production runs.
BOUNDS={
    'constraint_scaled':2e-8,
    'inner_cross_method_scaled':2e-8,
    'proper_chart_scaled':2e-8,
    'outer_ADM_last_change':2e-6,
    'outer_lapse_last_change':2e-7,
    'Einstein_FD_scaled':2e-7,
    'Ward_FD_scaled':2e-7,
    'tolerance_refinement_scaled':2e-8,
    'asymptotic_logD_slope_error':2e-3,
    'analytic_bound_slack':2e-11,
    'independent_mass_chart':2e-8,
}


def proper_rhs(x,y):
    r,rho,u,v=y
    den=r*r-1
    return np.array([u,v,v*((r*r+1)*u+r*v)/den,
                     -r*v*(2*u+r*v)/den])


def q_from_logD(r,logD):
    D=np.exp(logD)
    return (1-D)/(r+np.sqrt(r*r-1+D))


def inner_rhs(s,y):
    r,logD,proper,affine=y
    q=q_from_logD(r,logD)
    return np.array([q,-2*r*(r-2*q)/(r*r-1),
                     -np.exp(logD/2),np.exp(-s+logD/2)])


def inner_solution(r0,rtol=2e-12,method='DOP853'):
    sol=solve_ivp(inner_rhs,(0,S_MAX),[r0,0.,0.,0.],method=method,
                  rtol=rtol,atol=rtol*.01,dense_output=True)
    assert sol.success,sol.message
    return sol


def proper_solution(r0,rtol=2e-12,rstop=None,sstop=None,method='DOP853',max_step=np.inf):
    assert (rstop is None)!=(sstop is None)
    if rstop is not None:
        def event(x,y):return y[0]-rstop
        span=(0.,2*rstop+10*r0)
    else:
        def event(x,y):return y[1]+sstop
        span=(0.,-10*r0)
    event.terminal=True
    sol=solve_ivp(proper_rhs,span,[r0,0.,0.,1.],method=method,
                  rtol=rtol,atol=rtol*.01,dense_output=True,events=event,max_step=max_step)
    assert sol.success and len(sol.t_events[0])==1,sol.message
    return sol


def exterior_mass_chart(r0,outer):
    rstart=2*r0
    xstart=brentq(lambda x:outer.sol(x)[0]-rstart,0,outer.t[-1])
    rr,rho,u,v=outer.sol(xstart)
    mstart=rr*(1-u*u)/2
    def rhs(z,y):
        r=rstart*np.exp(z);M,rho=y
        u=np.sqrt(1-2*M/r)
        v=(2*M/r)/(r*u+np.sqrt(r*r*u*u+2*M/r))
        vp=-r*v*(2*u+r*v)/(r*r-1)
        return [r*(2*vp+v*v)/2,r*v/u]
    sol=solve_ivp(rhs,(0,np.log(OUTER_FACTORS[-1]/2)),[mstart,rho],
                  method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
    assert sol.success
    return sol


def asymptotic_completion(radius,mass,rho):
    # Directly derived large-r series, alpha=1:
    # M(r)=m+m/r^2+m^2/(6r^3)+(m^3/4+m)/r^4+O(r^-5).
    m=mass
    for _ in range(4):
        f=m*(1+radius**-2)+m*m/(6*radius**3)+(m**3/4+m)/radius**4-mass
        df=1+radius**-2+m/(3*radius**3)+(3*m*m/4+1)/radius**4
        m-=f/df
    rho_inf=rho-.5*np.log1p(-2*m/radius)+m/(3*radius**3)+11*m*m/(12*radius**4)
    return m,rho_inf


def invariants(r,u,v):
    den=r*r-1
    up=v*((r*r+1)*u+r*v)/den
    vp=-r*v*(2*u+r*v)/den
    # Stable forms of A=rho''+rho'^2 and the remaining orthonormal
    # sectional curvatures.  Constraint is separately checked.
    A=(u*u-1)/den
    B=u*v/r
    C=up/r
    E=(1-u*u)/(r*r)
    K=4*(A*A+2*B*B+2*C*C+E*E)
    R=-2*A-4*B-4*C+2*E
    energy=(2*vp+v*v)/(r*r)  # 8*pi*G times 4D energy density
    pressure=-v*v/(r*r)     # 8*pi*G times radial pressure
    return {'rpp':up,'rhopp':vp,'K':K,'Ricci':R,
            'source_energy':energy,'source_radial_pressure':pressure,
            'mass':r*(1-u*u)/2}


def symbolic_audit():
    r,u,v=sp.symbols('r u v',real=True)
    up=v*((r*r+1)*u+r*v)/(r*r-1)
    vp=-r*v*(2*u+r*v)/(r*r-1)
    constraint=u*u+2*r*u*v+v*v-1
    derivative=sp.factor(sp.diff(constraint,r)*u+sp.diff(constraint,u)*up+sp.diff(constraint,v)*vp)
    assert derivative==0
    # All Einstein components, reduced modulo the independent constraint.
    equations=[1-u*u-2*r*up-(2*vp+v*v),
               u*u-1+2*r*u*v+v*v,
               r*(vp+v*v)+u*v+up]
    quotients=[]
    for eq in equations:
        quotients.append(str(sp.factor(eq/constraint)) if eq!=0 else '0')
        assert sp.simplify(eq.subs(u*u,1-2*r*u*v-v*v))==0 or sp.factor(eq/constraint) in (-1,0,1)
    z=sp.symbols('z',real=True)
    rho=sp.Function('rho')(z)
    ep=2*sp.diff(rho,z,2)+sp.diff(rho,z)**2
    pp=-sp.diff(rho,z)**2
    assert sp.simplify(sp.diff(pp,z)+sp.diff(rho,z)*(ep+pp))==0
    assert sp.simplify(-ep+pp+2*(sp.diff(rho,z,2)+sp.diff(rho,z)**2))==0
    # In arbitrary spherical geometry the 1/r^2 factor cancels the two
    # extra radial divergence terms when the angular pressure is zero.
    rr=sp.Function('r')(z);ee=sp.Function('ee')(z);pr=sp.Function('pr')(z)
    ward4=sp.diff(pr/rr**2,z)+sp.diff(rho,z)*(ee+pr)/rr**2+2*sp.diff(rr,z)/rr*(pr/rr**2)
    assert sp.simplify(ward4-(sp.diff(pr,z)+sp.diff(rho,z)*(ee+pr))/rr**2)==0
    # Derive (rather than fit) the asymptotic completion from the mass ODE.
    zz,m,aa,bb,cc=sp.symbols('zz m aa bb cc')
    MM=m+aa*zz**2+bb*zz**3+cc*zz**4
    uu=sp.sqrt(1-2*MM*zz)
    vv=2*MM*zz**2/(uu+sp.sqrt(uu**2+2*MM*zz**3))
    us=sp.series(uu,zz,0,6).removeO();vs=sp.series(vv,zz,0,6).removeO()
    Mr=sp.series(-(2*zz*us*vs+vs**2)/(1-zz**2)+vs**2/2,zz,0,6).removeO()
    residual=sp.expand(-zz**2*sp.diff(MM,zz)-Mr)
    coefficients=sp.solve([residual.coeff(zz,k) for k in range(3,6)],[aa,bb,cc])
    assert coefficients=={aa:m,bb:m*m/6,cc:m**3/4+m}
    rhoprime=sp.series(vs/us,zz,0,6).removeO().subs(coefficients)
    tail=sp.expand(sp.integrate(rhoprime/zz**2,(zz,0,zz)))
    expected=m*zz+m*m*zz**2+(4*m**3/3+m/3)*zz**3+(2*m**4+11*m*m/12)*zz**4
    assert sp.simplify(tail-expected)==0
    return {'constraint_derivative_exactly_zero':True,
            'all_three_Einstein_components':True,
            'Einstein_residual_constraint_multipliers':quotients,
            'two_dimensional_trace_and_Ward_identity':True,
            'four_dimensional_uplift_Ward_identity':True,
            'asymptotic_mass_and_lapse_coefficients_derived':True}


def finite_difference_checks(r0,outer,inside):
    # Differentiate interpolated geometry, not the RHS, for Einstein tests.
    def values(x):return inside.sol(x) if x<0 else outer.sol(x)
    points=np.r_[np.linspace(inside.t[-1]*.65,inside.t[-1]*.12,6),
                 np.linspace(.15*r0,2*r0,8)]
    rows=[]
    for h in (1e-2,5e-3,2.5e-3):
        errs=[];wards=[]
        for x in points:
            ym2,ym1,y0,yp1,yp2=[values(x+j*h) for j in (-2,-1,0,1,2)]
            first=(ym2-8*ym1+8*yp1-yp2)/(12*h)
            second=(-yp2+16*yp1-30*y0+16*ym1-ym2)/(12*h*h)
            r,rho=y0[:2];u,v=first[:2];upp,vpp=second[:2]
            einstein=np.array([(1-u*u-2*r*upp)/(r*r),
                              (u*u-1+2*r*u*v)/(r*r),
                              vpp+v*v+(u*v+upp)/r])
            source=np.array([(2*vpp+v*v)/(r*r),-v*v/(r*r),0.])
            errs.append(float(np.max(abs(einstein-source)/(1+abs(einstein)+abs(source)))))
            # Independently differentiate pressure from the integrated v.
            pd=[-yy[3]**2 for yy in (ym2,ym1,y0,yp1,yp2)]
            pdot=(pd[0]-8*pd[1]+8*pd[3]-pd[4])/(12*h)
            vprime=first[3]
            exchange=2*y0[3]*vprime
            wards.append(abs(pdot+exchange)/(1+abs(pdot)+abs(exchange)))
        rows.append({'step':h,'Einstein_max_scaled_residual':max(errs),
                     'Ward_max_scaled_residual':float(max(wards))})
    return rows


def run_model(r0,pilot=False):
    s=np.linspace(0,S_MAX,721)
    inner=inner_solution(r0)
    data=inner.sol(s);r,logD,x,affine=data
    D=np.exp(logD);q=q_from_logD(r,logD)
    v=np.exp(-logD/2);u=-q*v
    inner_inv=invariants(r,u,v)
    # Same proper-coordinate system, separately integrated toward either side.
    outer=proper_solution(r0,rstop=OUTER_FACTORS[-1]*r0)
    inside=proper_solution(r0,sstop=3.,max_step=.015)
    outer_local=proper_solution(r0,rstop=8*r0,max_step=.04)
    mass_chart=exterior_mass_chart(r0,outer)
    targets=np.array(OUTER_FACTORS)*r0
    exterior_rows=[]
    for radius in targets:
        xx=brentq(lambda z:outer.sol(z)[0]-radius,0,outer.t[-1],xtol=2e-11)
        rr,rho,uu,vv=outer.sol(xx)
        mass=rr*(1-uu*uu)/2
        mass_direct,rho_direct=mass_chart.sol(np.log(rr/(2*r0)))
        madm,rho_inf=asymptotic_completion(rr,mass_direct,rho_direct)
        exterior_rows.append({'outer_radius':float(rr),'proper_x':float(xx),
                              'raw_mass':float(mass),'independent_mass':float(mass_direct),
                              'independent_rho':float(rho_direct),'ADM_estimator':float(madm),
                              'rho_at_radius':float(rho),'rho_infinity_estimator':float(rho_inf)})
    rho_inf=exterior_rows[-1]['rho_infinity_estimator']
    madm=exterior_rows[-1]['ADM_estimator']
    # Plot grid retains the smooth proper-distance neck.
    xx_outer=np.linspace(0,6*r0,481)
    yy=outer_local.sol(xx_outer)
    out_inv=invariants(yy[0],yy[2],yy[3])
    constraint_outer=yy[2]**2+2*yy[0]*yy[2]*yy[3]+yy[3]**2-1
    constraint_inner=u*u+2*r*u*v+v*v-1
    scale_inner=1+u*u+2*abs(r*u*v)+v*v
    independent=inner_solution(r0,method='Radau')
    cross_error=float(np.max(abs(independent.sol(s)-data)/(1+abs(data))))
    chart_rows=[]
    for ss in np.linspace(0,2.9,31):
        ref=inner.sol(ss);rr,ld,proper,aff=ref
        yyx=inside.sol(proper)
        expected=np.array([rr,-ss,-q_from_logD(rr,ld)*np.exp(-ld/2),np.exp(-ld/2)])
        chart_rows.append(float(np.max(abs(yyx-expected)/(1+abs(expected)))))
    refinements=[]
    for tol in TOLS:
        comparison=inner_solution(r0,rtol=tol)
        outcomparison=proper_solution(r0,rtol=tol,rstop=OUTER_FACTORS[-1]*r0)
        masscomparison=exterior_mass_chart(r0,outcomparison)
        lastmass,lastrho=masscomparison.y[:,-1]
        mm,rhoi=asymptotic_completion(OUTER_FACTORS[-1]*r0,lastmass,lastrho)
        refinements.append({'rtol':tol,'max_scaled_error_against_finest':float(np.max(abs(comparison.sol(s)-data)/(1+abs(data)))),
                            'ADM_mass':float(mm),'rho_infinity':float(rhoi),
                            'ADM_difference_from_finest':float(abs(mm-madm)),
                            'rho_infinity_difference_from_finest':float(abs(rhoi-rho_inf))})
    fd=finite_difference_checks(r0,outer_local,inside)
    lower_r2=r0*r0+s-(1-np.exp(-2*s))/2
    upper_r2=r0*r0+2*s
    lower_K=2*(1-np.exp(-2*s))**2*np.exp(8*s/3)/(r0*r0+2*s)**2
    upper_D=np.exp(-4*s/3)
    # Tail estimates are strict consequences of D<=exp(-4s/3).
    universal_proper_remaining_bound=1.5*np.exp(-2*S_MAX/3)
    proper_remaining_bound=1.5*np.sqrt(D[-1])
    affine_remaining_bound=.6*np.exp(-S_MAX)*np.sqrt(D[-1])*np.exp(-rho_inf)
    mass_chart_error=max(abs(z['raw_mass']-z['independent_mass']) for z in exterior_rows)
    rho_chart_error=max(abs(z['rho_at_radius']-z['independent_rho']) for z in exterior_rows)
    inv=inner_inv
    # Independent Einstein evaluation using the proper ODE, with quantities
    # normalized before huge inner values can cause loss of significance.
    vp=inv['rhopp'];up=inv['rpp']
    Gtt=(1-u*u-2*r*up)/(r*r)
    Gxx=(u*u-1+2*r*u*v)/(r*r)
    Gang=vp+v*v+(u*v+up)/r
    residuals=[np.max(abs(Gtt-inv['source_energy'])/(1+abs(Gtt)+abs(inv['source_energy']))),
               np.max(abs(Gxx-inv['source_radial_pressure'])/(1+abs(Gxx)+abs(inv['source_radial_pressure']))),
               np.max(abs(Gang)/(1+abs(vp)+v*v+abs(u*v/r)+abs(up/r)))]
    throat_expected_K=4*((1/(r0*r0-1))**2+2*(1/(r0*r0-1))**2+1/r0**4)
    slope=inner_rhs(S_MAX,data[:,-1])[1]
    checks={
        'proper_constraint':float(np.max(abs(constraint_outer)))<BOUNDS['constraint_scaled'],
        'inner_constraint_scaled':float(np.max(abs(constraint_inner)/scale_inner))<BOUNDS['constraint_scaled'],
        'DOP853_vs_Radau':cross_error<BOUNDS['inner_cross_method_scaled'],
        'independent_proper_chart':max(chart_rows)<BOUNDS['proper_chart_scaled'],
        'three_tolerances_converge':refinements[1]['max_scaled_error_against_finest']<refinements[0]['max_scaled_error_against_finest']/3,
        'refined_tolerance_accuracy':refinements[1]['max_scaled_error_against_finest']<BOUNDS['tolerance_refinement_scaled'],
        'outer_tolerance_convergence':refinements[1]['ADM_difference_from_finest']<refinements[0]['ADM_difference_from_finest']/3 and refinements[1]['rho_infinity_difference_from_finest']<refinements[0]['rho_infinity_difference_from_finest']/3,
        'outer_tolerance_accuracy':max(refinements[1]['ADM_difference_from_finest'],refinements[1]['rho_infinity_difference_from_finest'])<BOUNDS['tolerance_refinement_scaled'],
        'ADM_outer_radius_convergence':abs(exterior_rows[-1]['ADM_estimator']-exterior_rows[-2]['ADM_estimator'])<BOUNDS['outer_ADM_last_change'],
        'time_normalization_convergence':abs(exterior_rows[-1]['rho_infinity_estimator']-exterior_rows[-2]['rho_infinity_estimator'])<BOUNDS['outer_lapse_last_change'],
        'independent_exterior_mass_chart':mass_chart_error<BOUNDS['independent_mass_chart'] and rho_chart_error<BOUNDS['independent_mass_chart'],
        'all_Einstein_components':max(residuals)<1e-12,
        'Einstein_geometry_finite_difference':fd[-1]['Einstein_max_scaled_residual']<BOUNDS['Einstein_FD_scaled'],
        'stress_Ward_finite_difference':fd[-1]['Ward_max_scaled_residual']<BOUNDS['Ward_FD_scaled'],
        'Ward_FD_fourth_order_refinement':fd[-1]['Ward_max_scaled_residual']<fd[0]['Ward_max_scaled_residual']/100,
        'Einstein_FD_refinement':fd[-1]['Einstein_max_scaled_residual']<fd[0]['Einstein_max_scaled_residual']/20,
        'analytic_radius_lower_bound':float(np.min(r*r-lower_r2))>-BOUNDS['analytic_bound_slack'],
        'analytic_radius_upper_bound':float(np.max(r*r-upper_r2))<BOUNDS['analytic_bound_slack'],
        'analytic_D_upper_bound':float(np.max(D-upper_D))<BOUNDS['analytic_bound_slack'],
        'analytic_K_lower_bound':bool(np.all(inv['K']>=lower_K*(1-1e-11))),
        'inner_K_strictly_increases':bool(np.all(np.diff(inv['K'])>0)),
        'asymptotic_logD_slope':abs(slope+2)<BOUNDS['asymptotic_logD_slope_error'],
        'throat_lapse_positive':bool(0<np.exp(-2*rho_inf)<1),
        'throat_flaring_and_K':abs(inv['K'][0]-throat_expected_K)<1e-12 and inv['rpp'][0]>0,
        'negative_radial_null_energy':bool(np.all(inv['source_energy']+inv['source_radial_pressure']<0)),
    }
    checks={key:bool(value) for key,value in checks.items()}
    report={'throat_radius':r0,'ADM_mass':madm,'rho_infinity':rho_inf,
            'throat_lapse_squared':float(np.exp(-2*rho_inf)),
            'throat_radius_over_2ADM':r0/(2*madm),
            'throat_K':float(inv['K'][0]),
            'inner_endpoint_at_s_max':{'s':S_MAX,'r':float(r[-1]),'proper_x':float(x[-1]),
                'proper_remaining_upper_bound':proper_remaining_bound,
                'universal_proper_remaining_upper_bound':universal_proper_remaining_bound,
                'affine_interval_normalized':float(affine[-1]*np.exp(-rho_inf)),
                'affine_remaining_upper_bound':affine_remaining_bound,
                'logD':float(logD[-1]),'K':float(inv['K'][-1]),
                'K_lower_bound':float(lower_K[-1]),'logD_slope':float(slope)},
            'exterior_radius_convergence':exterior_rows,
            'independent_exterior_mass_max_difference':mass_chart_error,
            'independent_exterior_rho_max_difference':rho_chart_error,
            'constraint_outer_max':float(np.max(abs(constraint_outer))),
            'constraint_inner_scaled_max':float(np.max(abs(constraint_inner)/scale_inner)),
            'DOP853_Radau_max_scaled_difference':cross_error,
            'proper_chart_max_scaled_difference':max(chart_rows),
            'Einstein_component_scaled_residuals':list(map(float,residuals)),
            'finite_differences':fd,'tolerance_refinements':refinements,'checks':checks}
    arrays={'inner_s':s,'inner_r':r,'inner_logD':logD,'inner_q':q,'inner_proper_x':x,
            'inner_rho_normalized':-s-rho_inf,'inner_affine_normalized':affine*np.exp(-rho_inf),
            'inner_u':u,'inner_v':v,'inner_K_lower_bound':lower_K,
            'inner_r2_lower_bound':lower_r2,'inner_r2_upper_bound':upper_r2,
            'inner_D_upper_bound':upper_D,'outer_x':xx_outer,'outer_r':yy[0],
            'outer_rho_normalized':yy[1]-rho_inf,'outer_u':yy[2],'outer_v':yy[3]}
    arrays.update({'inner_'+k:val for k,val in inv.items()})
    arrays.update({'outer_'+k:val for k,val in out_inv.items()})
    return report,arrays


def plot_results(stem,models,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.23})
    colors=['#137c8b','#c95632','#6750a4']
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    for row,col in zip(models,colors):
        r0=row['throat_radius'];d=arrays[f'r{int(r0)}']
        mask=d['inner_s']<=4
        xx=np.r_[d['inner_proper_x'][mask][::-1],d['outer_x'][1:]]/r0
        rr=np.r_[d['inner_r'][mask][::-1],d['outer_r'][1:]]/r0
        rho=np.r_[d['inner_rho_normalized'][mask][::-1],d['outer_rho_normalized'][1:]]
        axs[0].plot(xx,rr,color=col,label=f'r₀={r0:g}')
        axs[1].semilogy(xx,np.exp(2*rho),color=col,label=f'r₀={r0:g}')
        axs[1].plot(0,row['throat_lapse_squared'],'o',color=col,ms=4)
    axs[0].set(xlabel='Radiale Eigendistanz x/r₀',ylabel='Flächenradius r/r₀',title='Regulärer Hals bei x=0',xlim=(-.65,4),ylim=(.95,4.8))
    axs[1].set(xlabel='Radiale Eigendistanz x/r₀',ylabel='e²ʳʰᵒ, normiert auf 1 im Unendlichen',title='Am Hals bleibt die statische Lapse positiv',xlim=(-.65,4),ylim=(1e-7,1.1))
    for ax in axs:ax.axvline(0,color='gray',ls=':',lw=1);ax.legend()
    fig.suptitle('Selbstkonsistentes Einstein–Polyakov-Modell · α=1 · Boulware-Zustand')
    save(fig,stem+'_geometry');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    for row,col in zip(models,colors):
        d=arrays[f"r{int(row['throat_radius'])}"];s=d['inner_s']
        axs[0].semilogy(s,d['inner_K'],color=col,label=f"K, r₀={row['throat_radius']:g}")
        axs[0].semilogy(s[1:],d['inner_K_lower_bound'][1:],color=col,ls='--',alpha=.65)
        axs[1].plot(s,d['inner_proper_x'],color=col,label=f"r₀={row['throat_radius']:g}")
    axs[0].set(xlabel='Innere Variable s=ρ₀−ρ',ylabel='Kretschmann-Skalar α²K',title='Krümmungswachstum; gestrichelt: bewiesene Untergrenze')
    axs[1].set(xlabel='Innere Variable s=ρ₀−ρ',ylabel='Radiale Eigendistanz x/√α',title='Endliche Eigendistanz zum singulären Rand')
    for ax in axs:ax.legend(fontsize=9)
    fig.suptitle('Halsregularität beseitigt die innere Singularität nicht')
    save(fig,stem+'_interior');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    for row,col in zip(models,colors):
        r0=row['throat_radius'];d=arrays[f'r{int(r0)}'];mask=d['inner_s']<=3
        xx=np.r_[d['inner_proper_x'][mask][::-1],d['outer_x'][1:]]/r0
        ene=np.r_[d['inner_source_energy'][mask][::-1],d['outer_source_energy'][1:]]
        nec=ene+np.r_[d['inner_source_radial_pressure'][mask][::-1],d['outer_source_radial_pressure'][1:]]
        axs[0].plot(xx,ene,color=col,label=f'ε, r₀={r0:g}')
        axs[0].plot(xx,nec,color=col,ls='--',alpha=.65)
        conv=row['exterior_radius_convergence'];R=[z['outer_radius']/r0 for z in conv[:-1]]
        err=[abs(z['independent_mass']-conv[-1]['ADM_estimator']) for z in conv[:-1]]
        axs[1].loglog(R,err,'o-',color=col,label=f'r₀={r0:g}')
    axs[0].set(yscale='symlog',xlim=(-.65,3),xlabel='Radiale Eigendistanz x/r₀',ylabel='8πGα × Quelldichte',title='Negative Energie; gestrichelt: ε+pᵣ')
    axs[1].set(xlabel='Äußerer Kontrollradius R/r₀',ylabel='|M(R) − ADM-Schätzer| / √α',title='Unabhängige Kontrolle des asymptotischen Abschlusses')
    for ax in axs:ax.legend(fontsize=9)
    fig.suptitle('Vollständiger effektiver Quellabschluss und Außenkonvergenz')
    save(fig,stem+'_source');plt.close(fig)
    return matplotlib.__version__


def save(fig,name):
    fig.savefig(ROOT/'figures'/(name+'.png'),dpi=170)
    fig.savefig(ROOT/'figures'/(name+'.svg'),metadata={'Date':None})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',action='store_true');parser.add_argument('--tag',default='20261018')
    args=parser.parse_args();start=time.perf_counter()
    if not args.tag.replace('_','').replace('-','').isalnum():parser.error('Invalid tag')
    if args.pilot:
        report,_=run_model(2.,pilot=True)
        print(json.dumps(report,indent=2),flush=True);return
    stem='bhback_numerics_'+args.tag;figstem='bhback_'+args.tag
    paths=[ROOT/'work'/(stem+ext) for ext in ('.json','.npz','.md')]
    paths += [ROOT/'figures'/(figstem+part+ext) for part in ('_geometry','_interior','_source') for ext in ('.png','.svg')]
    if any(p.exists() for p in paths):parser.error('Historical output exists; use a new tag')
    symbolic=symbolic_audit();models=[];arrays={}
    for r0 in RADII:
        print(f'Integrating throat radius {r0:g}...',flush=True)
        report,raw=run_model(r0);models.append(report);arrays[f'r{int(r0)}']=raw
    mpl=plot_results(figstem,models,arrays)
    flattened={f'{key}_{name}':a for key,values in arrays.items() for name,a in values.items()}
    assert all(np.all(np.isfinite(a)) for a in flattened.values())
    np.savez_compressed(paths[1],**flattened)
    checks={f'r{int(row["throat_radius"])}_{key}':bool(value) for row in models for key,value in row['checks'].items()}
    report={'model':{'alpha':1.,'throat_radii':RADII,'inner_s_max':S_MAX,
        'source':'2D Boulware Polyakov stress divided by 4*pi*r^2; zero angular pressure',
        'scope':'Static, self-consistent effective spherical model; no evaporation; not full 4D scalar RSET',
        'time_normalization':'rho infinity=0; initial throat rho=0 shifted afterward',
        'affine_normalization':'unit Killing energy E=1 after exterior time normalization'},
        'acceptance_bounds':BOUNDS,'symbolic':symbolic,'models':models,'checks':checks,
        'all_checks_pass':all(checks.values()),'provenance':{
        'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,
        'matplotlib':mpl,'platform':platform.platform(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arrays_sha256':hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        'elapsed_seconds':time.perf_counter()-start,
        'command':f'python work/bhback_numerics.py --tag {args.tag}'}}
    paths[0].write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    passed=sum(checks.values());total=len(checks)
    text=f'''# Statisches Einstein–Polyakov-Modell mit Rückwirkung

Lauf `{args.tag}`: **{passed}/{total}** vorab festgelegte Kontrollen bestanden.
Finale Ergebnisse dieses Laufs: JSON `{stem}.json`, Roharrays `{stem}.npz`.
Das Modell hat α=1 und Halsradien 2, 4 und 10. Die Boulware-Quelle ist der
zweidimensionale Polyakov-Tensor, geteilt durch 4πr², bei verschwindendem
Winkeldruck. Alle Einstein-Komponenten und die Quellenerhaltung werden geprüft.

Ein regulärer Hals ersetzt den klassischen Schwarzschild-Horizont. Die Lapse
bleibt dort positiv. Auf der inneren Seite wächst der Flächenradius erneut,
während die Krümmung an einem Rand in endlicher Eigendistanz divergiert.
Der Endpunkt selbst wird nicht numerisch durchquert; bewiesene Restschranken
kontrollieren die verbleibende Eigen- und Nullaffindistanz.

Eine unabhängige Integration in radialer Eigendistanz prüft den stabilen
inneren Variablenwechsel. DOP853 und Radau, drei Integratortoleranzen, vier
Außenradien und Finite-Differenzen der Metrik prüfen getrennte Fehlerquellen.
Der asymptotische ADM- und Zeitabschluss wird über Radiuskonvergenz kontrolliert.

Dies ist ein statisches, selbstkonsistentes effektives Modell. Es beschreibt
keine Unruh-Verdampfung und verwendet nicht den vollständigen vierdimensionalen
Skalar-Stresstensor. Zahlen sind theoretische Modellrechnungen, keine Messdaten.

Quell-SHA256: `{report['provenance']['source_sha256']}`.
Roharray-SHA256: `{report['provenance']['arrays_sha256']}`.
Reproduktion: `python work/bhback_numerics.py --tag NEUER_TAG`.
'''
    paths[2].write_text(text,encoding='utf-8')
    print(json.dumps({'passed':passed,'total':total,'models':models,'provenance':report['provenance']},indent=2),flush=True)
    if not report['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
