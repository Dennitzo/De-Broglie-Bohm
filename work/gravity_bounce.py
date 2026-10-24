"""Radiation-clock Bohmian minisuperspace with a Neumann half-line Hamiltonian.

Run: python work/gravity_bounce.py --tag 20261017
The rescaled constraint is -p_a**2/4+p_T=0. In the lapse N=a,
i*d_T psi=+d_a**2 psi/4 and da/dT=-d_a S/2. The Gaussian solution
is independently checked against cosine spectral propagation and ODEs.
This is a symmetry-reduced quantum geometry, not a full GR quantization.
Only new gravity_bounce_* artifacts are written; no PDF is generated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_legendre
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
T0 = 1.0
AB = 1.0
THRESHOLDS = {'spectral_wave':2e-11,'schrodinger_fd':2e-6,
              'continuity_fd':2e-7,'guidance_fd':2e-9,
              'rk4_trajectory':2e-6,'cosmic_trajectory':3e-9,
              'cosmic_clock':3e-10,'normalization':2e-11,
              'constraint':2e-12,'curvature_fd':3e-6,
              'continuity_effective_fd':2e-7,'jacobian':2e-13}


def psi(a,T,T0=T0):
    z=T0-1j*np.asarray(T)
    return (8*T0/np.pi)**.25/z**.5*np.exp(-np.asarray(a)**2/z)


def phase(a,T,T0=T0):
    return .5*np.arctan(np.asarray(T)/T0)-np.asarray(a)**2*T/(T0*T0+np.asarray(T)**2)


def guidance(T,a,T0=T0):
    return np.asarray(a)*T/(T0*T0+T*T)


def scale(T,ab=AB,T0=T0):
    return ab*np.sqrt(1+(np.asarray(T)/T0)**2)


def proper_time(T,ab=AB,T0=T0):
    T=np.asarray(T)
    return ab/(2*T0)*(T*np.sqrt(T0*T0+T*T)+T0*T0*np.arcsinh(T/T0))


def geometry(T,ab=AB,T0=T0):
    T=np.asarray(T)
    a=scale(T,ab,T0)
    A,B=ab*ab/(T0*T0),ab**4/(T0*T0)
    hubble=T/(a*(T0*T0+T*T))
    hubble_dot=-2*A/a**4+3*B/a**6
    acceleration=-A/a**4+2*B/a**6
    ricci=6*B/a**6
    kretschmann=12*(acceleration**2+hubble**4)
    rho_r=A/a**4
    rho_q=-B/a**6
    return dict(a=a,H=hubble,Hdot=hubble_dot,acceleration=acceleration,
                R=ricci,K=kretschmann,rho_r=rho_r,rho_q=rho_q,
                rho_total=rho_r+rho_q,p_total=rho_r/3+rho_q)


def rk4(Tstart,Tend,initial,step):
    count=int(round(abs(Tend-Tstart)/step))
    ts=np.linspace(Tstart,Tend,count+1)
    dt=(Tend-Tstart)/count
    out=np.empty((count+1,len(initial)))
    out[0]=initial
    for i,t in enumerate(ts[:-1]):
        q=out[i]
        k1=guidance(t,q)
        k2=guidance(t+dt/2,q+dt*k1/2)
        k3=guidance(t+dt/2,q+dt*k2/2)
        k4=guidance(t+dt,q+dt*k3)
        out[i+1]=q+dt*(k1+2*k2+2*k3+k4)/6
    return ts,out


def spectral_check(arrays):
    a=np.linspace(0,8,321)
    times=np.array([-4.,-2.,-.2,0.,.2,2.,4.])
    exact=np.stack([psi(a,t) for t in times])
    cutoff=12/np.sqrt(T0)
    initial_C=(8/(np.pi*T0))**.25
    rows=[]
    for order in (24,48,96,192):
        x,w=roots_legendre(order)
        p=(x+1)*cutoff/2
        weights=w*cutoff/2
        # Initial cosine transform; subsequent propagation uses only spectrum.
        coefficients=initial_C*np.sqrt(T0/2)*np.exp(-T0*p*p/4)
        propagated=np.sqrt(2/np.pi)*(np.exp(1j*times[:,None]*p[None,:]**2/4)*
                                     (weights*coefficients)[None,:]) @ np.cos(p[:,None]*a[None,:])
        rows.append({'nodes':order,'max_complex_error':float(np.max(abs(propagated-exact)))})
    arrays.update(spectral_a=a,spectral_T=times,spectral_psi=propagated,spectral_exact=exact)
    return {'method':'Half-line cosine spectral evolution exp(+i p^2 T/4)',
            'momentum_cutoff':float(cutoff),'phase_fitted':False,
            'normalization_adjusted':False,'convergence':rows}


def checks(arrays):
    spectral=spectral_check(arrays)
    agrid=np.array([.05,.2,.8,1.3,2.4])
    tgrid=np.array([-3.,-1.,-.2,0.,.2,1.,3.])
    rows=[]
    for h in (.002,.001,.0005):
        sch,cont,guide=0.,0.,0.
        for T in tgrid:
            ps=psi(agrid,T)
            dT=(psi(agrid,T+h)-psi(agrid,T-h))/(2*h)
            d2a=(psi(agrid+h,T)-2*ps+psi(agrid-h,T))/h**2
            # Absolute residual: tiny tails are not divided by tiny |psi|.
            sch=max(sch,float(np.max(abs(1j*dT-d2a/4))))
            rho_t=(abs(psi(agrid,T+h))**2-abs(psi(agrid,T-h))**2)/(2*h)
            jplus=abs(psi(agrid+h,T))**2*guidance(T,agrid+h)
            jminus=abs(psi(agrid-h,T))**2*guidance(T,agrid-h)
            cont=max(cont,float(np.max(abs(rho_t+(jplus-jminus)/(2*h)))))
            # Independent phase-gradient current from complex finite differences.
            deriv=(psi(agrid+h,T)-psi(agrid-h,T))/(2*h)
            g=-.5*np.imag(deriv/ps)
            guide=max(guide,float(np.max(abs(g-guidance(T,agrid)))))
        rows.append({'h':h,'schrodinger_absolute_residual':sch,
                     'continuity_absolute_residual':cont,'guidance_fd_error':guide})
    # High accuracy phase finite difference at a separate fixed step.
    h=1e-5
    phase_error=max(float(np.max(abs(-.5*(phase(agrid+h,T)-phase(agrid-h,T))/(2*h)-guidance(T,agrid)))) for T in tgrid)
    normalizations=[quad(lambda a:abs(psi(a,T))**2,0,np.inf,epsabs=2e-13,epsrel=2e-13)[0] for T in tgrid]
    boundary_current=max(abs(-.5*np.imag(np.conj(psi(0.,T))*(psi(1e-5,T)-psi(-1e-5,T))/(2e-5))) for T in tgrid)
    branch_change=float(np.max(abs(psi(agrid,1e-8)-psi(agrid,-1e-8))))
    abs_at_bounce=np.array([.2,.5,1.,1.8])
    trajectories=[]
    for step in (.08,.04,.02):
        ts,q=rk4(-4.,4.,scale(-4.)*abs_at_bounce,step)
        ref=scale(ts)[:,None]*abs_at_bounce[None,:]
        error=float(np.max(abs(q-ref)))
        trajectories.append({'dt':step,'max_absolute_error':error})
        arrays[f'rk4_T_{step:g}']=ts
        arrays[f'rk4_a_{step:g}']=q
    # Integrate in proper time with T as additional dependent clock variable.
    Tmax=4.
    tmax=float(proper_time(Tmax))
    def cosmic_rhs(t,y):
        a,T=y
        return [T/(T0*T0+T*T),1/a]
    Tref=np.linspace(-Tmax,Tmax,401)
    tref=proper_time(Tref)
    sol=solve_ivp(cosmic_rhs,(-tmax,tmax),[float(scale(-Tmax)),-Tmax],
                  method='DOP853',rtol=3e-12,atol=3e-13,t_eval=tref)
    assert sol.success
    cosmic_a_error=float(np.max(abs(sol.y[0]-scale(Tref))))
    cosmic_clock_error=float(np.max(abs(sol.y[1]-Tref)))
    arrays.update(cosmic_t=tref,cosmic_T_integrated=sol.y[1],cosmic_a_integrated=sol.y[0])
    hj_error,friedmann_error,energy_fd,curvature_fd,jac_error=0.,0.,0.,0.,0.
    derivative_h=2e-5
    for T in np.linspace(-4,4,81):
        a=scale(T)
        D=T0*T0+T*T
        pa=-2*a*T/D
        pT=T0/(2*D)-a*a*(T0*T0-T*T)/D**2
        Q=-T0/(2*D)+a*a*T0*T0/D**2
        hj_error=max(hj_error,abs(-pa*pa/4+pT+Q))
        g=geometry(T)
        friedmann_error=max(friedmann_error,abs(g['H']**2-g['rho_total']))
        gp,gm=geometry(T+derivative_h),geometry(T-derivative_h)
        rho_dot=(gp['rho_total']-gm['rho_total'])/(2*derivative_h*a)
        energy_fd=max(energy_fd,abs(rho_dot+3*g['H']*(g['rho_total']+g['p_total'])))
        Hdot_fd=(gp['H']-gm['H'])/(2*derivative_h*a)
        curvature_fd=max(curvature_fd,abs(6*(Hdot_fd+2*g['H']**2)-g['R']))
        initial=np.array([.1,.5,1.,2.])
        ratio=scale(T)
        jac_error=max(jac_error,float(np.max(abs(abs(psi(initial*ratio,T))**2*ratio-abs(psi(initial,0.))**2))))
    tclock_error=max(abs(proper_time(T)-quad(lambda u:scale(u),0,T,epsabs=2e-13,epsrel=2e-13)[0]) for T in (-4.,-1.,1.,4.))
    # Positive radiation clock operator: <p_T>=-<H>=1/(4T0).
    mean_A=quad(lambda a:(a*a/T0**2)*abs(psi(a,0))**2,0,np.inf,epsabs=2e-13)[0]
    quantum_clock_energy=quad(lambda a:(4*a*a/T0**2)*abs(psi(a,0))**2/4,0,np.inf,epsabs=2e-13)[0]
    gplot=geometry(np.linspace(-4,4,801))
    result={'spectral':spectral,'finite_differences':rows,'phase_guidance_error':phase_error,
            'normalization_max_error':float(max(abs(np.array(normalizations)-1))),
            'neumann_boundary_current':float(boundary_current),'wave_change_across_T_zero':branch_change,
            'rk4':trajectories,'proper_time_ODE_scale_error':cosmic_a_error,
            'proper_time_ODE_clock_error':cosmic_clock_error,'clock_quadrature_error':float(tclock_error),
            'quantum_HJ_constraint_error':float(hj_error),'modified_friedmann_error':float(friedmann_error),
            'effective_continuity_fd_error':float(energy_fd),'ricci_fd_error':float(curvature_fd),
            'density_jacobian_error':float(jac_error),'initial_mean_radiation_A':float(mean_A),
            'positive_radiation_clock_operator_energy':float(quantum_clock_energy),
            'minimum_scale':float(gplot['a'].min()),'maximum_R':float(gplot['R'].max()),
            'maximum_K':float(gplot['K'].max()),
            'NEC_at_bounce':float(geometry(0.)['rho_total']+geometry(0.)['p_total'])}
    flags={
        'cosine_spectral_evolution':spectral['convergence'][-1]['max_complex_error']<THRESHOLDS['spectral_wave'],
        'finite_difference_schrodinger':rows[-1]['schrodinger_absolute_residual']<THRESHOLDS['schrodinger_fd'],
        'finite_difference_continuity':rows[-1]['continuity_absolute_residual']<THRESHOLDS['continuity_fd'],
        'finite_difference_convergence':rows[-1]['schrodinger_absolute_residual']<rows[0]['schrodinger_absolute_residual']/10,
        'phase_gradient_guidance':phase_error<THRESHOLDS['guidance_fd'],
        'normalization':result['normalization_max_error']<THRESHOLDS['normalization'],
        'neumann_zero_boundary_flux':boundary_current == 0,
        'rk4_guidance':trajectories[-1]['max_absolute_error']<THRESHOLDS['rk4_trajectory'],
        'rk4_convergence':trajectories[-1]['max_absolute_error']<trajectories[0]['max_absolute_error']/100,
        'proper_time_independent_ODE':cosmic_a_error<THRESHOLDS['cosmic_trajectory'] and cosmic_clock_error<THRESHOLDS['cosmic_trajectory'],
        'proper_time_integral':tclock_error<THRESHOLDS['cosmic_clock'],
        'quantum_hamilton_jacobi':hj_error<THRESHOLDS['constraint'],
        'modified_friedmann':friedmann_error<THRESHOLDS['constraint'],
        'effective_source_continuity':energy_fd<THRESHOLDS['continuity_effective_fd'],
        'independent_curvature_derivative':curvature_fd<THRESHOLDS['curvature_fd'],
        'density_jacobian':jac_error<THRESHOLDS['jacobian'],
        'positive_minimum_scale':result['minimum_scale']==AB,
        'curvature_at_bounce':abs(result['maximum_R']-6/(AB*T0)**2)<1e-12 and abs(result['maximum_K']-12/(AB*T0)**4)<1e-12,
        'effective_NEC_violated_at_bounce':result['NEC_at_bounce']<0,
        'clock_energy_and_ensemble_asymptote':abs(mean_A-1/(4*T0))<2e-12 and abs(quantum_clock_energy-mean_A)<2e-12,
    }
    flags={name:bool(value) for name,value in flags.items()}
    result['acceptance_checks']=flags
    result['all_checks_pass']=all(flags.values())
    return result


def plot(tag,arrays):
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','svg.hashsalt':'bohm-gravity-bounce-v1'})
    T=np.linspace(-4,4,801)
    t=proper_time(T)
    g=geometry(T)
    agrid=np.linspace(0,5,301)
    density=abs(psi(agrid[:,None],T[None,:]))**2
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    ax=axes[0,0]
    m=ax.pcolormesh(T,agrid,density,shading='auto',cmap='Blues',rasterized=True,vmin=0,vmax=1.6)
    for ab in (.2,.4,.8,1.2):
        ax.plot(T,scale(T,ab),color='#b43b2d',lw=1.6)
    ax.set(xlabel='Strahlungsuhr / konforme Zeit T',ylabel='Skalenfaktor a',ylim=(0,5),
           title='Quantendichte und Bohmsche Geometrien')
    fig.colorbar(m,ax=ax,label=r'$|\psi(a,T)|^2$')
    ax=axes[0,1]
    ax.plot(t,g['a'],color='#17628d',lw=2.2,label='Bohmsche Geometrie')
    ax.plot(t,np.sqrt(2*AB*abs(t)/T0),'--',color='#b43b2d',label='Klassische Strahlungszweige')
    ax.axhline(AB,color='.45',ls=':',label=r'Minimum $a_b=1$')
    ax.set(xlabel='Kosmische Eigenzeit t',ylabel='Skalenfaktor a',title='Bounce und klassischer Grenzbereich')
    ax.legend(fontsize=9)
    ax=axes[1,0]
    ax.plot(t,g['rho_r'],label=r'$\rho_r=A/a^4$',color='#17628d')
    ax.plot(t,g['rho_q'],label=r'$\rho_Q=-B/a^6$',color='#b43b2d')
    ax.plot(t,g['rho_total'],label=r'$\rho_{\rm eff}=H^2$',color='#267d4a',lw=2)
    ax.plot(t,g['rho_total']+g['p_total'],'--',label=r'$\rho_{\rm eff}+p_{\rm eff}$',color='#7556a6')
    ax.axhline(0,color='.6',lw=.6)
    ax.set(xlabel='Kosmische Eigenzeit t',ylabel='Effektive Dichten',title='Erhaltene geometrische Quellenbilanz')
    ax.legend(fontsize=9)
    ax=axes[1,1]
    ax.semilogy(t,g['R'],color='#17628d',lw=2,label=r'$\mathcal{R}=6B/a^6$')
    ax.semilogy(t,np.sqrt(g['K']),color='#b43b2d',lw=2,label=r'$\sqrt{\mathcal{K}}$')
    ax.set(xlabel='Kosmische Eigenzeit t',ylabel='Krümmung (reskalierte Einheiten)',title='Endliche Krümmungsinvarianten')
    ax.legend(fontsize=9)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Bohmsche Quantengeometrie: flaches Strahlungs-Minisuperspace\n'
                 r'$T_0=1$, Referenzbahn $a_b=1$; effektive Quellen in Einheiten $8\pi G/3=1$',fontsize=15)
    png=ROOT/f'figures/gravity_bounce_{tag}.png'
    svg=ROOT/f'figures/gravity_bounce_{tag}.svg'
    fig.savefig(png,dpi=170,metadata={'Software':'gravity_bounce.py'})
    fig.savefig(svg,metadata={'Date':None})
    plt.close(fig)
    arrays.update(plot_T=T,plot_t=t,plot_a_grid=agrid,plot_density=density,
                  plot_scale=g['a'],plot_R=g['R'],plot_K=g['K'],
                  plot_rho_r=g['rho_r'],plot_rho_q=g['rho_q'],plot_p_total=g['p_total'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='20261017')
    parser.add_argument('--pilot',action='store_true',help='Print checks without saving any artifact')
    args=parser.parse_args()
    paths=[ROOT/f'work/gravity_bounce_{args.tag}.{ext}' for ext in ('json','npz','md')]
    paths += [ROOT/f'figures/gravity_bounce_{args.tag}.{ext}' for ext in ('png','svg')]
    if not args.pilot and any(p.exists() for p in paths):parser.error('Historical output exists; choose another tag.')
    arrays={}
    result=checks(arrays)
    result['model']={'T0':T0,'reference_ab':AB,'half_line':'a>0; Neumann boundary at a=0',
                     'clock':'T is radiation time; proper time dt=a dT',
                     'hamiltonian':'i dT psi=+1/4 da^2 psi',
                     'units':'rescaled minisuperspace; effective densities use 8piG/3=1',
                     'source':'Pinto-Neto and Fabris 2013, section 4.2.1; explicit radiation specialization'}
    result['acceptance_thresholds']=THRESHOLDS
    if args.pilot:
        print(json.dumps(result,indent=2))
        return
    plot(args.tag,arrays)
    np.savez_compressed(paths[1],**arrays)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    result['provenance']={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,
                          'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
                          'source_sha256':digest(Path(__file__)),'arrays_sha256':digest(paths[1]),
                          'png_sha256':digest(paths[3]),'svg_sha256':digest(paths[4])}
    paths[0].write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    lines=['# Radiation-clock minisuperspace bounce',
           f"Passed {sum(result['acceptance_checks'].values())}/{len(result['acceptance_checks'])} checks.",
           '',f"Spectral complex error: {result['spectral']['convergence'][-1]['max_complex_error']:.12g}",
           f"Proper-time ODE scale error: {result['proper_time_ODE_scale_error']:.12g}",
           f"Density-Jacobian error: {result['density_jacobian_error']:.12g}",
           f"Source SHA256: {result['provenance']['source_sha256']}"]
    paths[2].write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'all_checks_pass':result['all_checks_pass'],'checks':result['acceptance_checks'],
                      'metrics':str(paths[0]),'plot':str(paths[3])},indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
