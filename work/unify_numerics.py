"""Exact coupled geometry/scalar packet: independent numerical evolution.

Units hbar=c=1; beta=sigma=1, fiducial comoving volume one. The state is
f(phi+alpha)+f(phi-alpha), f(s)=1/(sigma-i*s). No QFT RSET is asserted.
Run from any directory: python work/unify_numerics.py --tag 20261018_release
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SIGMA = BETA = 1.0
ALPHA0 = 0.5


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def packet(alpha, phi, sigma=SIGMA):
    """Raw packet; j-normalized state is sigma/sqrt(pi) times this."""
    return 1 / (sigma - 1j * (phi + alpha)) + 1 / (sigma - 1j * (phi - alpha))


def phase(alpha, phi, sigma=SIGMA):
    s = sigma**2 + phi**2
    d = (sigma**2 + (phi + alpha)**2) * (sigma**2 + (phi - alpha)**2)
    pa = -4 * sigma * alpha * phi / d
    pp = sigma * (s**2 + 4 * alpha**2 * phi**2 - alpha**4) / (s * d)
    return pa, pp


def primitive(s, sigma=SIGMA):
    return s / (2*sigma*(sigma**2+s**2)) + np.arctan(s/sigma)/(2*sigma**2)


def stream(alpha, phi, sigma=SIGMA):
    u, v = phi + alpha, phi - alpha
    return primitive(u, sigma) - primitive(v, sigma) + 2*sigma*alpha / ((sigma**2+u**2)*(sigma**2+v**2))


def critical(phi, sigma=SIGMA):
    return np.sqrt(2*phi**2 + np.sqrt(sigma**4+2*sigma**2*phi**2+5*phi**4))


def level_alpha(phi, level):
    return brentq(lambda a: stream(a, phi)-level, 0, critical(phi), xtol=1e-14, rtol=1e-14)


def build_functions():
    """Derive Q from amplitude, independently of the phase/HJ identity."""
    a, f = sp.symbols('a f', real=True)
    d = (1+(f+a)**2)*(1+(f-a)**2)
    log_r = sp.log(4*(1+f*f)/d)/2
    q = (sp.diff(log_r,a,2)+sp.diff(log_r,a)**2-sp.diff(log_r,f,2)-sp.diff(log_r,f)**2)/2
    # Exact rational simplification preserves symmetry zeros (notably phi=0)
    # before large proper-time factors amplify floating-point cancellation.
    q = sp.factor(q)
    qa, qf = sp.factor(sp.diff(q,a)), sp.factor(sp.diff(q,f))
    pa = -4*a*f/d
    pp = ((1+f*f)**2+4*a*a*f*f-a**4)/((1+f*f)*d)
    lapse = sp.exp(-3*a)
    h, v = -lapse*pa, lapse*pp
    hd = sp.diff(h,a)*h + sp.diff(h,f)*v
    rho = v*v/2
    rq = sp.exp(-6*a)*q
    pq = sp.exp(-6*a)*(q-qa/3)
    drho = sp.factor(sp.diff(rho,a)*h+sp.diff(rho,f)*v)
    drq = sp.factor(sp.diff(rq,a)*h+sp.diff(rq,f)*v)
    expressions = [q,qa,qf,h,v,hd,rho,rq,pq,drho,drq]
    fn = sp.lambdify((a,f),expressions,'numpy',cse=True)
    return fn


def main(tag):
    receipt = ROOT/'work'/f'unify_numerics_{tag}.json'
    arrays = receipt.with_suffix('.npz')
    if receipt.exists() or arrays.exists():
        raise FileExistsError('Choose a fresh tag; evidence is append-only.')
    fn = build_functions()
    checks, metrics = {}, {}
    def check(name, value):
        checks[name] = bool(value)
    c = float(stream(ALPHA0,0))
    phi = np.linspace(-6,6,601)
    alpha = np.array([level_alpha(f,c) for f in phi])
    pa, pp = phase(alpha,phi)
    q,qa,qf,h,v,hd,rho,rq,pq,drho,drq = [np.asarray(x) for x in fn(alpha,phi)]
    ricci = 6*(hd+2*h*h)
    curvature = 12*((hd+h*h)**2+h**4)
    exchange = -np.exp(-6*alpha)*v*qf

    def relational_rhs(f,y):
        m,n = phase(y[0],f)
        return [-m/n, np.exp(3*y[0])/n]
    positive = phi[phi >= 0]
    guidance = solve_ivp(relational_rhs,(0,6),[ALPHA0,0],t_eval=positive,
                         method='DOP853',rtol=2e-12,atol=[2e-13,2e-12],max_step=.04)
    if not guidance.success:
        raise RuntimeError(guidance.message)
    t = np.r_[-guidance.y[1,1:][::-1], guidance.y[1]]
    guidance_error = float(np.max(np.abs(guidance.y[0]-alpha[300:])))
    metrics['implicit_vs_phase_ode_alpha_max_abs'] = guidance_error
    check('implicit_level_vs_phase_guidance',guidance_error < 2e-9)

    # Independent canonical momenta, evolved from the amplitude-derived Q.
    # Rescaled Hamiltonian C+Q, with phi used as independent variable.
    def canonical_rhs(f,y):
        a,ma,mf,_ = y
        _,dq_a,dq_f,*_ = fn(a,f)
        return [-ma/mf,-dq_a/mf,-dq_f/mf,np.exp(3*a)/mf]
    canonical = solve_ivp(canonical_rhs,(0,6),[ALPHA0,0,.6,0],t_eval=positive,
                          method='DOP853',rtol=3e-12,atol=2e-13,max_step=.035)
    if not canonical.success:
        raise RuntimeError(canonical.message)
    ca,cm,cn,ct = canonical.y
    phase_ca,phase_cf=phase(ca,positive)
    cq=fn(ca,positive)[0]
    metrics['canonical_alpha_max_abs']=float(np.max(np.abs(ca-alpha[300:])))
    metrics['canonical_phase_graph_max_abs']=float(max(np.max(np.abs(cm-phase_ca)),np.max(np.abs(cn-phase_cf))))
    metrics['canonical_constraint_max_abs']=float(np.max(np.abs((-cm*cm+cn*cn)/2+cq)))
    metrics['canonical_clock_scaled_max_abs']=float(np.max(np.abs(ct-guidance.y[1])/(1+np.abs(ct))))
    check('amplitude_canonical_flow_vs_phase_guidance', metrics['canonical_alpha_max_abs']<2e-8 and metrics['canonical_phase_graph_max_abs']<2e-9)
    check('canonical_constraint_preserved',metrics['canonical_constraint_max_abs']<2e-10)
    check('canonical_and_guidance_proper_clocks',metrics['canonical_clock_scaled_max_abs']<2e-8)

    # Direct proper-time integration, independent of the scalar-clock gauge.
    cut=np.searchsorted(positive,1.0)
    t_end=float(guidance.y[1,cut])
    def proper_rhs(time,y):
        ma,mf=phase(*y)
        factor=np.exp(-3*y[0])
        return [-factor*ma,factor*mf]
    proper=solve_ivp(proper_rhs,(0,t_end),[ALPHA0,0],method='DOP853',rtol=2e-12,atol=2e-13,
                     t_eval=guidance.y[1,:cut+1],max_step=t_end/100)
    if not proper.success:
        raise RuntimeError(proper.message)
    proper_error=float(max(np.max(np.abs(proper.y[0]-alpha[300:301+cut])),np.max(np.abs(proper.y[1]-positive[:cut+1]))))
    metrics['proper_vs_scalar_clock_max_abs']=proper_error
    check('positive_lapse_reparametrization',proper_error<2e-9)

    # Local equations checked on a 2D domain, not only one selected solution.
    aa,ff=np.meshgrid(np.linspace(-3.1,3.1,37),np.linspace(-2.7,2.7,35))
    m,n=phase(aa,ff)
    z,za,zf,hh,vv,hhp,rr,rrq,ppq,drr,drrq=fn(aa,ff)
    residuals={
        'hamilton_jacobi':(-m*m+n*n)/2+z,
        'friedmann':hh*hh-2*(rr+rrq),
        'raychaudhuri':hhp+3*(2*rr+rrq+ppq),
        'field_ward':drr+6*hh*rr+np.exp(-6*aa)*vv*zf,
        'quantum_ward':drrq+3*hh*(rrq+ppq)-np.exp(-6*aa)*vv*zf,
    }
    # Scale each residual by one plus its typical local term to avoid giving
    # exponential coordinate factors a misleading absolute tolerance.
    scales={
        'hamilton_jacobi':1+m*m+n*n+np.abs(z),
        'friedmann':1+hh*hh+2*np.abs(rr)+2*np.abs(rrq),
        'raychaudhuri':1+np.abs(hhp)+3*(2*np.abs(rr)+np.abs(rrq)+np.abs(ppq)),
        'field_ward':1+np.abs(drr)+6*np.abs(hh*rr)+np.abs(np.exp(-6*aa)*vv*zf),
        'quantum_ward':1+np.abs(drrq)+3*np.abs(hh*(rrq+ppq))+np.abs(np.exp(-6*aa)*vv*zf),
    }
    for key,res in residuals.items():
        err=float(np.max(np.abs(res)/scales[key]))
        metrics[key+'_scaled_max_abs']=err
        check(key,err<3e-12)

    # Independent finite-difference directional derivatives of both densities.
    # q +/- step*v(q) probes the actual first directional derivative at q.
    samples=np.array([[.5,.0],[.7,.3],[1.,.8],[1.7,1.1]])
    fd_errors=[]
    for step in [.04,.02,.01]:
        errors=[]
        for a0,f0 in samples:
            vals=fn(a0,f0)
            ha,va=vals[3:5]
            def dens(k):
                item=fn(a0+k*step*ha,f0+k*step*va)
                return np.array(item[6:8])
            deriv=(dens(-2)-8*dens(-1)+8*dens(1)-dens(2))/(12*step)
            errors.append(np.max(np.abs(deriv-np.array(vals[9:11]))))
        fd_errors.append(float(max(errors)))
    metrics['ward_density_fivepoint_errors']=fd_errors
    check('actual_derivative_convergence',fd_errors[-1]<1e-10 and fd_errors[0]/fd_errors[1]>12 and fd_errors[1]/fd_errors[2]>10)

    # Current/L2 integral checks at several scalar-clock slices.
    norms=[]
    for f0 in [0.,.5,2.,6.]:
        l2=quad(lambda a:abs(packet(a,f0))**2,-np.inf,np.inf,epsabs=2e-10,epsrel=2e-11,limit=400)[0]
        flux=quad(lambda a:abs(packet(a,f0))**2*phase(a,f0)[1],-np.inf,np.inf,epsabs=2e-10,epsrel=2e-11,limit=400)[0]
        norms.append({'phi':f0,'raw_l2':l2,'raw_flux':flux})
    metrics['norm_quadratures']=norms
    check('finite_cauchy_norm_and_positive_spectral_flux',all(abs(x['raw_l2']-2*np.pi)<2e-9 and abs(x['raw_flux']-np.pi)<2e-9 for x in norms))
    check('positive_frequency_not_pointwise_positive',phase(2.,0.)[1]<0 and np.all(pp>0))

    # Positive conditional flux: selected invariant tube, no global Born claim.
    q_values=np.linspace(.05,.95,13)
    quantiles=np.array([[level_alpha(f,float(q0*c)) for f in phi] for q0 in q_values])
    tube_errors=[]
    for i in [0,150,300,450,600]:
        value=quad(lambda a:abs(packet(a,phi[i]))**2*phase(a,phi[i])[1],0,alpha[i],epsabs=1e-11)[0]
        tube_errors.append(abs(value/c-1))
    metrics['conditional_tube_norm_max_abs']=float(max(tube_errors))
    quantile_error=max(np.max(np.abs(stream(row,phi)/c-q0)) for q0,row in zip(q_values,quantiles))
    metrics['conditional_quantile_level_max_abs']=float(quantile_error)
    check('conditional_positive_flux_equivariance',max(tube_errors)<2e-10 and quantile_error<2e-12 and np.all(np.diff(quantiles,axis=0)>0))

    vinf=brentq(lambda x:primitive(x)-(np.pi/4-c),-20,20,xtol=1e-14)
    metrics['asymptotic_v']=vinf
    farphi=np.array([100.,1000.,10000.])
    faralpha=np.array([level_alpha(f,c) for f in farphi])
    metrics['asymptotic_v_errors']=list(np.abs(farphi-faralpha-vinf))
    # Leading interference tail is O(1/phi), not O(1/phi^2).
    ratios=np.array(metrics['asymptotic_v_errors'][:-1])/np.array(metrics['asymptotic_v_errors'][1:])
    check('asymptotic_regular_branch',metrics['asymptotic_v_errors'][-1]<2e-4 and np.all((ratios>9.9)&(ratios<10.1)))
    metrics['bounce']={'alpha':ALPHA0,'a':float(np.exp(ALPHA0)),'phase_phi':float(pp[300]),'Hdot':float(hd[300]),'R':float(ricci[300]),'K':float(curvature[300]),'rho_field':float(rho[300]),'rho_quantum':float(rq[300]),'rho_plus_pressure':float((2*rho+rq+pq)[300])}
    check('regular_bounce_and_source_cancellation',np.min(alpha)>=ALPHA0-1e-13 and hd[300]>0 and abs(rho[300]+rq[300])<1e-14 and (2*rho+rq+pq)[300]<0)

    # Singular packet f(phi+alpha), same sigma and initial (alpha,phi).
    k=1/(1+ALPHA0**2)
    ts=np.exp(3*ALPHA0)/(3*k)
    singular_t=ts*(1-np.geomspace(1,1e-5,401))
    delta=ts-singular_t
    singular_a=(3*k*delta)**(1/3)
    singular_k=20/(27*delta**4)
    def chiral_rhs(time,y):
        return [-k*np.exp(-3*y[0]),k*np.exp(-3*y[0])]
    control=solve_ivp(chiral_rhs,(0,.98*ts),[ALPHA0,0],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True,max_step=ts/100)
    if not control.success:
        raise RuntimeError(control.message)
    tc=np.linspace(0,.98*ts,201)
    ac,fc=control.sol(tc)
    metrics['singular_control_volume_max_abs']=float(np.max(np.abs(np.exp(3*ac)-(np.exp(3*ALPHA0)-3*k*tc))))
    metrics['singular_time']=float(ts)
    check('same_equation_singular_control',metrics['singular_control_volume_max_abs']<2e-10 and np.max(np.abs(ac+fc-ALPHA0))<1e-12)

    data=dict(phi=phi,alpha=alpha,t=t,hubble=h,hdot=hd,ricci=ricci,kretschmann=curvature,rho_field=rho,rho_quantum=rq,pressure_quantum=pq,exchange=exchange,quantiles=quantiles,q_values=q_values,singular_t=singular_t,singular_a=singular_a,singular_k=singular_k,phase_alpha=pa,phase_phi=pp,quantum_potential=q,canonical_alpha=ca,canonical_pa=cm,canonical_pp=cn,canonical_t=ct)
    check('all_arrays_finite',all(np.all(np.isfinite(x)) for x in data.values()))
    np.savez_compressed(arrays,**data)
    sources=[Path(__file__).resolve(),ROOT/'inputs/unify_source_struyve2023.pdf',ROOT/'inputs/unify_source_pn2013.pdf']
    result={'model':'Joint WDW flat-FLRW gravity and homogeneous massless scalar; rational positive-frequency packet',
            'scope':'Reduced geometry/field quantum potential, not an absolute renormalized 4D RSET.',
            'parameters':{'sigma':SIGMA,'beta':BETA,'alpha0':ALPHA0,'stream_level':c,'stream_plateau':np.pi/2},
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
            'checks':checks,'all_passed':all(checks.values()),'metrics':metrics,
            'arrays':{'path':str(arrays.relative_to(ROOT)).replace('\\','/'),'sha256':sha(arrays)},
            'inputs':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sources}}
    receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'receipt':str(receipt),'all_passed':result['all_passed'],'checks':checks,'metrics':metrics},indent=2))
    if not result['all_passed']:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--tag',required=True)
    main(parser.parse_args().tag)
