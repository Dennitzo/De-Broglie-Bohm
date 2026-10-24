"""Absolute massive Dirac de Sitter source and order-reduced Schwarzschild response.

The exact de Sitter source is checked against subtracted Hankel mode integrals.
The Schwarzschild source is ONLY the local O(mu^-2) DeWitt approximation.
No absolute massive black-hole mode sum or singularity resolution is claimed.
"""
from pathlib import Path
import argparse, hashlib, json, math, time
from datetime import datetime, timezone
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp, quad
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def F(x):
    x = mp.mpf(x)
    if not x:
        return mp.mpf(11)
    return 11+10*x*x+120*x*x*(1+x*x)*(mp.log(x)-mp.re(mp.digamma(1+1j*x)))

def F_original(x):
    x = mp.mpf(x)
    return 11+130*x*x+120*x*x*(1+x*x)*(mp.log(x)-mp.re(mp.digamma(-1+1j*x)))

def F_series_bounds(x, nmax=2000):
    """Integral-test enclosure of the digamma series (mp arithmetic, not intervals)."""
    x = mp.mpf(x)
    s = mp.fsum(x*x/(n*(n*n+x*x)) for n in range(1, nmax+1))
    low = -mp.euler+s+mp.log1p(x*x/(nmax+1)**2)/2
    high = -mp.euler+s+mp.log1p(x*x/nmax**2)/2
    a = 11+10*x*x
    b = 120*x*x*(1+x*x)
    return a+b*(mp.log(x)-high), a+b*(mp.log(x)-low)

def bisect(fun, lo, hi, steps=170):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    flo, fhi = fun(lo), fun(hi)
    assert flo*fhi < 0
    for _ in range(steps):
        mid=(lo+hi)/2
        fm=fun(mid)
        if flo*fm > 0:
            lo,flo=mid,fm
        else:
            hi,fhi=mid,fm
    return (lo+hi)/2, [str(lo),str(hi)]

def mode_integrands(z, x):
    """H=1, source equations 85/86,89,91--93; pressure from actual mode."""
    z, x = mp.mpf(z), mp.mpf(x)
    nu=mp.mpf('.5')-1j*x
    pref=mp.sqrt(mp.pi*z)*mp.exp(mp.pi*x/2)/2
    u=1j*pref*mp.hankel1(nu,z)
    v=pref*mp.hankel1(nu-1,z)
    s1=2*mp.re(mp.conj(u)*v)
    s3=abs(u)**2-abs(v)**2
    bare_r=-2*(x*s3+z*s1)
    bare_p=-2*z*s1/3
    w2=z*z+x*x
    r0=-2*mp.sqrt(w2)
    p0=-2*z*z/(3*mp.sqrt(w2))
    terms2=[(x*x/4,mp.mpf(3)/2),(-x**4/4,mp.mpf(5)/2)]
    terms4=[(105*x**8/64,mp.mpf(11)/2),(-119*x**6/32,mp.mpf(9)/2),
            (165*x**4/64,mp.mpf(7)/2),(-x*x/2,mp.mpf(5)/2)]
    rn=mp.fsum(c*w2**(-p) for c,p in terms2+terms4)
    pn=mp.fsum(-2*p*z*z*c*w2**(-p-1)/3 for c,p in terms2+terms4)
    factor=z*z/(2*mp.pi**2)
    return factor*(bare_r-r0-rn),factor*(bare_p-p0-pn),abs(abs(u)**2+abs(v)**2-1)

def hankel_quadrature(x, order, cutoffs=(32,64,128)):
    """Finite cutoffs and order comparison; extrapolation is explicitly diagnostic."""
    x=mp.mpf(x)
    gx,gw=leggauss(order)
    edges=sorted(set([mp.mpf(0),x/8,x/2,x,2*x]+[mp.mpf(z) for z in [1,4,8,16,*cutoffs]]))
    edges=[z for z in edges if z<=max(cutoffs)]
    total=[mp.mpf(0),mp.mpf(0)]
    norm=mp.mpf(0);out={}
    for a,b in zip(edges[:-1],edges[1:]):
        vals=[mp.mpf(0),mp.mpf(0)]
        for xi,wi in zip(gx,gw):
            z=(a+b)/2+(b-a)*mp.mpf(float(xi))/2
            rr,pp,nn=mode_integrands(z,x)
            norm=max(norm,nn)
            weight=(b-a)*mp.mpf(float(wi))/2
            vals[0]+=weight*rr;vals[1]+=weight*pp
        total=[total[k]+vals[k] for k in range(2)]
        if b in cutoffs:
            out[str(int(b))]=[float(q) for q in total]
    # E=z sum c_k z^-2k, c_(k+1)=(1-2k)(k²+x²)c_k/[2(k+1)].
    # After adiabatic order four subtraction the mode density starts
    # 45 x²/(16 z^7); hence the integrated leading tail is O(cutoff^-4).
    # This asymptotic Richardson diagnostic is not a rigorous tail bound.
    z1,z2=cutoffs[-2:]
    rich=[(z2**4*out[str(z2)][k]-z1**4*out[str(z1)][k])/(z2**4-z1**4) for k in range(2)]
    return dict(order=order,cutoffs=out,richardson=rich,max_mode_norm_error=float(norm))

def massive_desitter(mode_check):
    mp.mp.dps=55
    xs=np.geomspace(.003,50,181)
    values=np.array([float(F(str(x))) for x in xs])
    checks={}
    recurrence=max(float(abs(F(str(x))-F_original(str(x)))) for x in [.03,.1,.5,1,10,50])
    checks['digamma_recurrence']=recurrence<1e-40
    series=[]
    for x in [.05,.1,.5,1,2]:
        lo,hi=F_series_bounds(str(x));v=F(str(x))
        assert lo<v<hi
        series.append(dict(x=x,low=str(lo),high=str(hi),value=str(v)))
    checks['independent_positive_series']=True
    root,bracket=bisect(F,'.4','.6')
    eta=mp.mpf(1)
    xr,bxr=bisect(lambda x: eta*F(x)-360*mp.pi*x*x,'.05','.15')
    hr=1/xr
    ell=mp.mpf('.03')
    equation=lambda h: 3*h*h-ell-eta*h**4*F(1/h)/(120*mp.pi)
    hl,bhl=bisect(equation,'.099','.101')
    root_records=[]
    for label,h,l,br in [('vacuum_Lambda0',hr,mp.mpf(0),bxr),('Lambda_positive_low_curvature',hl,ell,bhl)]:
        rho=h**4*F(1/h)/(960*mp.pi**2)
        residual=3*h*h-l-8*mp.pi*eta*rho
        root_records.append(dict(label=label,eta=str(eta),Lambda_over_mu2=str(l),H_over_mu=str(h),
             rho_over_Nmu4=str(rho),NGH2=str(eta*h*h),GH2_for_N100=str(eta*h*h/100),
             R_over_mu2=str(12*h*h),K_over_mu4=str(24*h**4),residual=str(residual),root_bracket=br,
             bracket_variable='x=mu/H' if l==0 else 'h=H/mu'))
    checks['selfconsistent_roots']=all(abs(mp.mpf(r['residual']))<mp.mpf('1e-45') for r in root_records)
    s20=mp.fsum(mp.mpf(1)/(n*(n*n+1)) for n in range(1,21))
    upper=21+240*(mp.mpf(29)/50-s20)
    checks['exact_existence_sign_bound']=upper<0
    small_limit=float(F('1e-8'))
    asym_x=mp.mpf(100)
    asym_error=float(abs(asym_x**4*(F(asym_x)+mp.mpf(31)/(21*asym_x**2))+mp.mpf(41)/42))
    checks['massless_limit']=abs(small_limit-11)<1e-11
    checks['massive_asymptotic']=asym_error<.001
    modes=[]
    if mode_check:
        for x in ['0.09050308440870979379','1.0','10.0']:
            exact=float(F(x)/(960*mp.pi**2))
            # The heavy field has a much larger ultraviolet remainder. Extend
            # the actual integral instead of loosening the acceptance threshold.
            scale=max(1,math.ceil(float(x)/5))
            cuts=tuple(scale*z for z in (128,256,512))
            coarse=hankel_quadrature(x,32,cuts)
            fine=hankel_quadrature(x,64,cuts)
            err=abs(fine['richardson'][0]-exact)
            pref=max(abs(exact),1e-8)
            record=dict(x=x,exact_rho_over_H4=exact,coarse=coarse,fine=fine,
                relative_rho_error=err/pref,relative_pressure_error=abs(fine['richardson'][1]+exact)/pref,
                relative_order_change=abs(fine['richardson'][0]-coarse['richardson'][0])/pref)
            modes.append(record)
            print('mode integral',x,json.dumps({k:v for k,v in record.items() if 'relative' in k}),flush=True)
        checks['subtracted_mode_integral']=all(r['relative_rho_error']<1e-7 and r['relative_pressure_error']<1e-7 for r in modes)
        checks['mode_quadrature_refinement']=all(r['relative_order_change']<2e-6 for r in modes)
        checks['mode_normalization']=all(r['fine']['max_mode_norm_error']<1e-40 for r in modes)
    hh=np.geomspace(.04,20,201)
    curve=np.array([float(mp.mpf(str(h))**4*F(1/mp.mpf(str(h)))/(120*mp.pi)) for h in hh])
    result=dict(model='exact massive Dirac Bunch-Davies source, fourth order adiabatic renormalization',
        source_equation='del Rio/Navarro-Salas/Torrenti2014 Eq96, signature converted to (-+++)',
        recurrence_error=recurrence,series_enclosures=series,sign_change_x=str(root),sign_bracket=bracket,
        F_at_one=str(F(1)),F_at_one_rational_upper_bound=str(upper),roots=root_records,
        large_mass_asymptotic_error=asym_error,modes=modes,checks=checks)
    return result,dict(ds_x=xs,ds_F=values,ds_h=hh,ds_quantum_term=curve)

def schwarzschild():
    # dimensionless radius x=r/M, source scale N/(mu² M⁶), epsilon=NG/(mu²M⁴).
    x=np.geomspace(2,40,501)
    tt=lambda r:(2384-1080*r)/(20160*math.pi**2*r**9)
    rr=lambda r:(504*r-784)/(20160*math.pi**2*r**9)
    tang=lambda r:(3536-1512*r)/(20160*math.pi**2*r**9)
    dmass=lambda r:(149-81*r)/(1890*math.pi*r**6)
    chi=lambda r:-11/(210*math.pi*r**6)
    massprime=lambda r:-4*math.pi*r*r*tt(r)
    chiprime=lambda r:11/(35*math.pi*r**7)
    # Independent integral from infinity and actual ODE, no fitted constants.
    integral_errors=[]
    for r in [2,2.2,3,5,10,40]:
        mi=-quad(massprime,r,np.inf,epsabs=1e-16,epsrel=1e-12)[0]
        ci=-quad(chiprime,r,np.inf,epsabs=1e-16,epsrel=1e-12)[0]
        integral_errors += [abs(mi-dmass(r)),abs(ci-chi(r))]
    sol=solve_ivp(lambda r,y:[massprime(r),chiprime(r)],(40,2),[dmass(40),chi(40)],
        method='DOP853',rtol=2e-12,atol=1e-17,dense_output=True,max_step=.2)
    ode=sol.sol(x)
    odeerr=float(np.max(abs(ode-np.array([dmass(x),chi(x)]))))
    # Derivative Ward check by an independent five-point stencil.
    ward=[]
    for h in [2e-3,1e-3,5e-4]:
        r=np.linspace(2.01,8,80)
        dp=(-rr(r+2*h)+8*rr(r+h)-8*rr(r-h)+rr(r-2*h))/(12*h)
        f=1-2/r;fp=2/r**2
        rest=dp+fp/(2*f)*(rr(r)-tt(r))+2/r*(rr(r)-tang(r))
        ward.append(float(np.max(abs(rest))))
    eps=1e-4;muM=10;NGM2=.01
    shift=-13*eps/(60480*math.pi)
    # Only a diagnostic of the local expansion, not an error bound on unknown state terms.
    local_curvature=np.sqrt(48)/(muM**2*x**3)
    checks=dict(infinity_quadratures=max(integral_errors)<1e-15,
        independent_metric_ode=sol.success and odeerr<2e-15,
        ward_finite_difference=ward[-1]<1e-12,
        ward_refinement=ward[-1]<ward[0]/10,
        perturbative_horizon_negative=shift<0,
        local_large_mass_parameter=float(max(local_curvature))<.01)
    info=dict(model='local DeWitt O(mu^-2), first-order Einstein response only; no full black-hole RSET',
       muM=muM,NG_over_M2=NGM2,epsilon=eps,N_example=100,G_over_M2=.0001,
       delta_horizon_over_M=shift,delta_mass_at_2M_over_M=eps*dmass(2),chi_at_2M=eps*chi(2),
       max_quadrature_error=max(integral_errors),metric_ode_error=odeerr,ward_errors=ward,
       max_sqrtK_over_mu2=float(max(local_curvature)),checks=checks)
    return info,dict(bh_x=x,bh_tt=tt(x),bh_rr=rr(x),bh_tang=tang(x),
        bh_delta_mass_over_epsilon=dmass(x),bh_chi_over_epsilon=chi(x),bh_metric_ode=ode,
        bh_sqrtK_over_mu2=local_curvature)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tag',required=True)
    ap.add_argument('--skip-mode-integrals',action='store_true')
    args=ap.parse_args()
    assert all(c.isalnum() or c in '_-' for c in args.tag)
    out=ROOT/f'work/absource_numerics_{args.tag}.json'
    arr=out.with_suffix('.npz')
    assert not out.exists() and not arr.exists()
    started=time.perf_counter()
    ds,a=massive_desitter(not args.skip_mode_integrals)
    bh,b=schwarzschild();a.update(b)
    assert all(np.isfinite(v).all() for v in a.values())
    np.savez_compressed(arr,**a)
    data=dict(schema='absolute-dirac-source-numerics-v1',utc=datetime.now(timezone.utc).isoformat(),
       source_sha256=sha(__file__),inputs={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in
       [ROOT/'inputs/absource_dirac_ds2014.html',ROOT/'inputs/absource_matyjasek2000.html']},
       arrays_file=str(arr.relative_to(ROOT)).replace('\\','/'),arrays_sha256=sha(arr),
       desitter=ds,schwarzschild=bh,mode_integrals_performed=not args.skip_mode_integrals,
       all_checks_pass=all(ds['checks'].values()) and all(bh['checks'].values()),elapsed_seconds=time.perf_counter()-started)
    out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:data[k] for k in ['all_checks_pass','elapsed_seconds']},indent=2))
    print(ds['checks']);print(bh['checks']);print(out)
    return 0 if data['all_checks_pass'] else 1

if __name__=='__main__':
    raise SystemExit(main())
