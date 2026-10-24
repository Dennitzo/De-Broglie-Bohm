"""Independent Einstein-tensor audit of the massive Dirac a3 backreaction.

This checks the algebraic order-mu^-2 Schwinger--DeWitt sector, not a
complete state-selected RSET. All four Einstein components are built from
the actual four-dimensional metric before a first-order expansion.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='20261018')
    args=parser.parse_args()
    output=ROOT/'work'/f'absource_schwarzschild_audit_{args.tag}.json'
    assert not output.exists(), 'Use a fresh output tag.'
    start=time.perf_counter()
    t,r,theta,phi=s.symbols('t r theta phi',real=True)
    M,mu,G,N=s.symbols('M mu G N',positive=True)
    eps=s.symbols('eps',real=True)
    coords=(t,r,theta,phi)
    F=s.Function('F')(r)
    chi=s.Function('chi')(r)
    diag=[-s.exp(2*chi)*F,1/F,r*r,r*r*s.sin(theta)**2]
    inv=[1/d for d in diag]
    metric=s.diag(*diag)
    gamma=[[[s.simplify(inv[i]*(s.diff(metric[i,k],coords[j])+s.diff(metric[i,j],coords[k])-s.diff(metric[j,k],coords[i]))/2)
             for k in range(4)] for j in range(4)] for i in range(4)]
    ricci=s.Matrix(4,4,lambda i,j:s.simplify(sum(s.diff(gamma[k][i][j],coords[k])-s.diff(gamma[k][i][k],coords[j])
        +sum(gamma[k][i][j]*gamma[l][k][l]-gamma[l][i][k]*gamma[k][j][l] for l in range(4)) for k in range(4))))
    scalar=s.simplify(sum(inv[i]*ricci[i,i] for i in range(4)))
    einstein=s.Matrix(4,4,lambda i,j:s.simplify(inv[i]*ricci[i,j]-(scalar/2 if i==j else 0)))
    f0=1-2*M/r
    P=M*M/(20160*s.pi**2*mu**2*r**9)
    tt=P*(2384*M-1080*r)
    rr=P*(504*r-784*M)
    th=P*(3536*M-1512*r)
    source=[tt,rr,th,th]
    dm=G*N*M*M*(149*M-81*r)/(1890*s.pi*mu**2*r**6)
    dc=-11*G*N*M*M/(210*s.pi*mu**2*r**6)
    df=-2*dm/r
    replace={F:f0+eps*df,chi:eps*dc}
    for order in (1,2):
        replace[s.diff(F,r,order)]=s.diff(f0+eps*df,r,order)
        replace[s.diff(chi,r,order)]=s.diff(eps*dc,r,order)
    checks=[]

    def check(name,value):
        value=s.factor(s.trigsimp(value))
        checks.append(dict(name=name,passed=value==0,residual=str(value)))

    for i in range(4):
        for j in range(4):
            evaluated=einstein[i,j].xreplace(replace)
            check(f'background_Einstein_{i}_{j}',evaluated.subs(eps,0))
            coefficient=s.diff(evaluated,eps).subs(eps,0)
            check(f'linear_Einstein_{i}_{j}',coefficient-(8*s.pi*G*N*source[i] if i==j else 0))
    check('stress_covariant_radial_Ward',s.diff(rr,r)+s.diff(f0,r)/(2*f0)*(rr-tt)+2*(rr-th)/r)
    check('mass_radial_equation',s.diff(dm,r)+4*s.pi*G*N*r*r*tt)
    check('lapse_radial_equation',s.diff(dc,r)-4*s.pi*G*N*r*(rr-tt)/f0)
    check('mass_boundary_at_infinity',s.limit(dm,r,s.oo))
    check('lapse_boundary_at_infinity',s.limit(dc,r,s.oo))
    check('horizon_regular_difference',(rr-tt)/f0-11*M*M/(140*s.pi**2*mu**2*r**8))
    horizon=2*M
    dmH=s.factor(dm.subs(r,horizon))
    dcH=s.factor(dc.subs(r,horizon))
    drH=s.factor(-df.subs(r,horizon)/s.diff(f0,r).subs(r,horizon))
    check('horizon_displacement',drH+13*G*N/(60480*s.pi*mu**2*M**3))
    check('horizon_constraint_linear',df.subs(r,horizon)+drH*s.diff(f0,r).subs(r,horizon))
    dkappa=s.factor((dcH*s.diff(f0,r).subs(r,horizon)+s.diff(df,r).subs(r,horizon)+drH*s.diff(f0,r,2).subs(r,horizon))/2)
    check('surface_gravity_correction',dkappa+G*N/(241920*s.pi*mu**2*M**5))
    dtemp=s.factor(dkappa/(2*s.pi))
    check('temperature_correction',dtemp+G*N/(483840*s.pi**2*mu**2*M**5))
    # Dimensional and dimensionless forms: mu has inverse-length dimension.
    x,z,lam=s.symbols('x z lam',positive=True)
    substitutes={r:M*x,mu:z/M,G:lam*M*M/N}
    dm_unit=s.factor((dm/M).subs(substitutes, simultaneous=True))
    chi_unit=s.factor(dc.subs(substitutes, simultaneous=True))
    check('dimensionless_mass',dm_unit-lam*(149-81*x)/(1890*s.pi*z*z*x**6))
    check('dimensionless_lapse',chi_unit+11*lam/(210*s.pi*z*z*x**6))
    ksqrt=4*s.sqrt(3)*M/r**3
    eta=s.factor(ksqrt/mu**2)
    check('curvature_scale_at_horizon',eta.subs(r,2*M)-s.sqrt(3)/(2*mu**2*M*M))
    expressions={
        'source_P':str(P),'source_mixed':[str(s.factor(q)) for q in source],
        'source_trace':str(s.factor(tt+rr+2*th)),
        'source_rr_minus_tt_over_f0':str(s.factor((rr-tt)/f0)),
        'generic_Einstein_mixed_diagonal':[str(s.factor(einstein[i,i])) for i in range(4)],
        'delta_m':str(dm),'chi':str(dc),'delta_m_horizon':str(dmH),'chi_horizon':str(dcH),
        'delta_r_horizon':str(drH),'delta_surface_gravity':str(dkappa),'delta_temperature':str(dtemp),
        'relative_surface_gravity_correction':str(s.factor(4*M*dkappa)),
        'dimensionless_delta_m_over_M':str(dm_unit),'dimensionless_chi':str(chi_unit),
        'eta_curvature':str(eta),'eta_horizon':str(eta.subs(r,2*M)),
        'eta_unity_radius':str((4*s.sqrt(3)*M/mu**2)**s.Rational(1,3)),
    }
    sources=['inputs/absource_matyjasek2000.html','inputs/absource_matyjasek2000.txt']
    report=dict(passed=all(c['passed'] for c in checks),check_count=len(checks),checks=checks,expressions=expressions,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        primary_source_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        sympy_version=s.__version__,elapsed_seconds=time.perf_counter()-start,
        scope='First order in the local renormalized a3/mu^2 effective-action source; physical metric evaluated at eps=1 only in perturbative validity regime. No global state or complete RSET is constructed.')
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(output),passed=report['passed'],checks=len(checks),expressions=expressions),indent=2))
    assert report['passed']


if __name__=='__main__':
    main()
