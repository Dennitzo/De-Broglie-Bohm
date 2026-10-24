"""Independent 4D curvature and conformal-Dirac trace obstruction.

The actual Christoffel and Riemann tensors are built from the diagonal
areal-radius metric; no invariant formula or predecessor code is imported.
This is a necessary-condition calculation, not an absolute massive RSET.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', default='20261018')
    args = parser.parse_args()
    output = ROOT/'work'/f'absource_symbolic_{args.tag}.json'
    assert not output.exists(), 'Use a fresh tag; previous receipts are immutable.'
    started = time.perf_counter()
    T, R, theta, phi = s.symbols('T R theta phi', real=True)
    M, a, GN = s.symbols('M a G', positive=True)
    b = s.symbols('b', real=True)
    coordinates = (T, R, theta, phi)
    f = 1-2*M/R
    h = f*(1-a*a/R**2)
    diag = [-f, 1/h, R**2, R**2*s.sin(theta)**2]
    inv = [1/x for x in diag]
    metric = s.diag(*diag)
    Gamma = [[[
        s.factor(inv[i]*(s.diff(metric[i,k], coordinates[j])
                    +s.diff(metric[i,j], coordinates[k])
                    -s.diff(metric[j,k], coordinates[i]))/2)
        for k in range(4)] for j in range(4)] for i in range(4)]
    Riemann = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    value = s.diff(Gamma[i][j][l],coordinates[k])-s.diff(Gamma[i][j][k],coordinates[l])
                    value += sum(Gamma[i][m][k]*Gamma[m][j][l]-Gamma[i][m][l]*Gamma[m][j][k] for m in range(4))
                    value = s.factor(s.trigsimp(value))
                    if value != 0:
                        Riemann[(i,j,k,l)] = value
    Ricci = s.Matrix(4,4,lambda i,j:s.factor(sum(Riemann.get((k,i,k,j),0) for k in range(4))))
    scal = s.factor(sum(inv[i]*Ricci[i,i] for i in range(4)))
    riccisq = s.factor(sum(inv[i]*inv[j]*Ricci[i,j]**2 for i in range(4) for j in range(4)))
    kretsch = s.factor(sum(diag[i]*inv[j]*inv[k]*inv[l]*value**2 for (i,j,k,l),value in Riemann.items()))
    weylsq = s.factor(kretsch-2*riccisq+scal**2/3)
    euler = s.factor(kretsch-4*riccisq+scal**2)
    # Direct scalar covariant Laplacian, not a pre-entered radial formula.
    box = s.factor(sum(inv[i]*(s.diff(scal,coordinates[i],2)-sum(Gamma[k][i][i]*s.diff(scal,coordinates[k]) for k in range(4))) for i in range(4)))
    anomaly_numerator = s.factor(18*weylsq-11*euler)
    anomaly = anomaly_numerator/(5760*s.pi**2)
    trace_residual = s.factor(scal+GN*(anomaly_numerator+b*box)/(720*s.pi))
    checks = []

    def check(name, expr):
        value=s.factor(s.trigsimp(expr))
        checks.append({'name':name,'passed':value==0,'residual':str(value)})

    for i in range(4):
        for j in range(4):
            check(f'ricci_symmetric_{i}_{j}',Ricci[i,j]-Ricci[j,i])
    check('ricci_scalar_known',scal-2*a*a*(3*M-R)/R**5)
    check('schwarzschild_ricci',scal.subs(a,0))
    check('schwarzschild_ricci_square',riccisq.subs(a,0))
    check('schwarzschild_kretschmann',kretsch.subs(a,0)-48*M*M/R**6)
    check('schwarzschild_boxR',box.subs(a,0))
    check('schwarzschild_dirac_anomaly',anomaly.subs(a,0)-7*M*M/(120*s.pi**2*R**6))
    check('regular_R4_scalar_limit',s.limit(R**4*scal,R,s.oo)+2*a*a)
    check('regular_R4_anomaly_limit',s.limit(R**4*anomaly,R,s.oo))
    check('regular_R4_boxR_limit',s.limit(R**4*box,R,s.oo))
    check('regular_R4_trace_residual_limit',s.limit(R**4*trace_residual,R,s.oo)+2*a*a)
    check('flat_space_all_invariants',sum(x.subs({M:0,a:0})**2 for x in (scal,riccisq,kretsch,weylsq,euler,box)))
    # A constant-curvature sign/normalization check independent of this metric:
    # R=12H^2, Ricci^2=36H^4, K=24H^4 => C^2=0,E4=24H^4.
    H=s.symbols('H',real=True)
    cc_trace=(18*(24-2*36+144/s.Integer(3))-11*(24-4*36+144))*H**4/(5760*s.pi**2)
    check('constant_curvature_dirac_trace',cc_trace+11*H**4/(240*s.pi**2))
    check('DHP2009_anomaly_convention',(anomaly_numerator-12*box)/(5760*s.pi**2)
          +(scal**2/s.Integer(1152)+box/s.Integer(480)-riccisq/s.Integer(720)-7*kretsch/s.Integer(5760))/s.pi**2)
    for i,expected in enumerate([-M*a*a/R**5,a*a*(3*M-2*R)/R**5,2*M*a*a/R**5,2*M*a*a/R**5]):
        check(f'ricci_component_{i}',inv[i]*Ricci[i,i]-expected)
    u=s.symbols('u',real=True)
    scalar_poly=-2*a*a*u**4+6*M*a*a*u**5
    anomaly_poly=336*M*M*u**6+224*M*a*a*u**7+(96*a**4-1008*M*M*a*a)*u**8-424*M*a**4*u**9+888*M*M*a**4*u**10
    box_poly=-24*a*a*u**6+184*M*a*a*u**7+(32*a**4-300*M*M*a*a)*u**8-230*M*a**4*u**9+360*M*M*a**4*u**10
    check('inverse_radius_scalar_polynomial',scalar_poly-scal.subs(R,1/u))
    check('inverse_radius_anomaly_polynomial',anomaly_poly-anomaly_numerator.subs(R,1/u))
    check('inverse_radius_box_polynomial',box_poly-box.subs(R,1/u))
    expressions={
        'metric':str(metric),'nonzero_riemann_components':len(Riemann),
        'Ricci_scalar':str(scal),'Ricci_squared':str(riccisq),
        'Kretschmann':str(kretsch),'Weyl_squared':str(weylsq),
        'Euler_density':str(euler),'Box_Ricci_scalar':str(box),
        'Dirac_anomaly_numerator_5760pi2':str(anomaly_numerator),
        'Einstein_trace_residual':str(trace_residual),
        'Ricci_mixed_diagonal':[str(s.factor(inv[i]*Ricci[i,i])) for i in range(4)],
        'R6_anomaly_limit':str(s.limit(R**6*anomaly,R,s.oo)),
        'R6_BoxR_limit':str(s.limit(R**6*box,R,s.oo)),
        'R4_trace_residual_limit':str(s.limit(R**4*trace_residual,R,s.oo)),
        'inverse_radius_polynomials':{'Ricci':str(scalar_poly),'anomaly_numerator':str(anomaly_poly),'BoxR':str(box_poly)},
        'throat_invariants':{name:str(s.factor(value.subs(R,a))) for name,value in [('Ricci',scal),('K',kretsch),('C2',weylsq),('E4',euler),('BoxR',box)]},
    }
    report={
        'passed':all(c['passed'] for c in checks),'check_count':len(checks),'checks':checks,
        'expressions':expressions,
        'conventions':{'signature':'-+++','Riemann':'R^i_jkl = d_k Gamma^i_jl - d_l Gamma^i_jk + Gamma^i_mk Gamma^m_jl - Gamma^i_ml Gamma^m_jk',
        'field':'one massless four-component Dirac field; c=1/20,a_E=11/360',
        'anomaly':'(18 C^2 - 11 E4 + b BoxR)/(5760 pi^2)',
        'b':'arbitrary fixed finite real coefficient after absorbing gravitational quadratic trace terms; DHP2009 base convention b=-12',
        'trace_equation':'R + G/(720 pi) (18C^2 -11E4 + b BoxR) =0; Lambda=0; no other matter',
        'scope':'necessary trace obstruction only; does not determine traceless state-dependent stress or massive condensate'},
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version':s.__version__,'elapsed_seconds':time.perf_counter()-started,
    }
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(output),'passed':report['passed'],'checks':len(checks),'expressions':expressions},indent=2))
    assert report['passed']


if __name__=='__main__':
    main()
