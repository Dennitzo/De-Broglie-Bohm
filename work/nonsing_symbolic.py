"""Actual 4D geometry, tetrad-Dirac reduction and currents for a black bounce.

No Schwarzschild radial operator is assumed: the warped PG spin connection
is constructed from the metric and tetrad before the 1/R rescaling.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def connection(metric,coordinates):
    inv=metric.inv()
    return [[[s.simplify(sum(inv[a,d]*(s.diff(metric[d,c],coordinates[b])
        +s.diff(metric[d,b],coordinates[c])-s.diff(metric[b,c],coordinates[d]))/2
        for d in range(4))) for c in range(4)] for b in range(4)] for a in range(4)]


def ricci_from_connection(G,coordinates):
    return s.Matrix(4,4,lambda a,b:s.simplify(sum(
        s.diff(G[c][a][b],coordinates[c])-s.diff(G[c][a][c],coordinates[b])
        +sum(G[c][c][d]*G[d][a][b]-G[c][b][d]*G[d][a][c]
             for d in range(4)) for c in range(4))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='20261018')
    args=parser.parse_args()
    output=ROOT/'work'/f'nonsing_symbolic_{args.tag}.json'
    assert not output.exists(),'Use a fresh tag.'
    start=time.perf_counter()
    checks=[]
    def simp(x):return s.simplify(s.trigsimp(s.expand(x)))
    def check(name,value,description):
        entries=list(value) if isinstance(value,s.MatrixBase) else [value]
        residuals=[simp(x) for x in entries]
        checks.append(dict(name=name,passed=all(x==0 for x in residuals),
                           residuals=[str(x) for x in residuals],description=description))
    T,r,theta,phi=s.symbols('T r theta phi',real=True)
    coordinates=(T,r,theta,phi)
    R=s.Function('R',positive=True)(r)
    F=s.Function('F',real=True)(r)
    v=s.Function('v',real=True)(r)
    M,a,mu=s.symbols('M a mu',positive=True,real=True)
    static=s.diag(-F,1/F,R**2,R**2*s.sin(theta)**2)
    Gs=connection(static,coordinates)
    riemann={}
    for aa in range(4):
        for bb in range(4):
            for cc in range(4):
                for dd in range(4):
                    value=s.diff(Gs[aa][dd][bb],coordinates[cc])-s.diff(Gs[aa][cc][bb],coordinates[dd])
                    value+=sum(Gs[aa][cc][ee]*Gs[ee][dd][bb]-Gs[aa][dd][ee]*Gs[ee][cc][bb] for ee in range(4))
                    riemann[aa,bb,cc,dd]=s.simplify(value)
    inv=static.inv()
    K=s.simplify(sum(inv[aa,aa]*inv[bb,bb]*inv[cc,cc]*inv[dd,dd]
          *(static[aa,aa]*value)**2 for (aa,bb,cc,dd),value in riemann.items()))
    expectedK=(s.diff(F,r,2)**2+2*(s.diff(F,r)*s.diff(R,r)/R)**2
       +8*(F*s.diff(R,r,2)+s.diff(F,r)*s.diff(R,r)/2)**2/R**2
       +4*(1-F*s.diff(R,r)**2)**2/R**4)
    check('actual_4d_Kretschmann',K-expectedK,'All 256 actual static-chart Riemann components are contracted with the inverse metric.')
    RR=s.sqrt(r*r+a*a);FF=1-2*M/RR
    substitutions={F:FF,R:RR,s.diff(F,r):s.diff(FF,r),s.diff(F,r,2):s.diff(FF,r,2),
                   s.diff(R,r):s.diff(RR,r),s.diff(R,r,2):s.diff(RR,r,2)}
    A=M*(2*r*r-a*a)/RR**5
    B=-M*r*r/RR**5
    C=(M*(2*a*a-r*r)-a*a*RR)/RR**5
    D=(2*M*r*r+a*a*RR)/RR**5
    Ksv=4*(A*A+2*B*B+2*C*C+D*D)
    check('black_bounce_K_components',expectedK.subs(substitutions,simultaneous=True)-Ksv,'Explicit curvature formula with a>0; no division by F remains.')
    check('throat_curvature',Ksv.subs(r,0)-(36*M*M/a**6-32*M/a**5+12/a**4),'Finite throat curvature.')
    check('Schwarzschild_limit',Ksv.subs(a,0).subs(r,s.Symbol('rp',positive=True))-48*M*M/s.Symbol('rp',positive=True)**6,'Schwarzschild invariant recovered away from the center.')
    print('Geometry contraction checked',flush=True)

    pg=s.Matrix([[-(1-v*v),v,0,0],[v,1,0,0],[0,0,R*R,0],[0,0,0,R*R*s.sin(theta)**2]])
    G=connection(pg,coordinates)
    ricci=ricci_from_connection(G,coordinates)
    scalar=s.simplify(sum(pg.inv()[aa,bb]*ricci[aa,bb] for aa in range(4) for bb in range(4)))
    einstein=s.simplify(ricci-scalar*pg/2)
    check('pg_determinant',pg.det()+R**4*s.sin(theta)**2,'The PG metric determinant is independent of the horizons.')
    check('pg_global_time_covector',pg.inv()[0,0]+1,'The chosen time covector has unit timelike norm.')
    for sign,name in ((-1,'minus'),(1,'plus')):
        k=s.Matrix([1,-v+sign,0,0])
        check('null_vector_'+name,(k.T*pg*k)[0],'Actual regular PG null vector.')
        contraction=s.simplify((k.T*einstein*k)[0])
        check('radial_null_Einstein_'+name,contraction+2*s.diff(R,r,2)*(-v+sign)**2/R,
              'Actual four-dimensional Einstein tensor, not a static observer expression.')
        check('black_bounce_NEC_'+name,contraction.subs({R:RR,s.diff(R,r):s.diff(RR,r),s.diff(R,r,2):s.diff(RR,r,2)},simultaneous=True)
              +2*a*a*(-v+sign)**2/RR**4,'Strictly negative for k_minus, all real r and a>0, v>0.')
    print('Regular PG Einstein/NEC checked',flush=True)

    # Spherical tetrad and constant Hamilton matrices.  Standard (+---)
    # Clifford matrices are used; multiplying the spacetime metric by -1
    # does not change its Levi-Civita connection.
    p1=s.Matrix([[0,1],[1,0]]);p2=s.Matrix([[0,-s.I],[s.I,0]]);p3=s.diag(1,-1);I=s.eye(2)
    alpha=[-s.kronecker_product(p2,I),s.kronecker_product(p1,p1),s.kronecker_product(p1,p2)]
    beta=s.kronecker_product(p3,I)
    gamma=[beta]+[beta*aa for aa in alpha]
    eta=s.diag(1,-1,-1,-1)
    for aa in range(4):
        for bb in range(aa,4):
            check(f'Clifford_{aa}{bb}',gamma[aa]*gamma[bb]+gamma[bb]*gamma[aa]-2*eta[aa,bb]*s.eye(4),'Constant Clifford matrices.')
    coframe=s.Matrix([[1,0,0,0],[v,1,0,0],[0,0,R,0],[0,0,0,R*s.sin(theta)]])
    frame=coframe.inv()
    check('actual_tetrad_metric',coframe.T*(-eta)*coframe-pg,'Actual orthonormal tetrad of the warped PG metric.')
    spin=[]
    for index in range(4):
        matrix=s.zeros(4)
        for aa in range(4):
            for bb in range(4):
                omega=sum(coframe[aa,nu]*(s.diff(frame[nu,bb],coordinates[index])
                    +sum(G[nu][index][lam]*frame[lam,bb] for lam in range(4))) for nu in range(4))
                omega=s.simplify(eta[aa,aa]*omega)
                matrix+=omega*gamma[aa]*gamma[bb]/4
        spin.append(s.simplify(matrix))
    Hconnection=-s.I*(spin[0]+(alpha[0]-v*s.eye(4))*spin[1]
                  +alpha[1]*spin[2]/R+alpha[2]*spin[3]/(R*s.sin(theta)))
    expected_connection=(-s.I*alpha[0]*s.diff(R,r)/R
                         -s.I*alpha[1]*s.cot(theta)/(2*R)
                         +s.I*(s.diff(v,r)/2+v*s.diff(R,r)/R)*s.eye(4))
    check('actual_tetrad_spin_connection',Hconnection-expected_connection,
          'Spin connection constructed from all tetrad derivatives and actual PG Christoffels.')
    eta_plus=s.exp(s.I*phi/2)*s.Matrix([s.cos(theta/2),s.I*s.sin(theta/2)])/s.sqrt(4*s.pi)
    DS=-s.I*(p1*(s.diff(eta_plus,theta)+s.cot(theta)*eta_plus/2)+p2*s.diff(eta_plus,phi)/s.sin(theta))
    check('actual_angular_eigenfunction',DS-eta_plus,'Intrinsic sphere Dirac eigenvalue +1.')
    u1,u2=s.Function('u1')(r),s.Function('u2')(r)
    u=s.Matrix([u1,u2]);Psi=s.kronecker_product(u,eta_plus)/R
    Hfull=(-s.I*(alpha[0]-v*s.eye(4))*s.diff(Psi,r)
          -s.I*alpha[1]*s.diff(Psi,theta)/R
          -s.I*alpha[2]*s.diff(Psi,phi)/(R*s.sin(theta))
          +mu*beta*Psi+Hconnection*Psi)
    Hu=-s.I*(-p2-v*I)*s.diff(u,r)+s.I*s.diff(v,r)*u/2+mu*p3*u+p1*u/R
    check('full_4d_radial_rescaling',Hfull-s.kronecker_product(Hu,eta_plus)/R,
          'Actual four-component Dirac equation: R-weighted radial amplitude removes every R-prime term.')
    # Unitary map: tau3 -> sigma1, tau1 -> -sigma2, -tau2 -> sigma3.
    W=s.Matrix([[1,s.I],[1,-s.I]])/s.sqrt(2)
    check('radial_basis_unitarity',W*W.conjugate().T-I,'Unitary conversion to the publication radial basis.')
    check('radial_basis_alpha',W*(-p2)*W.conjugate().T-p3,'Diagonal radial transport.')
    check('radial_basis_mass',W*p3*W.conjugate().T-p1,'Mass potential mu sigma1.')
    check('radial_basis_angular',W*p1*W.conjugate().T+p2,'Angular potential -kappa sigma2/R.')
    ar,ai,br,bi=s.symbols('ar ai br bi',real=True)
    chi=s.Matrix([ar+s.I*ai,br+s.I*bi]);uu=W.conjugate().T*chi
    full=s.kronecker_product(uu,eta_plus)
    q=(chi.conjugate().T*chi)[0];z=s.conjugate(chi[0])*chi[1]
    currents=[simp(4*s.pi*(full.conjugate().T*aa*full)[0]) for aa in alpha]
    check('actual_4d_density',4*s.pi*(full.conjugate().T*full)[0]-q,'Suppressed common R^-2 factor is restored in the physical density.')
    check('actual_4d_radial_current',currents[0]-(abs(chi[0])**2-abs(chi[1])**2),'Radial frame current.')
    check('actual_4d_polar_current',currents[1],'Polar current vanishes.')
    check('actual_4d_azimuthal_current',currents[2]+2*s.im(z)*s.sin(theta),'Azimuthal frame current gives phi_dot=-2Im(a*b)/(R q).')
    check('full_current_causal_norm',-q*q+sum(j*j for j in currents)
          +4*(s.re(z)**2+s.im(z)**2*s.cos(theta)**2),'Complete four-dimensional causal norm, including angular velocity.')
    total_scalar=sum(len(row['residuals']) for row in checks)
    payload=dict(schema='black-bounce-actual-tensor-dirac-v1',
       source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       checks=checks,check_count=len(checks),scalar_residual_count=total_scalar,
       all_checks_pass=all(row['passed'] for row in checks),
       formulas=dict(K=str(Ksv),K_at_throat=str(Ksv.subs(r,0)),
                     uniform_K_bound='4*(6*M**2/a**6+3*(2*M/a**3+1/a**2)**2)',
                     ingoing_NEC='G(k_minus,k_minus)=-2*a**2*(1+v)**2/R**4',
                     radial_H='-i*(sigma3-v I)*d_r + i*v_prime/2 I + mu*sigma1-kappa/R*sigma2'),
       scope='Prescribed regular black-bounce metric; Dirac test field; supporting exotic Einstein source is not derived from this test field.',
       seconds=time.perf_counter()-start)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(output),checks=len(checks),scalar_residuals=total_scalar,
                         passed=payload['all_checks_pass'],seconds=payload['seconds'])),flush=True)
    assert payload['all_checks_pass'],[row for row in checks if not row['passed']]


if __name__=='__main__':main()
