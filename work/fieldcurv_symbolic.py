"""Independent 4D Dirac stress and energy-work identities on warped PG space.

Builds the actual Levi-Civita and tetrad spin connections, then contracts
the symmetric four-component Dirac tensor.  No radial stress formula is
imported from another program.  Units c=hbar=1 and angular channel kappa=1.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def christoffels(g, x):
    inverse = g.inv()
    return [[[s.simplify(sum(inverse[a, d]*(s.diff(g[d, c], x[b])
        + s.diff(g[d, b], x[c])-s.diff(g[b, c], x[d]))/2
        for d in range(4))) for c in range(4)] for b in range(4)] for a in range(4)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', default='20261018')
    args = parser.parse_args()
    output = ROOT/'work'/f'fieldcurv_symbolic_{args.tag}.json'
    assert not output.exists(), 'Use a fresh tag; existing receipts are immutable.'
    started = time.perf_counter()
    checks = []
    real_derivatives = {}
    def simplify(x):
        inverse = {v:k for k,v in real_derivatives.items()}
        return s.simplify(s.trigsimp(s.expand(x.xreplace(real_derivatives)))).xreplace(inverse)
    def check(name, value, description):
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        values = [simplify(x) for x in entries]
        checks.append(dict(name=name, passed=all(x==0 for x in values),
            residuals=[str(x) for x in values], description=description))
        if not checks[-1]['passed']:
            print(name, checks[-1]['residuals'], flush=True)

    T, r, theta, phi = s.symbols('T r theta phi', real=True)
    x = (T, r, theta, phi)
    R = s.Function('R', positive=True)(r)
    v = s.Function('v', real=True)(r)
    real_derivatives.update({s.diff(f,r,n):s.Symbol(f'{f.func}_derivative_{n}',real=True)
                            for f in (R,v) for n in (1,2)})
    mu, M, a, GN = s.symbols('mu M a GN', positive=True, real=True)
    g = s.Matrix([[-1+v*v, v, 0, 0], [v, 1, 0, 0],
        [0, 0, R*R, 0], [0, 0, 0, R*R*s.sin(theta)**2]])
    Gamma = christoffels(g, x)
    coframe = s.Matrix([[1,0,0,0], [v,1,0,0], [0,0,R,0], [0,0,0,R*s.sin(theta)]])
    frame = coframe.inv()
    eta = s.diag(1,-1,-1,-1)
    sx = s.Matrix([[0,1],[1,0]])
    sy = s.Matrix([[0,-s.I],[s.I,0]])
    sz = s.diag(1,-1)
    ident = s.eye(2)
    alpha = [-s.kronecker_product(sy,ident), s.kronecker_product(sx,sx),
             s.kronecker_product(sx,sy)]
    beta = s.kronecker_product(sz,ident)
    gamma = [beta]+[beta*aa for aa in alpha]
    check('actual_tetrad', coframe.T*(-eta)*coframe-g, 'Actual warped PG metric with signature -+++.')
    spin = []
    for index in range(4):
        omega_matrix = s.zeros(4)
        for aa in range(4):
            for bb in range(4):
                omega = sum(coframe[aa, nu]*(s.diff(frame[nu, bb], x[index])
                    +sum(Gamma[nu][index][lam]*frame[lam, bb] for lam in range(4)))
                    for nu in range(4))
                omega_matrix += eta[aa, aa]*s.simplify(omega)*gamma[aa]*gamma[bb]/4
        spin.append(s.simplify(omega_matrix))
    print('Actual spin connection constructed',flush=True)

    ar, ai, br, bi, arp, aip, brp, bip, arpp, aipp, brpp, bipp = s.symbols(
        'ar ai br bi arp aip brp bip arpp aipp brpp bipp', real=True)
    variables = (ar, ai, br, bi)
    first = (arp, aip, brp, bip)
    second = (arpp, aipp, brpp, bipp)
    psi = s.Matrix([ar+s.I*ai, br+s.I*bi])
    psir = s.Matrix([arp+s.I*aip, brp+s.I*bip])
    psirr = s.Matrix([arpp+s.I*aipp, brpp+s.I*bipp])
    A = sz-v*ident
    V = mu*sx-sy/R
    psit = -A*psir+s.diff(v,r)*psi/2-s.I*V*psi
    psirt = -A*psirr+3*s.diff(v,r)*psir/2+s.diff(v,r,2)*psi/2-s.I*V*psir-s.I*s.diff(V,r)*psi
    W = s.Matrix([[1,s.I],[1,-s.I]])/s.sqrt(2)
    angular = s.exp(s.I*phi/2)*s.Matrix([s.cos(theta/2),s.I*s.sin(theta/2)])/s.sqrt(4*s.pi)
    field = s.kronecker_product(W.conjugate().T*psi, angular)/R
    partial = [s.kronecker_product(W.conjugate().T*psit,angular)/R,
        s.kronecker_product(W.conjugate().T*psir,angular)/R-s.diff(R,r)*field/R,
        s.diff(field,theta), s.diff(field,phi)]
    coordinate_D = [partial[i]+spin[i]*field for i in range(4)]
    D = [sum((frame[i,j]*coordinate_D[i] for i in range(4)),s.zeros(4,1)) for j in range(4)]
    dirac = s.I*sum((gamma[j]*D[j] for j in range(4)),s.zeros(4,1))-mu*field
    check('actual_full_4D_Dirac_equation',dirac,
          'Actual tetrad spin connection and full angular spinor obey the supplied radial PDE.')
    # Standard +--- gamma matrices give physical covariant components for
    # the -+++ metric through the following symmetric stress.  T00 is
    # calibrated below to positive energy for a flat positive-energy wave.
    C = [s.eye(4)]+[-aa for aa in alpha]
    def stress(i,j):
        return simplify(-s.im((field.conjugate().T*C[i]*D[j]
                            +field.conjugate().T*C[j]*D[i])[0])/2)
    raw = {(i,j):stress(i,j) for i in range(4) for j in range(i,4)}
    scale = 4*s.pi*R**2
    P = simplify(s.im((psi.conjugate().T*sz*psir)[0]))
    J = simplify(s.im((psi.conjugate().T*psir)[0]))
    s1 = simplify((psi.conjugate().T*sx*psi)[0])
    s2 = simplify((psi.conjugate().T*sy*psi)[0])
    potential = mu*s1-s2/R
    energy = P+potential
    tangent = -s2/(2*R)
    for name, ij, value in [('energy',(0,0),energy),('radial_covariant_momentum',(0,1),-J),
                           ('radial_pressure',(1,1),P),('polar_pressure',(2,2),tangent),
                           ('azimuthal_pressure',(3,3),tangent)]:
        check('actual_4D_stress_'+name,scale*raw[ij]-value,
              'Projection of the actual symmetric 4D stress; the result is angle independent.')
    q = simplify((psi.conjugate().T*psi)[0])
    s3 = simplify((psi.conjugate().T*sz*psi)[0])
    s1prime = 2*simplify(s.re((psi.conjugate().T*sx*psir)[0]))
    K2 = simplify(s.im((psi.conjugate().T*sy*psir)[0]))
    angular_energy = (-3*q+s.diff(R,r)*s1+R*s1prime)/(4*R)
    angular_radial = (s3+(R*s.diff(v,r)-v*s.diff(R,r))*s1-2*R*K2)/(4*R)
    def angular_rationalize(expression):
        # Exact half-angle substitution on the standard sphere chart
        # 0<theta<pi.  It avoids opaque heuristic trigonometric cancellation.
        h=s.Symbol('positive_half_angle_tangent',positive=True)
        half=s.Symbol('half_polar_angle',real=True)
        expanded=s.expand_trig(expression.subs(theta,2*half))
        return s.factor(s.together(expanded.subs({
            s.sin(half):h/s.sqrt(1+h*h),s.cos(half):1/s.sqrt(1+h*h),
            s.tan(half):h,s.cot(half):1/h},simultaneous=True)))
    for name,ij,value in [('polar_energy',(0,2),0),('azimuthal_energy',(0,3),s.sin(theta)*angular_energy),
                         ('radial_polar',(1,2),0),('radial_azimuthal',(1,3),s.sin(theta)*angular_radial),
                         ('polar_azimuthal',(2,3),0)]:
        check('actual_4D_stress_'+name,angular_rationalize(scale*raw[ij]-value),
              'Additional actual symmetric stress projection, with exact chartwise half-angle simplification.')
    all_reduced={'00':energy,'01':-J,'02':s.S.Zero,'03':s.sin(theta)*angular_energy,
                 '11':P,'12':s.S.Zero,'13':s.sin(theta)*angular_radial,
                 '22':tangent,'23':s.S.Zero,'33':tangent}
    check('on_shell_stress_trace',-energy+P+2*tangent+mu*s1,
          'T^a_a=-mu*bar(Psi)*Psi for signature -+++, after angular integration.')
    print('Actual tetrad stress projections checked',flush=True)

    time_values = tuple(simplify(x) for p in psit for x in (s.re(p),s.im(p)))
    time_first = tuple(simplify(x) for p in psirt for x in (s.re(p),s.im(p)))
    def dt(expression):
        return sum(s.diff(expression,z)*zt for z,zt in zip(variables,time_values)) + sum(
            s.diff(expression,z)*zt for z,zt in zip(first,time_first))
    def dr(expression):
        return s.diff(expression,r)+sum(s.diff(expression,z)*zr for z,zr in zip(variables,first)) + sum(
            s.diff(expression,z)*zr for z,zr in zip(first,second))
    work = s.diff(v,r)*P+2*v*s.diff(R,r)*tangent/R
    force = s.diff(v,r)*J+2*s.diff(R,r)*tangent/R
    killing = energy-v*J
    flux = (1+v*v)*J-v*(energy+P)
    check('Euler_energy_work_balance',dt(energy)+dr(J-v*energy)-work,
          'Direct derivatives from the actual PDE, retaining spatial gradients and all connection terms.')
    check('radial_momentum_force_balance',dt(J)+dr(P-v*J)-force,
          'Direct PDE verification of the sphere-integrated radial momentum equation.')
    check('Killing_energy_continuity',dt(killing)+dr(flux),
          'Stationary Killing energy has exactly zero local source.')
    Hpsi = -s.I*A*psir+s.I*s.diff(v,r)*psi/2+V*psi
    check('Killing_energy_Hamiltonian_density',s.re((psi.conjugate().T*Hpsi)[0])-killing,
          'Real Hermitian Hamiltonian density, without replacement by mu times probability density.')
    check('Killing_energy_Hamiltonian_flux',s.re((s.I*psi.conjugate().T*A*psit)[0])-flux,
          'Local Hamiltonian flux independently matches the geometric Killing current.')
    current = simplify((psi.conjugate().T*A*psi)[0])
    check('probability_continuity',dt(q)+dr(current),'Actual PDE probability continuity remains distinct from energy conservation.')
    axial_density = q/2-dr(R*s1)/6
    axial_flux = 2*R*(v*angular_energy+angular_radial)/3
    check('actual_axial_Killing_density',-2*R*angular_energy/3-axial_density,
          'Exact sphere integral of T^T_phi; integral sin(theta)^2 dOmega=8*pi/3.')
    check('axial_Killing_flux_improvement',axial_flux-current/2-R*dt(s1)/6,
          'Actual angular momentum flux differs from half the probability flux by a total-time improvement.')
    check('axial_Killing_continuity',dt(axial_density)+dr(axial_flux),
          'Full axial stress components have a separately conserved angular-momentum current.')

    # Independent geometric reduced-divergence calculation, before the PDE
    # identities above are used.  The time/radial metric sector is enough
    # after angular integration; angular divergence terms integrate to zero.
    ee, jj, pp, tt = s.symbols('ee jj pp tt',real=True)
    frame_cov = s.Matrix([[ee,-jj],[-jj,pp]])
    cf2 = coframe[:2,:2]
    reduced_cov = cf2.T*frame_cov*cf2
    h = g[:2,:2]
    mixed = s.simplify(h.inv()*reduced_cov)
    contravariant = s.simplify(h.inv()*reduced_cov*h.inv())
    check('geometric_reduced_momentum_density',mixed[0,1]-jj,'Area-weighted T^T_r.')
    check('geometric_reduced_momentum_flux',mixed[1,1]-(pp-v*jj),'Area-weighted T^r_r.')
    check('geometric_radial_metric_force',sum(contravariant[i,j]*s.diff(h[i,j],r)/2
          for i in range(2) for j in range(2))-s.diff(v,r)*jj,
          'The metric derivative contribution in the reduced covariant radial divergence.')
    check('geometric_Killing_density',-mixed[0,0]-(ee-v*jj),'Minus area-weighted T^T_T.')
    check('geometric_Killing_flux',-mixed[1,0]-((1+v*v)*jj-v*(ee+pp)),
          'Minus area-weighted T^r_T.')
    print('Energy, force and Killing identities checked',flush=True)

    # Source requirement is taken from the actual four-dimensional Einstein
    # tensor, evaluated independently here for generic R(r),v(r).
    ricci = s.Matrix(4,4,lambda aa,bb:s.simplify(sum(
        s.diff(Gamma[c][aa][bb],x[c])-s.diff(Gamma[c][aa][c],x[bb])
        +sum(Gamma[c][c][d]*Gamma[d][aa][bb]-Gamma[c][bb][d]*Gamma[d][aa][c]
             for d in range(4)) for c in range(4))))
    inverse = g.inv()
    scalar = s.simplify(sum(inverse[i,j]*ricci[i,j] for i in range(4) for j in range(4)))
    einstein = s.simplify(ricci-scalar*g/2)
    source_null = {}
    for sign,name in ((1,'plus'),(-1,'minus')):
        k = s.Matrix([1,-v+sign,0,0])
        geometric = s.simplify((k.T*einstein*k)[0])
        check('actual_geometric_null_'+name,geometric+2*s.diff(R,r,2)*(sign-v)**2/R,
              'Actual regular-PG four-dimensional Einstein contraction.')
        field_null = energy+P-2*sign*J
        check('actual_field_null_'+name,scale*(raw[0,0]+2*sign*raw[0,1]+raw[1,1])-field_null,
              'Actual Dirac radial null contraction; it has no assumed positivity.')
        source_null[name] = str(geometric/(8*s.pi*GN)-field_null/scale)
    rpos=s.Symbol('rpos',positive=True)
    expressions={}
    for name,rad,vel in [('Schwarzschild',r,s.sqrt(2*M/r)),
                         ('black_bounce',s.sqrt(r*r+a*a),s.sqrt(2*M/s.sqrt(r*r+a*a)))]:
        sub={R:rad,v:vel,s.diff(R,r):s.diff(rad,r),s.diff(v,r):s.diff(vel,r),
             s.diff(R,r,2):s.diff(rad,r,2)}
        expressions[name]=dict(Euler_work=str(s.simplify(work.subs(sub,simultaneous=True))),
            momentum_force=str(s.simplify(force.subs(sub,simultaneous=True))),
            area_density='(P+mu*s1-kappa*s2/R)/(4*pi*R**2)',
            radial_NEC_geometry=str(s.simplify((-2*s.diff(R,r,2)*(1+v)**2/R).subs(sub,simultaneous=True))))
        if name=='Schwarzschild':
            check('actual_Schwarzschild_vacuum',einstein.subs(sub,simultaneous=True).subs(r,rpos),
                  'Every component of the actual 4D Einstein tensor vanishes for r>0.')
            check('Schwarzschild_explicit_Euler_work',work.subs(sub,simultaneous=True)
                  +vel*P/(2*r)+vel*s2/r**2,
                  'Explicit work coefficients on the singular Schwarzschild background.')
        else:
            check('regular_explicit_Euler_work',work.subs(sub,simultaneous=True)
                  +vel*r*P/(2*rad**2)+vel*r*s2/rad**3,
                  'Explicit regular-background work coefficients, finite at r=0 for a>0.')
    check('regular_throat_zero_local_work',work.subs({s.diff(v,r):0,s.diff(R,r):0}),
          'At the reflection-symmetric throat both local Euler work coefficients vanish; nearby they need not.')
    plane_energy=s.Symbol('plane_energy',positive=True)
    check('flat_plane_energy_sign_calibration',-s.im((psi.conjugate().T*(-s.I*plane_energy*psi))[0])-plane_energy*q,
          'The actual T00 prescription gives E times density for exp(-i E T), fixing the stress sign.')
    payload=dict(schema='field-curvature-actual-4D-Dirac-stress-v1',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        checks=checks, check_count=len(checks),
        scalar_residual_count=sum(len(c['residuals']) for c in checks),
        all_checks_pass=all(c['passed'] for c in checks),
        conventions=dict(metric_signature='-+++',gamma_signature='+---',
            stress='T_ab=-(1/2) Im[Psi^dagger C_a D_b Psi+Psi^dagger C_b D_a Psi], C=(I,-alpha)',
            angular_channel='kappa=1,j=m_j=1/2',
            reduction='All reduced projections equal 4*pi*R**2 times the physical local projection'),
        formulas=dict(P='Im(psi^dagger sigma3 partial_r psi)',J='Im(psi^dagger partial_r psi)',
            s_i='psi^dagger sigma_i psi',W='mu*s1-kappa*s2/R',
            Euler_energy='e=P+W',radial_momentum='-4*pi*R**2*T(n,e_r)=J',
            radial_pressure='P',tangential_pressure='pt=-kappa*s2/(2*R)',
            Euler_energy_balance='partial_T e+partial_r(J-v*e)=v_prime*P+2*v*R_prime/R*pt',
            radial_momentum_balance='partial_T J+partial_r(P-v*J)=v_prime*J+2*R_prime/R*pt',
            Killing_energy='epsilon=e-v*J',Killing_flux='F=(1+v**2)*J-v*(e+P)',
            Killing_balance='partial_T epsilon+partial_r F=0',
            field_radial_null='4*pi*R**2*T(k_plusminus,k_plusminus)=2*P+W-minusplus2*J',
            angular_energy='4*pi*R**2*T03=sin(theta)*[-3*q+R_prime*s1+R*partial_r_s1]/(4*R)',
            angular_radial='4*pi*R**2*T13=sin(theta)*[s3+(R*v_prime-v*R_prime)*s1-2*R*K2]/(4*R)',
            K2='Im(psi^dagger sigma2 partial_r psi)',
            axial_density='ell=q/2-partial_r(R*s1)/6',
            axial_flux='Fell=j/2+R*partial_T(s1)/6',
            axial_balance='partial_T ell+partial_r Fell=0; integral ell dr=1/2 if norm=1 and R*s1 vanishes at both boundaries',
            extra_null_budget=source_null),
        all_reduced_symmetric_components={key:str(value) for key,value in all_reduced.items()},
        models=expressions,
        scope='Classical bilinear test-spinor tensor, also the finite excitation difference under an explicitly specified admissible Hadamard reference; not the absolute renormalized vacuum tensor and not a self-consistent Einstein solution.',
        seconds=time.perf_counter()-started)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(output),checks=len(checks),passed=payload['all_checks_pass'],
        scalar_residuals=payload['scalar_residual_count'],seconds=payload['seconds'])),flush=True)
    assert payload['all_checks_pass'],[c['name'] for c in checks if not c['passed']]


if __name__=='__main__':main()
