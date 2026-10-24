"""Free entangled Gaussian: state integrability, guidance holonomy and equivariance.

Dimensionless hbar=m1=m2=1. Finite configuration-space diagnostic, not a
relativistic spacetime or a coupled gravitational solution.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
A0 = np.array([[2., 1.], [1., 2.]])
SIGMA0 = np.linalg.inv(A0)/2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kernel(times, a=A0):
    return np.linalg.inv(np.linalg.inv(a)+1j*np.diag(times))


def covariance(times, a=A0):
    return np.linalg.inv(kernel(times, a).real)/2


def clock(s, eta):
    return np.array([s+eta*np.sin(np.pi*s)/np.pi,
                     s-eta*np.sin(np.pi*s)/np.pi])


def rates(s, eta):
    return np.array([1+eta*np.cos(np.pi*s), 1-eta*np.cos(np.pi*s)])


def flow(eta, a=A0, tolerance=2e-12, samples=301):
    def rhs(s, flat):
        return (-np.diag(rates(s, eta))@kernel(clock(s, eta), a).imag
                @flat.reshape(2, 2)).ravel()
    solution = solve_ivp(rhs, (0., 1.), np.eye(2).ravel(), method='DOP853',
                         rtol=tolerance, atol=tolerance/30,
                         t_eval=np.linspace(0, 1, samples))
    assert solution.success
    return solution.t, solution.y.T.reshape(-1, 2, 2)


def leg(times, coordinate, target, initial, a=A0):
    start = times.copy()
    def rhs(t, flat):
        ts = start.copy(); ts[coordinate] = t
        velocity = np.zeros((2, 2)); velocity[coordinate] = -kernel(ts, a).imag[coordinate]
        return (velocity@flat.reshape(2, 2)).ravel()
    solution = solve_ivp(rhs, (start[coordinate], target), initial.ravel(),
                         method='DOP853', rtol=8e-13, atol=2e-14)
    assert solution.success
    start[coordinate] = target
    return start, solution.y[:, -1].reshape(2, 2)


def rectangle(order, endpoint, a=A0):
    ts = np.zeros(2); f = np.eye(2)
    for coordinate in order:
        ts, f = leg(ts, coordinate, endpoint[coordinate], f, a)
    return f


def symbolic_checks():
    a, b, c, t, u, x, y = sp.symbols('a b c t u x y', real=True)
    delta = a*c-b*b
    den = 1-delta*t*u+sp.I*(a*t+c*u)
    w = sp.Matrix([[a+sp.I*delta*u, b], [b, c+sp.I*delta*t]])/den
    a0 = sp.Matrix([[a, b], [b, c]])
    sigma = (a0.inv()+sp.diag(t,u)*a0*sp.diag(t,u))/2
    checks = {}
    def zero(name, expression):
        checks[name] = bool(all(sp.simplify(item)==0 for item in expression)) if isinstance(expression, sp.MatrixBase) else bool(sp.simplify(expression)==0)
    zero('kernel_inverse', w*(a0.inv()+sp.I*sp.diag(t,u))-sp.eye(2))
    for i, time in enumerate([t, u]):
        e = sp.zeros(2); e[i,i] = 1
        zero(f'riccati_time_{i+1}', sp.I*w.diff(time)-w*e*w)
        zero(f'covariance_transport_{i+1}', sigma.diff(time)+e*sp.im(w)*sigma+sigma*sp.im(w)*e)
        # With partial_t log(normalization)=-i W_ii/2 these are both
        # coefficients of i partial_t Psi = -partial_qi^2 Psi/2.
        zero(f'schrodinger_constant_{i+1}', sp.I*(-sp.I*w[i,i]/2)-w[i,i]/2)
        zero(f'schrodinger_quadratic_{i+1}', -sp.I*w.diff(time)/2+w*e*w/2)
    zero('covariance_inverse_real_kernel', 2*sp.re(w)*sigma-sp.eye(2))
    zero('determinant_identity', sp.det(sigma)-(sp.re(den)**2+sp.im(den)**2)/(4*delta))
    q = sp.Matrix([x,y]); v = -sp.im(w)*q
    curvature = sp.Matrix([-v[0].diff(u),v[1].diff(t)]).subs({t:0,u:0})
    k = sp.Matrix([[-b*b,-b*c],[a*b,b*b]])
    zero('initial_curvature',curvature-k*q)
    zero('curvature_preserves_initial_covariance',k*a0.inv()+a0.inv()*k.T)
    zero('curvature_square',k*k+b*b*delta*sp.eye(2))
    zero('curvature_determinant',k.det()-b*b*delta)
    # Mixed derivative commutation and Hamiltonian commutation hold on this
    # Schwartz Gaussian; direct polynomial coefficient check on initial R.
    r=sp.exp(-(a*x*x+2*b*x*y+c*y*y)/2)
    zero('free_hamiltonians_commute_on_initial_gaussian',sp.diff(r,x,2,y,2)-sp.diff(r,y,2,x,2))
    assert all(checks.values()), checks
    return checks


def main(tag):
    assert tag.replace('_','').replace('-','').isalnum()
    out=ROOT/f'work/found_numerics_{tag}.json'
    arrayfile=out.with_suffix('.npz')
    assert not out.exists() and not arrayfile.exists(), 'Use a fresh tag.'
    checks=symbolic_checks()
    arrays={}; records=[]; q0=np.array([.7,-.2])
    for name,eta in [('synchronous',0.),('clock_plus',.6),('clock_minus',-.6)]:
        s, fs=flow(eta)
        endpoint=fs[-1]; sigmas=np.array([covariance(clock(v,eta)) for v in s])
        pushed=fs@SIGMA0@fs.transpose(0,2,1)
        error=float(np.max(abs(pushed-sigmas)))
        determinant_error=float(np.max(abs(np.linalg.det(fs)**2*np.linalg.det(SIGMA0)-np.linalg.det(sigmas))))
        _,tight=flow(eta,tolerance=2e-13)
        convergence=float(np.max(abs(tight[-1]-endpoint)))
        arrays.update({name+'_flow':fs,name+'_times':np.array([clock(v,eta) for v in s]),
                       name+'_path':fs@q0,name+'_covariance':sigmas})
        records.append(dict(name=name,eta=eta,endpoint=(endpoint@q0).tolist(),flow=endpoint.tolist(),
                            covariance_error=error,determinant_error=determinant_error,
                            tolerance_comparison=convergence))
        checks[name+'_covariance_transport']=error<2e-10
        checks[name+'_determinant_transport']=determinant_error<2e-10
        checks[name+'_tolerance_refinement']=convergence<2e-10
        checks[name+'_positive_rates']=bool(np.min([rates(v,eta) for v in s])>=.399999999999)
    arrays['s']=s;arrays['q0']=q0;arrays['a0']=A0;arrays['sigma0']=SIGMA0
    endpoint_gap=float(np.linalg.norm(arrays['clock_plus_path'][-1]-arrays['clock_minus_path'][-1]))
    checks['strictly_positive_clock_endpoint_difference']=endpoint_gap>.1
    # Relative endpoint maps preserve the INITIAL Gaussian. Whitening yields
    # a rotation rather than identity, not a change of endpoint Born density.
    evals,evecs=np.linalg.eigh(SIGMA0)
    sqrt0=(evecs*np.sqrt(evals))@evecs.T
    invsqrt0=(evecs/np.sqrt(evals))@evecs.T
    relative=np.linalg.solve(arrays['clock_minus_flow'][-1],arrays['clock_plus_flow'][-1])
    rotation=invsqrt0@relative@sqrt0
    orth_error=float(np.max(abs(rotation@rotation.T-np.eye(2))))
    angle=float(np.arctan2(rotation[1,0],rotation[0,0]))
    checks['relative_holonomy_is_rotation']=orth_error<3e-10 and abs(np.linalg.det(rotation)-1)<3e-10
    arrays['relative_map']=relative;arrays['whitened_rotation']=rotation
    f12=rectangle([0,1],np.ones(2));f21=rectangle([1,0],np.ones(2))
    arrays['rectangle12_flow']=f12;arrays['rectangle21_flow']=f21
    checks['rectangular_endpoint_different']=np.linalg.norm((f12-f21)@q0)>.1
    checks['rectangular_same_covariance']=np.max(abs(f12@SIGMA0@f12.T-f21@SIGMA0@f21.T))<3e-11
    eps=np.array([.08,.04,.02,.01,.005,.0025])
    k=np.array([[-1.,-2.],[2.,1.]])
    scaled=[]
    for h in eps:
        fxy=rectangle([0,1],np.array([h,h]));fyx=rectangle([1,0],np.array([h,h]))
        scaled.append((fxy-fyx)/(h*h))
    scaled=np.array(scaled);errors=np.max(abs(scaled-k),axis=(1,2))
    arrays['loop_epsilon']=eps;arrays['loop_curvature']=scaled;arrays['loop_error']=errors;arrays['curvature_matrix']=k
    checks['small_rectangle_curvature']=errors[-1]<1e-4 and np.all(errors[:-1]/errors[1:]>3.8)
    _,sepp=flow(.6,a=np.diag([2.,3.]));_,sepm=flow(-.6,a=np.diag([2.,3.]))
    sep_error=float(np.max(abs(sepp[-1]-sepm[-1])))
    checks['unentangled_path_independence']=sep_error<1e-10
    lambdas,rotation0=np.linalg.eigh(A0)
    ts=np.linspace(-10,10,801)
    scales=np.sqrt(1+ts[:,None]**2*lambdas[None,:]**2)
    sync=np.array([(rotation0*row)@rotation0.T for row in scales])
    arrays['global_time']=ts;arrays['global_sync_flow']=sync;arrays['global_sync_path']=sync@q0
    arrays['global_sync_scales']=scales
    sync_end=(rotation0*np.sqrt(1+lambdas**2))@rotation0.T
    checks['synchronous_exact_solution']=np.max(abs(sync_end-arrays['synchronous_flow'][-1]))<1e-10
    grid=np.linspace(-4,4,161); xx,yy=np.meshgrid(grid,grid)
    points=np.stack([xx,yy],axis=-1); w=kernel([1,1])
    density=np.sqrt(np.linalg.det(w.real))/np.pi*np.exp(-np.einsum('...i,ij,...j->...',points,w.real,points))
    arrays['density_x']=xx;arrays['density_y']=yy;arrays['endpoint_density']=density
    assert all(checks.values()), checks
    assert all(np.isfinite(v).all() for v in arrays.values())
    np.savez_compressed(arrayfile,**arrays)
    result=dict(schema='foundations-gaussian-clock-v1',time_utc=datetime.now(timezone.utc).isoformat(),
                source_sha256=sha(__file__),arrays_file=str(arrayfile.relative_to(ROOT)).replace('\\','/'),
                arrays_sha256=sha(arrayfile),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,
                checks={key:bool(value) for key,value in checks.items()},all_passed=True,
                clock_results=records,positive_clock_endpoint_distance=endpoint_gap,
                relative_rotation_angle_radians=angle,relative_rotation_orthogonality_error=orth_error,
                separable_endpoint_error=sep_error,small_loop_errors=errors.tolist(),
                rectangle12=f12.tolist(),rectangle21=f21.tolist(),
                scope='Finite free two-coordinate Gaussian. Exact analytic state evolution, numerical configuration transport; no gravity or spacetime foliation construction.')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({key:result[key] for key in ['all_passed','positive_clock_endpoint_distance','relative_rotation_angle_radians','clock_results','small_loop_errors']},indent=2))
    print('Checks:',len(checks))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    main(parser.parse_args().tag)
