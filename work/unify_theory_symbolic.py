"""Independent exact checks of the shared WDW/Einstein quantum-source model."""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
a, f, pa, pf = s.symbols('alpha phi p_alpha p_phi', real=True)
b, lam, sig = s.symbols('beta lambda sigma', positive=True)
Q = s.Function('Q')(a, f)
E = (-pa**2 + pf**2)/2 + Q
Hq = b*s.exp(-3*a)*E
variables = (a, f, pa, pf)
flow = (s.diff(Hq, pa), s.diff(Hq, pf), -s.diff(Hq, a), -s.diff(Hq, f))
def dt(expr):
    return sum(s.diff(expr, x)*v for x, v in zip(variables, flow))

checks = []
def check(name, expr):
    value = s.simplify(s.expand(expr))
    if value != 0:
        value = s.factor(value)
    checks.append({'name': name, 'exact_zero': value == 0, 'residual': str(value)})
    print(name, value == 0, flush=True)
    if value != 0:
        raise AssertionError((name, value))

H = flow[0]
fd = flow[1]
rm = fd**2/(2*b)
rq = b*s.exp(-6*a)*Q
pq = b*s.exp(-6*a)*(Q-s.diff(Q, a)/3)
exchange = b*s.exp(-6*a)*fd*s.diff(Q, f)
check('canonical_full_Hamiltonian_constant', dt(Hq))
check('rescaled_constraint_transport', dt(E)-3*H*E)
check('scalar_equation_from_Hamiltonian', dt(fd)+3*H*fd+b**2*s.exp(-6*a)*s.diff(Q,f))
check('matter_exchange', dt(rm)+6*H*rm+exchange)
check('quantum_exchange', dt(rq)+3*H*(rq+pq)-exchange)
check('total_source_conserved', dt(rm+rq)+3*H*(2*rm+rq+pq))
check('friedmann_constraint_exact_equivalence', H**2-2*b*(rm+rq)+2*b**2*s.exp(-6*a)*E)
check('acceleration_on_constraint', dt(H)+3*b*(2*rm+rq+pq)-3*b**2*s.exp(-6*a)*E)
check('quantum_pressure_from_volume_variation', -s.diff(b*s.exp(-3*a)*Q,a)/(3*s.exp(3*a))-pq)

# Hamilton-Jacobi graph preservation for the general state (constant reduced metric).
S = s.Function('S')(a, f)
Qs = (s.diff(S,a)**2-s.diff(S,f)**2)/2
sub = {Q: Qs}
on_graph = {pa:s.diff(S,a),pf:s.diff(S,f)}
check('Hamilton_Jacobi_graph_alpha', (flow[2]-dt(s.diff(S,a))).subs(sub).doit().subs(on_graph))
check('Hamilton_Jacobi_graph_phi', (flow[3]-dt(s.diff(S,f))).subs(sub).doit().subs(on_graph))

# The packet identities are checked independently in the separate audit.
Delta=s.symbols('Delta',positive=True)
Hs=-1/(3*Delta)
Hsd=-1/(3*Delta**2)
check('singular_control_Ricci',6*(Hsd+2*Hs**2)+2/(3*Delta**2))
check('singular_control_Kretschmann',12*((Hsd+Hs**2)**2+Hs**4)-20/(27*Delta**4))

report = {'all_passed':all(x['exact_zero'] for x in checks),'exact_check_count':len(checks),
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'sympy_version':s.__version__,'checks':checks,
          'scope':'Independent finite-dimensional identities, no continuum quantum-gravity proof.'}
out=ROOT/'work/unify_theory_symbolic_20261018.json'
out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
print(json.dumps({'all_passed':report['all_passed'],'exact_check_count':len(checks),'output':str(out)}))
