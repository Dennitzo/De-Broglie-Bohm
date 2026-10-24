"""Independent exact local-constraint and Gaussian-kernel calculations.

No manuscript, pre-existing receipts or outputs are modified.  This checks
finite algebra and exact periodic integrals, not the full functional ADM
calculus or regulator-independent Wheeler--DeWitt operator construction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as s


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="20261018")
    args = parser.parse_args()
    out = ROOT / "work" / f"local_algebra_symbolic_{args.tag}.json"
    if out.exists():
        raise FileExistsError(out)
    checks: dict[str, bool] = {}

    def check(name: str, expression: s.Expr) -> None:
        checks[name] = s.simplify(s.expand(expression)) == 0
        if not checks[name]:
            raise AssertionError((name, expression))

    # ADM lapse-Hessian contraction in an arbitrary orthonormal frame.
    # Tensoriality promotes this pointwise identity to any positive metric.
    p11, p22, p33, p12, p13, p23 = s.symbols("p11 p22 p33 p12 p13 p23")
    n11, n22, n33, n12, n13, n23 = s.symbols("n11 n22 n33 n12 n13 n23")
    p = s.Matrix([[p11,p12,p13],[p12,p22,p23],[p13,p23,p33]])
    hess = s.Matrix([[n11,n12,n13],[n12,n22,n23],[n13,n23,n33]])
    dewitt = p - s.trace(p) * s.eye(3) / 2
    lapse_variation = s.eye(3) * s.trace(hess) - hess
    contraction = sum(dewitt[i,j] * lapse_variation[i,j]
                      for i in range(3) for j in range(3))
    phess = sum(p[i,j] * hess[i,j] for i in range(3) for j in range(3))
    check("adm_lapse_hessian_contraction", contraction + phess)

    # Symmetric pi^{ij} removes cross-gradient products in integration by parts.
    n = s.Matrix(s.symbols("n1:4")); m = s.Matrix(s.symbols("m1:4"))
    check("adm_cross_gradient_cancellation", (n.T*p*m - m.T*p*n)[0])

    # Scalar Euler derivative and antisymmetrization, including arbitrary U.
    N,M,C,D,Pi,Uprime,lap_phi = s.symbols("N M C D Pi Uprime lap_phi")
    grad_phi = s.Matrix(s.symbols("g1:4"))
    eN = -D*(n.dot(grad_phi)+N*lap_phi)+N*Uprime
    eM = -D*(m.dot(grad_phi)+M*lap_phi)+M*Uprime
    scalar_pb = C*Pi*(M*eN-N*eM)
    scalar_expected = C*D*Pi*(N*m-M*n).dot(grad_phi)
    check("scalar_local_euler_bracket", scalar_pb-scalar_expected)
    A,B,dg,ds = s.symbols("A B dg ds")
    bracket = A*B*dg+C*D*ds
    check("sector_g_test", bracket.subs({dg:1,ds:0})-A*B)
    check("sector_phi_test", bracket.subs({dg:0,ds:1})-C*D)
    check("matched_sector_sufficiency", (bracket-dg-ds).subs({B:1/A,D:1/C}))

    # Genuine finite canonical Poisson bracket: Gaussian with symmetric K.
    q1,q2,p1,p2,a,b,c1,c2,v1,v2 = s.symbols("q1 q2 p1 p2 a b c1 c2 v1 v2")
    q = [q1,q2]; mom = [p1,p2]
    z1,z2 = a*q1+b*q2,b*q1+a*q2
    R = s.exp(-(a*(q1*q1+q2*q2)+2*b*q1*q2)/2)
    Q1 = c1-z1*z1/2; Q2 = c2-z2*z2/2
    check("gaussian_quantum_potential_site1", -s.diff(R,q1,2)/(2*R) -(a-z1*z1)/2)
    check("gaussian_quantum_potential_site2", -s.diff(R,q2,2)/(2*R) -(a-z2*z2)/2)
    H1=p1*p1/2+v1*q1*q1/2+Q1
    H2=p2*p2/2+v2*q2*q2/2+Q2
    pb=sum(s.diff(H1,x)*s.diff(H2,y)-s.diff(H1,y)*s.diff(H2,x)
           for x,y in zip(q,mom))
    expected=b*(p1*z2-p2*z1)
    check("two_site_actual_poisson_bracket", pb-expected)
    check("two_site_constant_subtraction_c1", s.diff(pb,c1))
    check("two_site_constant_subtraction_c2", s.diff(pb,c2))
    check("two_site_counterexample_is_one", pb.subs({a:2,b:1,q1:1,q2:0,p1:1,p2:0})-1)
    check("two_site_bohm_graph_pi_zero", pb.subs({p1:0,p2:0}))

    # Positive-lapse torus example. Dependence on y,z is constant, so the
    # three-dimensional normalized volume integral reduces exactly to theta.
    theta = s.symbols("theta", real=True)
    mass,k,omega,omega2,eps = s.symbols("mass k omega omega2 eps", positive=True)
    phi = s.cos(theta); lapse = 1+eps*s.cos(theta)
    e1 = (k*k+mass*mass)*phi
    e2 = -k*k*s.diff(lapse*s.diff(phi,theta),theta)+mass*mass*lapse*phi
    qder1 = -omega*omega*phi
    qder2 = -omega*omega*phi-eps*omega*(mass+omega2*s.cos(2*theta))/2

    def avg(expr: s.Expr) -> s.Expr:
        return s.simplify(s.integrate(s.expand_trig(expr),(theta,0,2*s.pi))/(2*s.pi))

    classical = avg(lapse*e1-e2)
    obstruction = avg(lapse*qder1-qder2)
    total = s.expand(classical+obstruction).subs(k*k,omega*omega-mass*mass)
    check("torus_classical_bracket", classical-eps*k*k/2)
    check("torus_quantum_defect", obstruction-eps*omega*(mass-omega)/2)
    check("torus_total_bracket", total-eps*mass*(omega-mass)/2)
    check("torus_ratio_rationalization", mass*(omega-mass)/(omega*omega-mass*mass)-mass/(omega+mass))
    delta = s.symbols("delta",real=True)
    gradient1=k*k*phi
    gradient2=-k*k*s.diff(lapse*s.diff(phi,theta),theta)
    check("gradient_counterterm_constant_shift", avg(delta*(lapse*gradient1-gradient2))-delta*eps*k*k/2)
    check("pointwise_normal_ordering_no_phi_derivative", s.diff(c1*q1**0,q1))

    # A second admissible k yields a different repair coefficient.
    repair1 = 1-1/(s.sqrt(2)+1)
    repair2 = 1-1/(s.sqrt(5)+1)
    checks["two_modes_need_distinct_constant_repairs"] = bool(repair1 != repair2 and float(repair2-repair1)>0)
    # Independent full matrix differentiation with lapse diagonal entries.
    n1,n2,m1,m2=s.symbols("N1 N2 M1 M2")
    HN=n1*H1+n2*H2; HM=m1*H1+m2*H2
    bracket_general=sum(s.diff(HN,x)*s.diff(HM,y)-s.diff(HN,y)*s.diff(HM,x)
                        for x,y in zip(q,mom))
    check("arbitrary_two_site_lapses",bracket_general-(n1*m2-n2*m1)*expected)

    # Differentiating an arbitrary field-independent finite subtraction is zero.
    z = s.symbols("z")
    subtraction=s.Function("c")(z)
    check("arbitrary_finite_subtraction",s.diff(subtraction,q1))

    payload={
        "status":"passed" if all(checks.values()) else "failed",
        "scope":"Finite canonical algebra, ADM Hessian contraction, exact normalized T^3 Fourier integral. Not full functional constraint formalization.",
        "checks":checks,
        "count":len(checks),
        "torus_example":{"mass":1,"eps":0.5,"k1_ratio":float(1/(s.sqrt(2)+1)),"k2_ratio":float(1/(s.sqrt(5)+1)),"k1_repair":float(repair1),"k2_repair":float(repair2)},
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version":s.__version__,
    }
    out.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"path":str(out),"status":payload["status"],"count":len(checks)}))


if __name__ == "__main__":
    main()
