"""Independent exact identities behind the foundations research notes.

This verifies finite matrix/derivative identities, not the full functional
Noether theorem or a renormalized quantum gravity theory.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
checks: list[dict] = []


def check(name: str, expr) -> None:
    terms = list(expr) if isinstance(expr, sp.MatrixBase) else [expr]
    residual = [sp.simplify(sp.together(e)) for e in terms]
    ok = all(e == 0 for e in residual)
    checks.append({"name": name, "passed": ok, "residual": [str(e) for e in residual]})
    if not ok:
        raise AssertionError((name, residual))


a, b, t1, t2, x1, x2 = sp.symbols("a b t1 t2 x1 x2", real=True)
A0 = sp.Matrix([[a, b], [b, a]])
D = sp.diag(t1, t2)
W = A0.inv() + sp.I * D
A = W.inv()
E1, E2 = sp.diag(1, 0), sp.diag(0, 1)
check("gaussian_matrix_inverse", W * A - sp.eye(2))
check("partial_time1_riccati", A.diff(t1) + sp.I * A * E1 * A)
check("partial_time2_riccati", A.diff(t2) + sp.I * A * E2 * A)
det = (sp.eye(2) + sp.I * D * A0).det()
check("normalization_time1", -det.diff(t1) / (2 * det) + sp.I * A[0, 0] / 2)
check("normalization_time2", -det.diff(t2) / (2 * det) + sp.I * A[1, 1] / 2)

x = sp.Matrix([x1, x2])
R0 = sp.exp(-(x.T * A0 * x)[0] / 2)
Q1 = -R0.diff(x1, 2) / (2 * R0)
Q2 = -R0.diff(x2, 2) / (2 * R0)
check("partial_quantum_potential1", Q1 - (a - (a * x1 + b * x2) ** 2) / 2)
check("partial_quantum_potential2", Q2 - (a - (b * x1 + a * x2) ** 2) / 2)
curvature = sp.Matrix([Q2.diff(x1), -Q1.diff(x2)])
check("guidance_curvature_at_initial_time", curvature - sp.Matrix([-b**2*x1-a*b*x2, a*b*x1+b**2*x2]))
check("curvature_concrete_value", curvature.subs({a: 2, b: 1, x1: 1, x2: 0}) - sp.Matrix([-1, 2]))
check("factorized_curvature_zero", curvature.subs(b, 0))

# Euler derivative of N |grad log N|^2 in a flat spatial chart.
y = sp.symbols("y0:3", real=True)
N = sp.Function("N")(*y)
gradN = [N.diff(yi) for yi in y]
normgradN = sum(gi**2 for gi in gradN)
L_N = normgradN / N
EL_N = sp.diff(L_N, N) - sum(sp.diff(sp.diff(L_N, gi), yi) for gi, yi in zip(gradN, y))
a_sq = normgradN / N**2
div_a = sum(sp.diff(gi / N, yi) for gi, yi in zip(gradN, y))
check("lapse_acceleration_euler_derivative", EL_N + a_sq + 2 * div_a)

# The timelike normalized clock gradient varies by a spatial projector.
q = sp.Matrix(sp.symbols("q0:4", real=True))
eta = sp.diag(-1, 1, 1, 1)
X = -(q.T * eta * q)[0]
u = -q / sp.sqrt(X)
projector = sp.eye(4) + u * (eta * u).T
check("clock_unit_norm", (u.T * eta * u)[0] + 1)
check("clock_gradient_variation", u.jacobian(q) + projector / sp.sqrt(X))
check("clock_projector_annihilates_normal", projector * u)

# Explicit canonical scalar-clock Ward identity in flat spacetime.
z = sp.symbols("z0:4", real=True)
T = sp.Function("T")(*z)
gradT = sp.Matrix([T.diff(zi) for zi in z])
grad_up = eta * gradT
grad_sq = (gradT.T * grad_up)[0]
stress_mixed = grad_up * gradT.T - sp.eye(4) * grad_sq / 2
boxT = sum(eta[i, i] * T.diff(z[i], 2) for i in range(4))
for nu in range(4):
    divT = sum(stress_mixed[mu, nu].diff(z[mu]) for mu in range(4))
    check(f"canonical_clock_ward_component_{nu}", divT - boxT * gradT[nu])

sources = [Path(__file__), ROOT / "work/found_theory_notes.md", ROOT / "work/found_theory_smearing.md"]
report = {
    "status": "passed",
    "method": "independent SymPy exact matrix and derivative identities",
    "sympy": sp.__version__,
    "check_count": len(checks),
    "checks": checks,
    "bound_files": [{"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
    "scope": "Finite identities only; full variational geometry and continuum arguments are analytic proofs in the notes.",
}
output = ROOT / "work/found_theory_symbolic_20261018.json"
output.write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps({"status": "passed", "checks": len(checks), "report": str(output)}))
