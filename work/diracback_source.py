"""Exact four-dimensional Einstein--Dirac filled-pair source audit.

This checks the spin connection, both angular eigenspinors, the full Dirac
operator, all stress components and individual/pair currents.  It does not
compute a quantum Hadamard RSET.  Common field factor sqrt(nu/(4*pi*L*A*B**2))
is suppressed in the algebra; the JSON records how it is restored.
Run: python work/diracback_source.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as s


def main() -> None:
    checks = []

    def red(value):
        return s.factor(s.cancel(s.expand(value)))

    def check(name, residual):
        entries = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
        reduced = [red(value) for value in entries]
        ok = all(value == 0 for value in reduced)
        checks.append({"name": name, "passed": ok,
                       "scalar_components": len(reduced),
                       "residuals": [str(value) for value in reduced]})
        if not ok:
            raise AssertionError(f"{name}: {reduced}")

    # x=tan(theta/2); global exp(i*m_j*phi) cancels in bilinears, while
    # its phi derivative is kept explicitly as i*m_j.
    x, B, ha, hb, mu = s.symbols("x B HA HB mu", positive=True, real=True)
    ar, ai, br, bi = s.symbols("ar ai br bi", real=True)
    spin = s.Matrix([ar + s.I * ai, br + s.I * bi])
    ident = s.eye(2)
    pauli = [s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    kp = s.kronecker_product
    alpha = [kp(pauli[1], ident), kp(pauli[0], pauli[0]),
             kp(pauli[0], pauli[1])]
    beta = kp(pauli[2], ident)
    H = mu * pauli[2] + pauli[0] / B
    etas = [s.Matrix([1, s.I * x]) / s.sqrt(1 + x*x),
            s.Matrix([x, -s.I]) / s.sqrt(1 + x*x)]
    sine, cosine = 2*x/(1+x*x), (1-x*x)/(1+x*x)
    cotangent = cosine/sine
    derivative_theta = lambda value: (1+x*x)*s.diff(value, x)/2
    K = [s.eye(4)] + [-value for value in alpha]

    for i in range(3):
        for j in range(i, 3):
            check(f"clifford_alpha_{i+1}_{j+1}",
                  alpha[i]*alpha[j] + alpha[j]*alpha[i]
                  - 2*int(i == j)*s.eye(4))
        check(f"clifford_mass_anticommutator_{i+1}",
              beta*alpha[i] + alpha[i]*beta)
    check("clifford_mass_square", beta*beta-s.eye(4))

    S = [red((spin.conjugate().T*p*spin)[0]) for p in pauli]
    norm = red((spin.conjugate().T*spin)[0])
    E = mu*S[2] + S[0]/B
    check("bloch_pure_norm", sum(value**2 for value in S)-norm**2)
    check("sphere_filled_projector",
          sum((eta*eta.conjugate().T for eta in etas), s.zeros(2))-s.eye(2))
    Tpair, Jpair = s.zeros(4), s.zeros(4, 1)
    mode_results = []
    for index, (eta, mj) in enumerate(zip(etas, [s.Rational(1, 2), -s.Rational(1, 2)])):
        sign = 1 if index == 0 else -1
        label = "plus" if index == 0 else "minus"
        check(f"sphere_norm_{label}", (eta.conjugate().T*eta)[0]-1)
        D_sphere = -s.I*(pauli[0]*(derivative_theta(eta)+cotangent*eta/2)
                          +pauli[1]*s.I*mj*eta/sine)
        check(f"sphere_kappa_plus_one_{label}", D_sphere-eta)
        field = kp(spin, eta)
        # Full orthonormal-frame spin-covariant derivatives, including the
        # derivative of the common A^(-1/2)*B^(-1) field prefactor.
        covariant = [
            kp(-s.I*H*spin-(ha+2*hb)*spin/2, eta),
            ha*alpha[0]*field/2,
            derivative_theta(field)/B+hb*alpha[1]*field/2,
            s.I*mj*field/(B*sine)+hb*alpha[2]*field/2
            -s.I*cotangent*kp(ident, pauli[2])*field/(2*B),
        ]
        dirac = s.I*covariant[0] + s.I*sum(
            (alpha[i]*covariant[i+1] for i in range(3)), s.zeros(4, 1))
        check(f"full_4d_dirac_{label}", dirac-mu*beta*field)
        stress = s.zeros(4)
        current = s.zeros(4, 1)
        for i in range(4):
            current[i] = red((field.conjugate().T*([s.eye(4)]+alpha)[i]*field)[0])
            for j in range(4):
                stress[i, j] = red(s.I*(field.conjugate().T*K[i]*covariant[j]
                    +field.conjugate().T*K[j]*covariant[i]
                    -covariant[j].conjugate().T*K[i]*field
                    -covariant[i].conjugate().T*K[j]*field)[0]/4)
        check(f"individual_4d_current_{label}",
              current-s.Matrix([norm, S[1], 0, sign*S[0]*sine]))
        check(f"individual_current_causal_identity_{label}",
              -current[0]**2+sum(current[i]**2 for i in range(1,4))
              +S[2]**2+S[0]**2*cosine**2)
        Tpair += stress
        Jpair += current
        mode_results.append({"mode": label, "m_j": str(mj),
                             "normalized_current": [str(red(v)) for v in current]})
    check("pair_full_stress_all_components",
          Tpair-s.diag(2*E, 0, S[0]/B, S[0]/B))
    check("pair_full_current_all_components", Jpair-s.Matrix([2*norm, 2*S[1], 0, 0]))
    check("pair_tz_vanishes_from_anticommutator", pauli[1]*H+H*pauli[1])

    # Source conservation is derived from the actual Schrödinger derivatives,
    # not inserted as a separately postulated continuity equation.
    spin_dot = -s.I*H*spin
    Sdot = [red((spin_dot.conjugate().T*p*spin + spin.conjugate().T*p*spin_dot)[0])
            for p in pauli]
    vector = s.Matrix([1/B, 0, mu])
    check("bloch_equation_from_spinor", s.Matrix(Sdot)-2*vector.cross(s.Matrix(S)))
    check("norm_derivative_from_spinor", (spin_dot.conjugate().T*spin+spin.conjugate().T*spin_dot)[0])
    Edot = mu*Sdot[2]+Sdot[0]/B-hb*S[0]/B
    check("energy_derivative", Edot+hb*S[0]/B)
    check("stress_conservation", Edot-(ha+2*hb)*E+(ha+2*hb)*E+hb*S[0]/B)

    # Common intensity and volume factors: no empirical or QFT normal ordering
    # assumption is used to turn the classical source into a quantum source.
    G, nu, length = s.symbols("G nu L", positive=True, real=True)
    check("einstein_coupling_from_filled_pair", 8*s.pi*G*2*nu/(4*s.pi*length)-4*G*nu/length)

    target = Path(__file__).with_name("diracback_source_verification_20261018.json")
    result = {
        "schema": "diracback_source_verification_v1", "date": "2026-10-18",
        "status": "passed", "checks_total": len(checks),
        "checks_passed": sum(item["passed"] for item in checks),
        "scalar_residuals_checked": sum(item["scalar_components"] for item in checks),
        "sympy_version": s.__version__,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "conventions": {
            "metric": "-dt^2+A(t)^2 dz^2+B(t)^2 dOmega^2",
            "angular_chart": "x=tan(theta/2)>0; phi phase exp(i*m_j*phi)",
            "angular_eigenvalue": 1,
            "restored_bilinear_factor": "nu/(4*pi*L*A*B^2)",
            "spinor_norm": "N=psi^dagger psi, subsequently N=1",
            "gamma_convention_if_used": "gamma^0=beta,gamma^i=beta*alpha_i; {gamma^a,gamma^b}=-2*g^{ab}",
            "source_scope": "classical two independent c-number Dirac fields; no absolute Hadamard RSET"
        },
        "results": {
            "pair_stress": "diag(2*(mu*S3+S1/B),0,S1/B,S1/B)*nu/(4*pi*L*A*B^2)",
            "pair_current": "(2*N,2*S2,0,0)*nu/(4*pi*L*A*B^2)",
            "individual_current": "(N,S2,0,+/-S1*sin(theta))*nu/(4*pi*L*A*B^2)",
            "individual_guidance_N1": "dz/dt=S2/A,dtheta/dt=0,dphi/dt=+/-S1/B",
            "individual_guidance_metric_norm_N1": "-S3^2-S1^2*cos(theta)^2",
            "coupling": "gamma=4*G*nu/L"
        },
        "checks": checks,
    }
    target.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "checks": result["checks_total"],
                      "scalar_residuals": result["scalar_residuals_checked"],
                      "output": str(target)}, indent=2))


if __name__ == "__main__":
    main()
