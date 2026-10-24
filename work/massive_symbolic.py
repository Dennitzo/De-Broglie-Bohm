"""Exact SymPy audit of massive PG Dirac reduction and Bohm currents.

Run from the repository root:
    python work/massive_symbolic.py

The audit constructs Cartesian Pauli/Dirac matrices, applies the full 3+1
operator to the kappa=+1 spinor, and checks exact symbolic residuals. It is
not a numerical PDE solver, a Lean proof, or a quantum backreaction model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    script = Path(__file__).resolve()
    root = script.parent.parent
    output = args.output or root / "work/massive_symbolic_20261018.json"
    checks = []

    def reduce(expr):
        return sp.simplify(sp.trigsimp(sp.expand(expr)))

    def check(name, expression, explanation):
        entries = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
        residuals = [reduce(entry) for entry in entries]
        passed = all(entry == 0 for entry in residuals)
        checks.append(
            {
                "name": name,
                "passed": passed,
                "residuals": [str(entry) for entry in residuals],
                "explanation": explanation,
            }
        )

    r, M, mu = sp.symbols("r M mu", positive=True, real=True)
    p, kappa = sp.symbols("p kappa", real=True)
    theta, phi = sp.symbols("theta phi", real=True)
    ar, ai, br, bi = sp.symbols("ar ai br bi", real=True)
    dar, dai, dbr, dbi = sp.symbols("dar dai dbr dbi", real=True)
    a, b = ar + sp.I * ai, br + sp.I * bi
    da, db = dar + sp.I * dai, dbr + sp.I * dbi
    psi = sp.Matrix([a, b])
    psi_r = sp.Matrix([da, db])
    z = sp.conjugate(a) * b
    q = sp.conjugate(a) * a + sp.conjugate(b) * b
    v = sp.sqrt(2 * M / r)
    vp = sp.diff(v, r)
    identity = sp.eye(2)
    sigma = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    zero = sp.zeros(2)
    alpha = [
        zero.row_join(s).col_join(s.row_join(zero)) for s in sigma
    ]
    beta = sp.diag(1, 1, -1, -1)
    W = (sigma[0] + sigma[2]) / sp.sqrt(2)
    n = sp.Matrix(
        [sp.sin(theta) * sp.cos(phi),
         sp.sin(theta) * sp.sin(phi),
         sp.cos(theta)]
    )
    e_theta = sp.Matrix(
        [sp.cos(theta) * sp.cos(phi),
         sp.cos(theta) * sp.sin(phi),
         -sp.sin(theta)]
    )
    e_phi = sp.Matrix([-sp.sin(phi), sp.cos(phi), 0])

    def matrix_projection(matrices, vector):
        return sum((vector[i] * matrices[i] for i in range(3)),
                   sp.zeros(matrices[0].rows))

    sigma_r = matrix_projection(sigma, n)
    alpha_r = matrix_projection(alpha, n)
    alpha_theta = matrix_projection(alpha, e_theta)
    alpha_phi = matrix_projection(alpha, e_phi)

    check("pg_shift_derivative", vp + v / (2 * r),
          "The incoming PG shift is v=sqrt(2M/r), v'=-v/(2r).")
    check("full_shift_half_divergence",
          sp.diff(r**2 * v, r) / (2 * r**2) - 3 * v / (4 * r),
          "The full Cartesian spatial half-divergence equals 3v/(4r).")
    check("basis_unitarity", W.conjugate().T * W - identity,
          "The constant Hadamard basis transformation is unitary.")
    check("basis_radial_alpha", W * sigma[0] * W - sigma[2],
          "The radial alpha matrix sigma1 transforms into sigma3.")
    check("basis_mass_beta", W * sigma[2] * W - sigma[0],
          "The mass matrix sigma3 transforms into sigma1.")
    check("basis_angular_sign", W * sigma[1] * W + sigma[1],
          "The kappa sigma2/r term changes sign under this basis change.")

    # eta is a Cartesian constant spinor. Its normalization is restored
    # explicitly in the full operator check; currents below suppress the
    # common factor 1/(4*pi*r^2).
    eta = sp.Matrix([1, 0])

    def sigma_dot_L(spinor):
        # L=-i n cross (e_theta d_theta + e_phi/sin(theta) d_phi).
        angular_L = [
            -sp.I * (
                e_phi[i] * sp.diff(spinor, theta)
                - e_theta[i] * sp.diff(spinor, phi) / sp.sin(theta)
            )
            for i in range(3)
        ]
        return sum((sigma[i] * angular_L[i] for i in range(3)),
                   sp.zeros(2, 1))

    check("angular_kappa_plus_one", sigma_dot_L(eta),
          "(sigma.L+1) eta=+eta for eta=(1,0); kappa=+1.")
    eta_other = sigma_r * eta
    check("angular_partner_kappa_minus_one",
          sigma_dot_L(eta_other) + 2 * eta_other,
          "(sigma.L+1) sigma_r eta=-sigma_r eta, with actual angular derivatives.")
    check("angular_mj_plus_half",
          -sp.I * sp.diff(eta, phi) + sigma[2] * eta / 2 - eta / 2,
          "J_z eta=eta/2, so m_j=+1/2 in this convention.")
    check("pointwise_angular_norm",
          (eta_other.conjugate().T * eta_other)[0] - 1,
          "The lower angular partner has the same pointwise norm.")

    # Actual Cartesian 4x4-alpha currents of the reconstructed 4-spinor.
    u = W * psi
    Psi_current = (u[0] * eta).col_join(u[1] * eta_other)
    J = sp.Matrix([
        sp.expand_complex((Psi_current.conjugate().T * A * Psi_current)[0])
        for A in alpha
    ])
    density4 = (Psi_current.conjugate().T * Psi_current)[0]
    radial = sp.conjugate(a) * a - sp.conjugate(b) * b
    check("cartesian_density", density4 - q,
          "Direct four-spinor density, with the common 1/(4*pi*r^2) suppressed.")
    check("cartesian_radial_current", n.dot(J) - radial,
          "Direct projection of Psi^dagger alpha_i Psi onto e_r.")
    check("cartesian_polar_current", e_theta.dot(J),
          "Direct projection onto e_theta is zero in this angular sector.")
    check("cartesian_azimuthal_current",
          e_phi.dot(J) + 2 * sp.im(z) * sp.sin(theta),
          "Direct projection onto e_phi is -2 Im(a* b) sin(theta).")
    check("full_four_dimensional_causal_norm",
          -q**2 + J.dot(J)
          + 4 * (sp.re(z)**2 + sp.im(z)**2 * sp.cos(theta)**2),
          "Actual Cartesian spatial current gives the full PG tangent norm.")
    check("radial_two_dimensional_causal_norm",
          -q**2 + radial**2
          + 4 * sp.conjugate(a) * a * sp.conjugate(b) * b,
          "The radial 1+1 norm has a different strictness condition.")

    # Apply the true Cartesian spatial derivative written in spherical
    # coordinates; alpha_i remain the constant Cartesian Dirac matrices.
    # This tests the operator, not just a postulated current identity.
    U1 = sp.Function("U1")(r)
    U2 = sp.Function("U2")(r)
    u_radial = sp.Matrix([U1, U2])
    eta_normalized = eta / sp.sqrt(4 * sp.pi)
    eta_other_normalized = sigma_r * eta_normalized
    Psi = (U1 * eta_normalized).col_join(U2 * eta_other_normalized) / r
    H4Psi = (
        -sp.I * (
            alpha_r * sp.diff(Psi, r)
            + alpha_theta * sp.diff(Psi, theta) / r
            + alpha_phi * sp.diff(Psi, phi) / (r * sp.sin(theta))
        )
        + mu * beta * Psi
        + sp.I * v * (sp.diff(Psi, r) + 3 * Psi / (4 * r))
    )
    Hu = (
        -sp.I * (sigma[0] - v * identity) * sp.diff(u_radial, r)
        + sp.I * vp * u_radial / 2
        + mu * sigma[2] * u_radial
        + sigma[1] * u_radial / r
    )
    reconstructed_Hu = (
        Hu[0] * eta_normalized
    ).col_join(Hu[1] * eta_other_normalized) / r
    check("actual_3plus1_angular_operator_separation",
          H4Psi - reconstructed_Hu,
          "All four components of full PG H4 applied to Psi agree with the kappa=+1 radial operator.")
    check("radial_half_density_connection",
          r * v * (sp.diff(U1 / r, r) + 3 * U1 / (4 * r**2))
          - v * sp.diff(U1, r) - vp * U1 / 2,
          "Densitization by 1/r converts the full 3v/(4r) term into v'/2.")

    B = mu * sigma[0] - kappa * sigma[1] / r
    potential_dt = -sp.I * B * psi
    component_dt = [
        sp.conjugate(psi[i]) * potential_dt[i]
        + sp.conjugate(potential_dt[i]) * psi[i] for i in range(2)
    ]
    transfer = 2 * mu * sp.im(z) + 2 * kappa * sp.re(z) / r
    check("mass_angular_transfer_upper", component_dt[0] - transfer,
          "The actual potential matrix transfers +2mu Im(z)+2kappa Re(z)/r into the upper density.")
    check("mass_angular_transfer_lower", component_dt[1] + transfer,
          "The lower component receives the exact opposite transfer.")
    check("local_potential_hermiticity", B - B.conjugate().T,
          "Both mass and angular coupling are Hermitian local terms.")
    check("local_potential_square",
          B * B - (mu**2 + kappa**2 / r**2) * identity,
          "This identity gives the exact sine/cosine potential propagator.")
    A = sigma[2] - v * identity
    psi_dt = -A * psi_r + vp * psi / 2 - sp.I * B * psi
    q_dt = (psi_dt.conjugate().T * psi + psi.conjugate().T * psi_dt)[0]
    j_r = (
        psi_r.conjugate().T * A * psi
        + psi.conjugate().T * A * psi_r
        - vp * psi.conjugate().T * psi
    )[0]
    check("full_density_continuity", q_dt + j_r,
          "Direct symbolic time/space derivatives of q and j cancel with the connection term.")
    check("horizon_flux",
          (psi.conjugate().T * (sigma[2] - identity) * psi)[0]
          + 2 * sp.conjugate(b) * b,
          "At v=1, j=-2|b|^2, independently of mass and kappa.")

    # Derivatives of the exact lower and upper interior arrival-time bounds.
    g = r - 2 * sp.sqrt(2 * M * r) + 4 * M * sp.log(1 + sp.sqrt(r / (2 * M)))
    G = -r - 2 * sp.sqrt(2 * M * r) - 4 * M * sp.log(1 - sp.sqrt(r / (2 * M)))
    check("lower_arrival_primitive", sp.diff(g, r) - 1 / (1 + v),
          "g'=1/(1+v), the fastest future radial null infall time.")
    check("upper_arrival_primitive", sp.diff(G, r) - 1 / (v - 1),
          "G'=1/(v-1), used only on the real interval 0<r<2M.")

    # Geodesic Hamiltonian is checked against the actual inverse PG metric.
    hgeo = sp.sqrt(p**2 + mu**2) - v * p
    metric = sp.Matrix([[-(1-v**2), v], [v, 1]])
    momentum = sp.Matrix([-hgeo, p])
    check("geodesic_mass_shell",
          (momentum.T * metric.inv() * momentum)[0] + mu**2,
          "p_mu g^{mu nu}p_nu=-mu^2 with positive frame energy H+vp.")
    rdot = sp.diff(hgeo, p)
    pdot = -sp.diff(hgeo, r)
    check("geodesic_radial_equation",
          rdot - p / sp.sqrt(p**2 + mu**2) + v,
          "Hamilton's radial equation in PG time.")
    check("geodesic_momentum_equation", pdot - vp * p,
          "Hamilton's momentum equation has sign +v' p.")
    check("geodesic_energy_conservation",
          sp.diff(hgeo, r) * rdot + sp.diff(hgeo, p) * pdot,
          "The autonomous geodesic Hamiltonian is conserved.")

    # Concrete main-case local polarization; it is not a Killing-frequency
    # projection of the compact spatial packet.
    E0 = sp.sqrt(5)
    chi = sp.Matrix([sp.sqrt((E0-2)/(2*E0)), sp.sqrt((E0+2)/(2*E0))])
    check("initial_local_spinor_norm", (chi.T * chi)[0] - 1,
          "The p0=-2, mu=1 local polarization is normalized.")
    check("initial_local_radial_velocity", (chi.T * sigma[2] * chi)[0] + 2/E0,
          "Its initial local radial velocity equals the comparison geodesic's p0/E0.")
    check("initial_azimuthal_velocity", sp.im(chi[0] * chi[1]),
          "Real initial components give zero initial azimuthal velocity.")

    report = {
        "schema": "massive-pg-symbolic-audit-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(script.relative_to(root)).replace("\\", "/"),
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__,
        "metric_signature": "-+++",
        "units": "c=hbar=1; M is the geometric Schwarzschild mass",
        "operator": "-i(sigma3-vI)d_r + i v'/2 I + mu sigma1 - kappa sigma2/r",
        "angular_convention": "(sigma.L+1)eta=kappa eta, eta=(1,0)/sqrt(4pi), kappa=+1, mj=+1/2",
        "domain": "r>0, M>0; upper arrival primitive only 0<r<2M; polar coordinates away from their axis",
        "check_count": len(checks),
        "passed_count": sum(item["passed"] for item in checks),
        "all_passed": all(item["passed"] for item in checks),
        "checks": checks,
        "scope": [
            "Exact symbolic residual identities for the specified differential operator and angular sector.",
            "No numerical PDE convergence claim and no Lean kernel proof are made by this script.",
            "No backreaction, renormalized source, positive Killing-frequency projection, or central extension is computed.",
            "The causal norm identities do not imply strict timelikeness merely from mu>0.",
        ],
        "primary_sources": [
            {"archive": "inputs/massive_source_lasenby_2005.html",
             "url": "https://arxiv.org/html/gr-qc/0209090v2",
             "equations": ["18-23", "26-27", "37-40"]},
            {"archive": "inputs/massive_source_doran_2005.html",
             "url": "https://arxiv.org/html/gr-qc/0503019v1",
             "equations": ["9-15", "23-24"]},
            {"archive": "inputs/massive_source_dolan_dempsey_2015.html",
             "url": "https://arxiv.org/html/1504.03190v2",
             "equations": ["123"]},
        ],
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "check_count": report["check_count"],
        "passed_count": report["passed_count"],
        "all_passed": report["all_passed"],
        "failed": [item["name"] for item in checks if not item["passed"]],
    }, ensure_ascii=False))
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
