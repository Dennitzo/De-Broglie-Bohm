import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "work" / "bohmian_flrw_model.json"
MARKDOWN_PATH = ROOT / "work" / "Bohmian_FLRW_Results.md"


def hamilton_constraint(a, phi, pa, p, *, V=0.0, k=0.0, Lambda=0.0):
    """Reduced FLRW Hamiltonian constraint in units kappa = 1."""
    if a <= 0:
        raise ValueError("scale factor must be positive")
    return (
        -pa * pa / (24.0 * a)
        + a * p * p / 2.0
        + a**3 * V
        - k * a
        + Lambda * a**3 / 3.0
    )


def classical_rhs(state, *, V=0.0, k=0.0, Lambda=0.0):
    a, phi, pa, p = state
    da = -pa / (12.0 * a)
    dphi = a * p
    dpa = -(pa * pa / (24.0 * a * a) + p * p / 2.0 + 3.0 * a * a * V - k + Lambda * a * a)
    dp = 0.0
    return (da, dphi, dpa, dp)


def rk4_classical(state, h, *, V=0.0, k=0.0, Lambda=0.0):
    k1 = classical_rhs(state, V=V, k=k, Lambda=Lambda)
    s2 = tuple(y + 0.5 * h * d for y, d in zip(state, k1))
    k2 = classical_rhs(s2, V=V, k=k, Lambda=Lambda)
    s3 = tuple(y + 0.5 * h * d for y, d in zip(state, k2))
    k3 = classical_rhs(s3, V=V, k=k, Lambda=Lambda)
    s4 = tuple(y + h * d for y, d in zip(state, k3))
    k4 = classical_rhs(s4, V=V, k=k, Lambda=Lambda)
    return tuple(
        y + h * (d1 + 2.0 * d2 + 2.0 * d3 + d4) / 6.0
        for y, d1, d2, d3, d4 in zip(state, k1, k2, k3, k4)
    )


def integrate_classical(state0, T, h, *, V=0.0, k=0.0, Lambda=0.0):
    steps = max(1, int(round(T / h)))
    dt = T / steps
    state = state0
    states = [state]
    for _ in range(steps):
        state = rk4_classical(state, dt, V=V, k=k, Lambda=Lambda)
        states.append(state)
    return states


def guidance_rhs(state, b, C):
    a, _ = state
    return (-b / 12.0, a * C)


def rk4_guidance(state, h, b, C):
    k1 = guidance_rhs(state, b, C)
    s2 = tuple(y + 0.5 * h * d for y, d in zip(state, k1))
    k2 = guidance_rhs(s2, b, C)
    s3 = tuple(y + 0.5 * h * d for y, d in zip(state, k2))
    k3 = guidance_rhs(s3, b, C)
    s4 = tuple(y + h * d for y, d in zip(state, k3))
    k4 = guidance_rhs(s4, b, C)
    return tuple(
        y + h * (d1 + 2.0 * d2 + 2.0 * d3 + d4) / 6.0
        for y, d1, d2, d3, d4 in zip(state, k1, k2, k3, k4)
    )


def integrate_guidance(state0, T, h, b, C):
    steps = max(1, int(round(T / h)))
    dt = T / steps
    state = state0
    states = [state]
    for _ in range(steps):
        state = rk4_guidance(state, dt, b, C)
        states.append(state)
    return states


def gaussian_quantum_potential(a, phi, hbar, a_center, phi_center, sigma_a, sigma_phi):
    """Q = -hbar^2 (d_a^2 + d_phi^2)R/(2R), R = exp(-F)."""
    fa = (a - a_center) / (2.0 * sigma_a * sigma_a)
    fphi = (phi - phi_center) / (2.0 * sigma_phi * sigma_phi)
    laplacian_over_R = (
        fa * fa
        + fphi * fphi
        - 1.0 / (2.0 * sigma_a * sigma_a)
        - 1.0 / (2.0 * sigma_phi * sigma_phi)
    )
    return -0.5 * hbar * hbar * laplacian_over_R


def quantum_constraint(a, phi, b, C, hbar, model):
    q = gaussian_quantum_potential(
        a, phi, hbar, model["a_center"], model["phi_center"], model["sigma_a"], model["sigma_phi"]
    )
    return hamilton_constraint(a, phi, b * a, C), q


def rate(e0, e1):
    if e0 <= 0.0 or e1 <= 0.0:
        return float("nan")
    return math.log(e0 / e1) / math.log(2.0)


def max_classical_constraint_drift(states, *, V=0.0, k=0.0, Lambda=0.0):
    h_values = [
        hamilton_constraint(a, phi, pa, p, V=V, k=k, Lambda=Lambda)
        for a, phi, pa, p in states
    ]
    return max(abs(h - h_values[0]) for h in h_values)


def run():
    checks = []
    data = {"model": {}, "results": {}, "checks": checks}

    def add_check(name, passed, details):
        checks.append({"name": name, "passed": bool(passed), "details": str(details)})
        print(("PASS | " if passed else "FAIL | ") + name + " | " + str(details))

    sqrt12 = math.sqrt(12.0)

    # 1. Exact massless scalar solution, k = Lambda = 0.
    state0 = (1.0, 0.2, -sqrt12, -1.0)
    T = 0.4
    states = integrate_classical(state0, T, 1.0e-3)
    drift_massless = max_classical_constraint_drift(states)
    a_num, phi_num, pa_num, p_num = states[-1]
    a_exact = 1.0 + T / (2.0 * math.sqrt(3.0))
    phi_exact = 0.2 - 1.0 * T - T * T / (4.0 * math.sqrt(3.0))
    pa_exact = -sqrt12 - T
    final_err_massless = max(
        abs(a_num - a_exact), abs(phi_num - phi_exact), abs(pa_num - pa_exact), abs(p_num + 1.0)
    )
    add_check(
        "Klassisch: masseloses Skalarfeld mit exakter Loesung",
        drift_massless < 1.0e-12 and final_err_massless < 1.0e-12,
        f"Constraint-Drift={drift_massless:.3e}; Endpunktfehler={final_err_massless:.3e}",
    )

    # 2. de-Sitter-artiger Test mit konstantem Potential V = 1/2.
    V = 0.5
    state0 = (1.0, 0.0, -sqrt12, 0.0)
    T = 0.5
    h_list = (0.05, 0.025, 0.0125)
    errors = []
    drifts = []
    for h in h_list:
        traj = integrate_classical(state0, T, h, V=V)
        errors.append(abs(traj[-1][0] - math.exp(math.sqrt(12.0) / 12.0 * T)))
        drifts.append(max_classical_constraint_drift(traj, V=V))
    rate1 = rate(errors[0], errors[1])
    rate2 = rate(errors[1], errors[2])
    add_check(
        "Klassisch: de-Sitter-artiger Constraint mit RK4-Konvergenz",
        drifts[-1] < 1.0e-5 and rate1 > 3.5 and rate2 > 3.5,
        f"Constraint-Drift={drifts[-1]:.3e}; RK4-Raten={rate1:.2f}, {rate2:.2f}",
    )

    # 3. Ebenezelle auf dem Superraum: Q = 0, exakt gefuehrte Trajektorie.
    C = -1.0
    b_plane = -sqrt12
    plane_states = integrate_guidance((1.0, 0.0), T, 1.0e-3, b_plane, C)
    plane_residuals = [
        hamilton_constraint(a, phi, b_plane * a, C) for a, phi in plane_states
    ]
    plane_max_residual = max(abs(r) for r in plane_residuals)
    plane_a = plane_states[-1][0]
    plane_a_exact = 1.0 + sqrt12 / 12.0 * T
    add_check(
        "Bohm: Ebenezelle mit Q=0 und klassischer Fuhrung",
        plane_max_residual < 1.0e-12 and abs(plane_a - plane_a_exact) < 1.0e-12,
        f"max H_Q={plane_max_residual:.3e}; a-Endfehler={abs(plane_a - plane_a_exact):.3e}",
    )

    # 4. Glattes gaussisches Wellenpaket und klassischer Grenzfall.
    model = {
        "a0": 1.0,
        "phi0": 0.0,
        "a_center": 1.0,
        "phi_center": 0.0,
        "sigma_a": 0.5,
        "sigma_phi": 0.8,
        "C": C,
        "T": 0.2,
        "h": 2.0e-4,
        "hbar_values": [0.05, 0.025, 0.0125],
        "constraint": "H_C=-p_a^2/(24a)+a p_phi^2/2+a^3 V-k a+Lambda a^3/3",
        "quantum_constraint": "H_Q=H_C(a,phi,b*a,C)+Q",
        "guidance": "dot a=-b/12, dot phi=a*C",
        "quantum_potential": "Q=-hbar^2 (partial_a^2+partial_phi^2)R/(2R), R=exp(-F)",
    }
    data["model"] = model
    hbar_runs = []
    b_classical = -math.sqrt(12.0 * C * C)

    for hbar in model["hbar_values"]:
        Q0 = gaussian_quantum_potential(
            model["a0"], model["phi0"], hbar, model["a_center"], model["phi_center"], model["sigma_a"], model["sigma_phi"]
        )
        inside = 12.0 * C * C + 24.0 * Q0 / model["a0"]
        if inside <= 0.0:
            hbar_runs.append({"hbar": hbar, "failed": "no real initial momentum"})
            continue
        b = -math.sqrt(inside)
        initial_h, initial_q = quantum_constraint(
            model["a0"], model["phi0"], b, C, hbar, model
        )
        initial_residual = abs(initial_h + initial_q)

        quantum_states = integrate_guidance((model["a0"], model["phi0"]), model["T"], model["h"], b, C)
        classical_states = integrate_guidance(
            (model["a0"], model["phi0"]), model["T"], model["h"], b_classical, C
        )

        max_q = 0.0
        max_drift = 0.0
        traj_err = 0.0
        for q_state, c_state in zip(quantum_states, classical_states):
            a, phi = q_state
            h_val, q_val = quantum_constraint(a, phi, b, C, hbar, model)
            max_q = max(max_q, abs(q_val))
            max_drift = max(max_drift, abs(h_val + q_val))
            traj_err = max(traj_err, math.hypot(a - c_state[0], phi - c_state[1]))

        hbar_runs.append(
            {
                "hbar": hbar,
                "b": b,
                "b_classical": b_classical,
                "initial_residual": initial_residual,
                "max_Q": max_q,
                "max_Hq_drift": max_drift,
                "trajectory_error_vs_classical": traj_err,
            }
        )

    data["results"]["hbar_runs"] = hbar_runs
    initial_residual_max = max(run["initial_residual"] for run in hbar_runs)
    q0, q1, q2 = (run["max_Q"] for run in hbar_runs)
    err0, err1, err2 = (run["trajectory_error_vs_classical"] for run in hbar_runs)
    drift0, drift1, drift2 = (run["max_Hq_drift"] for run in hbar_runs)
    q_rate1 = rate(q0, q1)
    q_rate2 = rate(q1, q2)
    traj_rate1 = rate(err0, err1)
    traj_rate2 = rate(err1, err2)
    drift_rate1 = rate(drift0, drift1)
    drift_rate2 = rate(drift1, drift2)

    add_check(
        "Bohm: gaussischer Wellenpaket-Anfang erfuellt H_Q=0",
        initial_residual_max < 1.0e-12,
        f"max Anfangsresiduum={initial_residual_max:.3e}",
    )
    add_check(
        "Bohm: Drift der quantenkorrigierten Zwangsfunktion",
        drift2 < 1.0e-4 and drift_rate1 > 1.8 and drift_rate2 > 1.8,
        f"Drift={drift0:.3e}, {drift1:.3e}, {drift2:.3e}; Raten={drift_rate1:.2f}, {drift_rate2:.2f}",
    )
    add_check(
        "Bohm: klassischer Grenzfall des Quantenpotentials",
        q_rate1 > 1.8 and q_rate2 > 1.8,
        f"max |Q|={q0:.3e}, {q1:.3e}, {q2:.3e}; Raten={q_rate1:.2f}, {q_rate2:.2f}",
    )
    add_check(
        "Bohm: Trajektorie naehert klassischer Loesung mit hbar^2",
        traj_rate1 > 1.8 and traj_rate2 > 1.8,
        f"Trajektoriefehler={err0:.3e}, {err1:.3e}, {err2:.3e}; Raten={traj_rate1:.2f}, {traj_rate2:.2f}",
    )

    all_pass = all(check["passed"] for check in checks)
    data["summary"] = {
        "all_pass": all_pass,
        "checks": len(checks),
        "failed": sum(1 for check in checks if not check["passed"]),
        "classical_constraint_drift": {
            "massless": drift_massless,
            "desitter": drifts[-1],
        },
        "plane_wave_residual": plane_max_residual,
    }

    JSON_PATH.write_text(json.dumps(data, indent=2, ensure_ascii=True), encoding="utf-8")

    md = []
    md.append("# Bohmsche FLRW-Dynamik: reduziertes Modell\n")
    md.append(
        "Dieses Modell ist eine explizite Minisuperspace-Reduktion. Es quantisiert nicht die volle Geometrie, "
        "definiert keinen vollständigen Wheeler-DeWitt-Operator und beansprucht keine globale Existenzlösung.\n"
    )
    md.append("## Modell\n")
    md.append(
        "Der Hamilton-Zwang lautet in Einheiten $\\kappa=1$\n"
        "$$\nH_C(a,\\phi,p_a,p_\\phi)=-\\frac{p_a^2}{24a}+\\frac{a p_\\phi^2}{2}+a^3V(\\phi)-ka+\\frac{\\Lambda a^3}{3}.\n$$\n"
        "Für die Bohmsche Superraumführung wird ein Polaransatz $\\Psi=R e^{iS/\\hbar}$ gewählt. Mit "
        "$p_a=\\partial_aS$ und $p_\\phi=\\partial_\\phiS$ ist die Führung\n"
        "$$\\dot a=-p_a/(12a),\\qquad \\dot\\phi=a p_\\phi.$$\n"
        "Das Quantenpotential benutzt die ausdrücklich einfache Superraum-Laplace-Struktur "
        "$\\Delta_D=\\partial_a^2+\\partial_\\phi^2$:\n"
        "$$Q=-\\frac{\\hbar^2}{2}\\frac{\\Delta_D R}{R},\\qquad H_Q=H_C+Q.$$\n"
        "Diese Wahl ist keine Behauptung einer vollständigen Faktorordnung der Wheeler-DeWitt-Gleichung.\n"
    )
    md.append("## Tatsächliche Prüfergebnisse\n")
    md.append("| Prüfung | Ergebnis |\n|---|---:|\n")
    md.append(f"| Masseloser Skalar: Constraint-Drift | ${drift_massless:.3e}$ |\n")
    md.append(f"| Masseloser Skalar: Endpunktfehler | ${final_err_massless:.3e}$ |\n")
    md.append(f"| De-Sitter-Art: Constraint-Drift | ${drifts[-1]:.3e}$ |\n")
    md.append(f"| De-Sitter-Art: RK4-Konvergenzraten | ${rate1:.2f}$ und ${rate2:.2f}$ |\n")
    md.append(f"| Ebenezelle: maximales $H_Q$ | ${plane_max_residual:.3e}$ |\n")
    md.append(f"| Wellenpaket: maximales Anfangsresiduum | ${initial_residual_max:.3e}$ |\n\n")
    md.append("## Klassischer Grenzfall des Wellenpakets\n")
    md.append(
        "| $\\hbar$ | max Q | Trajektoriefehler gegen klassisch | Drift von $H_Q$ |\n"
        "|---:|---:|---:|---:|\n"
    )
    for run in hbar_runs:
        md.append(
            f"| {run['hbar']} | ${run['max_Q']:.3e}$ | "
            f"${run['trajectory_error_vs_classical']:.3e}$ | ${run['max_Hq_drift']:.3e}$ |\n"
        )
    md.append(
        f"\nDie Konvergenzraten sind ${q_rate1:.2f}$ und ${q_rate2:.2f}$ für $|Q|$ sowie "
        f"${traj_rate1:.2f}$ und ${traj_rate2:.2f}$ für den Trajektoriefehler. Dies stützt im definierten "
        "Modell den $\\hbar^2$-Grenzfall.\n"
    )
    md.append("## Grenzen\n")
    md.append(
        "Der Bericht prüft nur die festgelegten klassischen Hamilton-Gleichungen und die festgelegte "
        "Superraumführung. Kontinuität und Äquivarianz der Gauss-Amplitude sind hier nicht "
        "als vollständige Wheeler-DeWitt-Aussage geprüft. Die Regularisierung des Operators, das Maß auf "
        "dem Superraum und globale Trajektoriensaussagen bleiben offen.\n"
    )
    MARKDOWN_PATH.write_text("".join(md), encoding="utf-8")
    print("\n==== SUMMARY ====")
    print(f"all_pass = {all_pass}")
    print(f"checks   = {len(checks)}")
    print(f"json     = {JSON_PATH}")
    print(f"markdown = {MARKDOWN_PATH}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(run())
