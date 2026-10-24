"""Doppelspalt (2D, Bohm) - korrigierte, ausfuehrbare Version.

Nutzt das reparierte numerische Modell aus bohm_numerics.py:
  - Zellzentriertes periodisches Gitter (FFT, spektrale Ableitungen).
  - Strang-Split-Schroedinger-Propagation (zweite Ordnung).
  - RK4-Teilchenstufen an tatsaechlichen Halbzeit-Wellen (kein linearer Blend).
  - Bohm-Fuehrung v = Im(conj(psi)*grad(psi))/abs(psi)^2 ohne eps-Add auf psi.
  - Geschwindigkeitsbegrenzung nur optional (--vmax, gezahlt, invalidiert Abnahme).
  - Detektorstrom getrennt vorzeichenbehaftet/positiv/negativ integriert.
  - Randmasse und Normdrift gemessen (periodisches Gebiet).

Erzeugt:
  plots/double_slit_v3_density.png  : |psi|^2 am Ende + Bohmsche Trajektorien.
  plots/double_slit_v3_detector.png : Erstankunft vs. positiver Fluss (Bins).
  double_slit_v3_metrics.json       : Metriken + Arrays (NPZ).

Hinweis: Ein einzelner Lauf ist kein Beweis quantitativer Korrektheit. Die
Genauigkeits-, Konvergenz- und Seed-Pruefungen bleiben offen (Publikation).
"""
import argparse
import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bohm_numerics import (Grid, Evolution, fixed_wall, initial_wave, sample_gaussian,
    step_count, positive, overlap_weights, relative_l2, write_result, metadata, validity,
    output_paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Lx", type=float, default=60.0)
    parser.add_argument("--Ly", type=float, default=80.0)
    parser.add_argument("--Nx", type=int, default=1201)
    parser.add_argument("--Ny", type=int, default=801)
    parser.add_argument("--T", type=float, default=6.0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--Np", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--det_x", type=float, default=10.0)
    parser.add_argument("--n_bins", type=int, default=60)
    parser.add_argument("--vmax", type=float, default=500.0,
                        help="0 deaktiviert Begrenzung; >0 wird gezahlt und invalidiert die Abnahme.")
    parser.add_argument("--sigma_x", type=float, default=0.15)
    parser.add_argument("--sigma_y", type=float, default=0.30)
    parser.add_argument("--wall_v0", type=float, default=120.0)
    parser.add_argument("--edge_width", type=float, default=1.0)
    parser.add_argument("--boundary_tolerance", type=float, default=1e-8)
    parser.add_argument("--out", default="double_slit_v3_metrics.json")
    parser.add_argument("--traj_every", type=int, default=5)
    args = parser.parse_args()

    output_paths(args.out)
    positive(args.det_x, "det_x")
    if not -args.Lx / 2 < args.det_x < args.Lx / 2 or args.n_bins < 2:
        raise ValueError("det_x muss strikt im Gebiet liegen und n_bins >= 2")
    total = step_count(args.T, args.dt)
    grid = Grid(args.Lx, args.Ly, args.Nx, args.Ny)
    potential = fixed_wall(grid, args.sigma_x, args.sigma_y, args.wall_v0)
    evolution = Evolution(grid, potential, args.dt, args.vmax,
                          edge_width=args.edge_width, boundary_tolerance=args.boundary_tolerance)
    psi = initial_wave(grid)
    positions, proposals = sample_gaussian(grid, args.Np, args.seed)
    derivatives = grid.derivatives(psi)
    evolution.inspect_wave(psi)

    # Detektor und Bins
    bin_edges = np.linspace(-args.Ly / 2, args.Ly / 2, args.n_bins + 1)
    weights = overlap_weights(grid.y_edges, bin_edges)
    right_widths = overlap_weights(grid.x_edges, np.array([args.det_x, args.Lx / 2]))[0]
    probe = np.column_stack((np.full(grid.Ny, args.det_x), grid.y))

    def current(wave, derivative):
        wp = grid.interpolate(wave, probe)
        dp = grid.interpolate(derivative, probe)
        return np.imag(np.conj(wp) * dp)  # hbar=m=1

    initial_right_mass = float(np.sum(abs(psi) ** 2 * right_widths[None, :]) * grid.dy)
    initially_right = positions[:, 0] >= args.det_x
    first_hit = np.zeros(args.Np, dtype=bool)
    first_y = np.full(args.Np, np.nan)
    forward_counts = np.zeros(args.Np, dtype=np.int64)
    backward_counts = np.zeros(args.Np, dtype=np.int64)
    signed_profile = np.zeros(grid.Ny)
    positive_profile = np.zeros(grid.Ny)
    negative_profile = np.zeros(grid.Ny)

    # Trajektorien-Speicher (alle traj_every Schritte)
    traj_store = []

    previous_current = current(psi, derivatives[0])
    for index in range(total):
        new, new_derivatives, updated = evolution.step(psi, derivatives, positions)

        # Kreuzungszaehlung an der Detektorebene
        right = (positions[:, 0] <= args.det_x) & (updated[:, 0] > args.det_x)
        left = (positions[:, 0] > args.det_x) & (updated[:, 0] <= args.det_x)
        forward_counts += right
        backward_counts += left
        new_first = right & ~first_hit & ~initially_right
        ids = np.flatnonzero(new_first)
        if len(ids):
            fraction = (args.det_x - positions[ids, 0]) / (updated[ids, 0] - positions[ids, 0])
            first_y[ids] = positions[ids, 1] + fraction * (updated[ids, 1] - positions[ids, 1])
            first_hit[ids] = True

        next_current = current(new, new_derivatives[0])
        signed_profile += args.dt * (previous_current + next_current) / 2
        positive_profile += args.dt * (np.maximum(previous_current, 0) + np.maximum(next_current, 0)) / 2
        negative_profile += args.dt * (np.maximum(-previous_current, 0) + np.maximum(-next_current, 0)) / 2

        psi, derivatives, positions = new, new_derivatives, updated
        previous_current = next_current
        if index % args.traj_every == 0:
            traj_store.append(positions.copy())

    # --- Auswertung ---
    signed_mass = weights @ signed_profile
    positive_mass = weights @ positive_profile
    negative_mass = weights @ negative_profile
    hist, _ = np.histogram(first_y[first_hit], bins=bin_edges)
    observed_mass = hist / args.Np
    n_hits = int(np.count_nonzero(first_hit))
    n_in_aperture = int(np.sum(hist))
    mass_positive = float(np.sum(positive_mass))
    obs_norm = observed_mass / np.sum(observed_mass) if n_in_aperture else None
    exp_norm = positive_mass / mass_positive if mass_positive > 0 else None
    shape_error = relative_l2(obs_norm, exp_norm) if obs_norm is not None and exp_norm is not None else None
    abs_mass_error = relative_l2(observed_mass, positive_mass) if n_in_aperture and mass_positive > 0 else None

    issues = evolution.issues()
    if initial_right_mass > args.boundary_tolerance or np.any(initially_right):
        issues.append("initial_probability_to_right_of_detector")
    if int(np.sum(backward_counts)) or int(np.sum(np.maximum(forward_counts - 1, 0))) or \
       float(np.sum(negative_mass)) > args.boundary_tolerance:
        issues.append("positive_current_not_established_as_first_arrival_reference")
    if n_in_aperture != n_hits:
        issues.append("first_arrivals_outside_detector_aperture")
    if shape_error is None:
        issues.append("conditional_shape_undefined_without_hits_and_reference_mass")

    diagnostics = {**evolution.diag, "sampling_proposals": proposals,
        "n_particles": args.Np, "n_hits": n_hits, "n_hits_in_aperture": n_in_aperture,
        "hit_prob_empirical": n_in_aperture / args.Np,
        "all_crossing_hit_prob_empirical": n_hits / args.Np,
        "flux_total_mass": mass_positive,
        "hit_prob_diff": n_in_aperture / args.Np - mass_positive,
        "net_flux_mass": float(np.sum(signed_mass)), "pos_flux_mass": mass_positive,
        "negative_flux_mass": float(np.sum(negative_mass)), "initial_right_mass": initial_right_mass,
        "mass_right_end": float(np.sum(abs(psi) ** 2 * right_widths[None, :]) * grid.dy),
        "n_cross_right_total": int(np.sum(forward_counts)),
        "n_cross_left_total": int(np.sum(backward_counts)),
        "n_repeated_right_crossings": int(np.sum(np.maximum(forward_counts - 1, 0))),
        "n_multi_cross": int(np.count_nonzero(forward_counts + backward_counts > 1)),
        "det_x_actual": args.det_x,
        "norm_end": float(np.sum(abs(psi) ** 2) * grid.dx * grid.dy),
        "edge_width": args.edge_width}

    # --- Plots ---
    os.makedirs("plots", exist_ok=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_w = bin_edges[1] - bin_edges[0]

    # Plot 1: Dichte + Trajektorien (Uebersicht und vergroesserter Durchlassbereich)
    density = abs(psi.T) ** 2
    vmax_d = float(np.percentile(density, 98))
    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    ax = axes[0]
    ax.imshow(density, extent=[-args.Lx / 2, args.Lx / 2, -args.Ly / 2, args.Ly / 2],
              origin='lower', cmap='Blues', alpha=0.85, aspect='auto', vmax=vmax_d)
    for k in range(min(args.Np, 800)):
        xs = np.array([p[k, 0] for p in traj_store])
        ys = np.array([p[k, 1] for p in traj_store])
        ax.plot(xs, ys, 'r-', lw=0.5, alpha=0.45)
    ax.axvline(0, color='k', lw=3)
    ax.axvline(args.det_x, color='g', lw=2, ls='--')
    ax.set_xlim(-args.Lx / 2, args.Lx / 2); ax.set_ylim(-args.Ly / 2, args.Ly / 2)
    ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
    ax.set_title('Doppelspalt (Bohm): $|\\psi|^2$ und Trajektorien (Uebersicht)')
    # Vergroesserung: Durchlassbereich x>0 und Interferenz
    ax = axes[1]
    xzoom = slice(int(np.searchsorted(grid.x, 0.0)), int(np.searchsorted(grid.x, args.Lx / 2)))
    dzoom = abs(psi[:, xzoom].T) ** 2
    ax.imshow(dzoom, extent=[grid.x[xzoom][0], grid.x[xzoom][-1], -args.Ly / 2, args.Ly / 2],
              origin='lower', cmap='Blues', alpha=0.85, aspect='auto',
              vmax=float(np.percentile(dzoom, 99)))
    for k in range(min(args.Np, 800)):
        xs = np.array([p[k, 0] for p in traj_store])
        ys = np.array([p[k, 1] for p in traj_store])
        ax.plot(xs, ys, 'r-', lw=0.5, alpha=0.45)
    ax.axvline(args.det_x, color='g', lw=2, ls='--')
    ax.set_xlim(grid.x[xzoom][0], grid.x[xzoom][-1]); ax.set_ylim(-args.Ly / 2, args.Ly / 2)
    ax.set_xlabel('$x$ (Durchlassbereich)'); ax.set_ylabel('$y$')
    ax.set_title('Vergroesserung: Durchlassbereich und Interferenz')
    fig.tight_layout(); fig.savefig('plots/double_slit_v3_density.png', dpi=150); plt.close(fig)

    # Plot 2: Erstankunft vs. positiver Fluss
    fig, ax = plt.subplots(figsize=(8, 5))
    if obs_norm is not None:
        ax.bar(bin_centers, obs_norm, width=bin_w * 0.8, alpha=0.6,
               label='Bohm-Erstankunft (Histogramm)')
    if exp_norm is not None:
        ax.plot(bin_centers, exp_norm, 'r-', lw=2, label='Positiver Fluss (Referenz)')
    ax.set_xlabel('$y$ (Erstankunft)'); ax.set_ylabel('normierte Wahrscheinlichkeit')
    title = f'Doppelspalt: Erstankunft vs. Fluss'
    if shape_error is not None:
        title += f'  (rel. L2={shape_error:.4f})'
    ax.set_title(title)
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig('plots/double_slit_v3_detector.png', dpi=150); plt.close(fig)

    # --- Ergebnis ---
    arrays = {"bin_edges": bin_edges, "observed_absolute_mass": observed_mass,
              "positive_flux_mass": positive_mass, "negative_flux_mass": negative_mass,
              "signed_flux_mass": signed_mass, "first_hit_y": first_y[first_hit]}
    if obs_norm is not None:
        arrays["obs_norm"] = obs_norm
    if exp_norm is not None:
        arrays["exp_norm"] = exp_norm
    result = {"kind": "detector", **metadata(args, __file__),
        "observable": "detector_first_arrival_vs_forward_crossing_flux",
        "reference": "Zeitintegrierter positiver Wahrscheinlichkeitsfluss an x=det_x; als Erstankunftsreferenz nur bei nachgewiesen vernachlaessigbarem Rueckfluss, Wiederkreuzungen und Randeinfluss.",
        "metrics": {"arrival_rel_l2_error": shape_error,
                    "arrival_absolute_mass_rel_l2_error": abs_mass_error,
                    "arrival_grid_change": None, "arrival_monte_carlo_std": None},
        "validity": validity(issues), "diagnostics": diagnostics}
    write_result(args.out, result, arrays)
    print(json.dumps({"metrics": result["metrics"], "validity": result["validity"],
                      "diagnostics": diagnostics}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
