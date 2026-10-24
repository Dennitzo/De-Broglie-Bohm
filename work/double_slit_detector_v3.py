"""First-arrival and signed/positive/negative current diagnostics, kept distinct.

Positive current counts forward crossings. Interpreting it as first-hit density
requires backflow, recrossing and boundary checks; no exact detector law is assumed.
"""
import argparse
import numpy as np
from bohm_numerics import (Grid, Evolution, common_parser, fixed_wall, initial_wave,
    sample_gaussian, metadata, output_paths, overlap_weights, positive, relative_l2,
    step_count, validity, write_result)

REFERENCE = "Zeitintegrierter positiver Wahrscheinlichkeitsfluss an x=10; als Erstankunftsreferenz nur bei nachgewiesen vernachlaessigbarem Rueckfluss, Wiederkreuzungen und Randeinfluss."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    common_parser(parser, 20000)
    parser.add_argument("--det_x", type=float, default=10.0)
    parser.add_argument("--n_bins", type=int, default=40)
    parser.add_argument("--sigma_x", type=float, default=0.15)
    parser.add_argument("--sigma_y", type=float, default=0.30)
    parser.add_argument("--wall_v0", type=float, default=120.0)
    parser.add_argument("--backflow_tolerance", type=float, default=1e-6)
    parser.add_argument("--out", default="detector_v3_corrected_metrics.json")
    args = parser.parse_args()
    output_paths(args.out)
    total = step_count(args.T, args.dt)
    positive(args.backflow_tolerance, "backflow_tolerance", allow_zero=True)
    if not -args.Lx/2 < args.det_x < args.Lx/2 or args.n_bins < 2:
        raise ValueError("Detector must lie strictly inside the domain and n_bins >= 2")
    grid = Grid(args.Lx, args.Ly, args.Nx, args.Ny)
    potential = fixed_wall(grid, args.sigma_x, args.sigma_y, args.wall_v0)
    evolution = Evolution(grid, potential, args.dt, args.vmax,
                          edge_width=args.edge_width, boundary_tolerance=args.boundary_tolerance)
    psi = initial_wave(grid)
    positions, proposals = sample_gaussian(grid, args.Np, args.seed)
    derivatives = grid.derivatives(psi)
    evolution.inspect_wave(psi)
    probe = np.column_stack((np.full(grid.Ny, args.det_x), grid.y))
    bin_edges = np.linspace(-args.Ly/2, args.Ly/2, args.n_bins+1)
    weights = overlap_weights(grid.y_edges, bin_edges)
    right_widths = overlap_weights(grid.x_edges, np.array([args.det_x, args.Lx/2]))[0]

    def current(wave, derivative):
        wp = grid.interpolate(wave, probe)
        dp = grid.interpolate(derivative, probe)
        return np.imag(np.conj(wp)*dp)  # hbar=m=1 in this experiment

    initial_right_mass = float(np.sum(abs(psi)**2 * right_widths[None, :]) * grid.dy)
    initially_right = positions[:, 0] >= args.det_x
    first_hit = np.zeros(args.Np, dtype=bool)
    first_y = np.full(args.Np, np.nan)
    first_time = np.full(args.Np, np.nan)
    forward_counts = np.zeros(args.Np, dtype=np.int64)
    backward_counts = np.zeros(args.Np, dtype=np.int64)
    signed_profile = np.zeros(grid.Ny)
    positive_profile = np.zeros(grid.Ny)
    negative_profile = np.zeros(grid.Ny)
    previous_current = current(psi, derivatives[0])
    for index in range(total):
        new, new_derivatives, updated = evolution.step(psi, derivatives, positions)
        right = (positions[:, 0] <= args.det_x) & (updated[:, 0] > args.det_x)
        left = (positions[:, 0] > args.det_x) & (updated[:, 0] <= args.det_x)
        forward_counts += right
        backward_counts += left
        new_first = right & ~first_hit & ~initially_right
        ids = np.flatnonzero(new_first)
        if len(ids):
            fraction = (args.det_x-positions[ids, 0]) / (updated[ids, 0]-positions[ids, 0])
            first_y[ids] = positions[ids, 1]+fraction*(updated[ids, 1]-positions[ids, 1])
            first_time[ids] = (index+fraction)*args.dt
            first_hit[ids] = True
        next_current = current(new, new_derivatives[0])
        signed_profile += args.dt * (previous_current+next_current) / 2
        positive_profile += args.dt * (np.maximum(previous_current, 0)+np.maximum(next_current, 0)) / 2
        negative_profile += args.dt * (np.maximum(-previous_current, 0)+np.maximum(-next_current, 0)) / 2
        psi, derivatives, positions = new, new_derivatives, updated
        previous_current = next_current
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
    issues = evolution.issues()
    reverse_total = int(np.sum(backward_counts))
    repeated_forward = int(np.sum(np.maximum(forward_counts-1, 0)))
    negative_total = float(np.sum(negative_mass))
    if initial_right_mass > args.backflow_tolerance or np.any(initially_right):
        issues.append("initial_probability_to_right_of_detector")
    if reverse_total or repeated_forward or negative_total > args.backflow_tolerance:
        issues.append("positive_current_not_established_as_first_arrival_reference")
    if n_in_aperture != n_hits:
        issues.append("first_arrivals_outside_detector_aperture")
    if shape_error is None:
        issues.append("conditional_shape_undefined_without_hits_and_reference_mass")
    diagnostics = {**evolution.diag, "sampling_proposals": proposals,
        "n_particles": args.Np, "n_hits": n_hits, "n_hits_in_aperture": n_in_aperture,
        "hit_prob_empirical": n_in_aperture/args.Np, "all_crossing_hit_prob_empirical": n_hits/args.Np,
        "flux_total_mass": mass_positive, "hit_prob_diff": n_in_aperture/args.Np-mass_positive,
        "net_flux_mass": float(np.sum(signed_mass)), "pos_flux_mass": mass_positive,
        "negative_flux_mass": negative_total, "initial_right_mass": initial_right_mass,
        "mass_right_end": float(np.sum(abs(psi)**2 * right_widths[None, :])*grid.dy),
        "n_cross_right_total": int(np.sum(forward_counts)), "n_cross_left_total": reverse_total,
        "n_repeated_right_crossings": repeated_forward,
        "n_multi_cross": int(np.count_nonzero(forward_counts+backward_counts > 1)),
        "det_x_actual": args.det_x, "norm_end": float(np.sum(abs(psi)**2)*grid.dx*grid.dy),
        "crossing_localization": "linear within each step; unresolved within-step recrossings require time refinement"}
    arrays = {"bin_edges": bin_edges, "observed_absolute_mass": observed_mass,
              "positive_flux_mass": positive_mass, "negative_flux_mass": negative_mass,
              "signed_flux_mass": signed_mass, "first_hit_time": first_time[first_hit],
              "first_hit_y": first_y[first_hit]}
    if obs_norm is not None:
        arrays["obs_norm"] = obs_norm
    if exp_norm is not None:
        arrays["exp_norm"] = exp_norm
    reference = REFERENCE if args.det_x == 10 else REFERENCE.replace("x=10;", f"x={args.det_x:g};")
    result = {"kind": "detector", **metadata(args, __file__),
        "observable": "detector_first_arrival_vs_forward_crossing_flux", "reference": reference,
        "metrics": {"arrival_rel_l2_error": shape_error,
                    "arrival_absolute_mass_rel_l2_error": relative_l2(observed_mass, positive_mass),
                    "arrival_probability_abs_error": abs(diagnostics["hit_prob_diff"]),
                    "arrival_grid_change": None, "arrival_monte_carlo_std": None},
        "validity": validity(issues), "diagnostics": diagnostics}
    write_result(args.out, result, arrays)


if __name__ == "__main__":
    main()
