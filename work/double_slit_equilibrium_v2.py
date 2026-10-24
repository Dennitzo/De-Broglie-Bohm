"""Fixed-time Born diagnostic. Formal source repair; no new result claimed.

Raw bin mass and equally smoothed mass are distinct reported observables.
Dense fields go to immutable compressed NPZ, never into tool-result JSON.
"""
import argparse
import numpy as np
from bohm_numerics import (Grid, Evolution, common_parser, fixed_wall, initial_wave,
    sample_gaussian, metadata, output_paths, positive, relative_l2, smooth_mass,
    step_count, validity, write_result)

REFERENCE = "Numerische Schrodinger-Dichte auf identischen Zellen und zu exakt bezeichneten Zeiten; raeumliche, zeitliche und statistische Konvergenz separat erforderlich."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    common_parser(parser, 50000)
    parser.add_argument("--sigma_kde", type=float, default=0.2)
    parser.add_argument("--wall_v0", type=float, default=120.0)
    # Explicit physical width, matching old fine-grid effective width. Never
    # widened automatically. This is not the old nominal width of 0.05.
    parser.add_argument("--wall_sigx", type=float, default=0.15)
    parser.add_argument("--wall_sigy", type=float, default=0.30)
    parser.add_argument("--check_times", type=float, nargs="+", default=[3.0, 4.0, 5.0])
    parser.add_argument("--out", default="equilibrium_v2_corrected_metrics.json")
    args = parser.parse_args()
    output_paths(args.out)
    positive(args.sigma_kde, "sigma_kde", allow_zero=True)
    total = step_count(args.T, args.dt)
    if len(set(args.check_times)) != len(args.check_times):
        raise ValueError("check_times must be unique")
    requested = {step_count(t, args.dt, "measurement time"): t for t in args.check_times}
    if any(index > total for index in requested):
        raise ValueError("A required measurement time lies beyond T")
    grid = Grid(args.Lx, args.Ly, args.Nx, args.Ny)
    potential = fixed_wall(grid, args.wall_sigx, args.wall_sigy, args.wall_v0)
    evolution = Evolution(grid, potential, args.dt, args.vmax,
                          edge_width=args.edge_width, boundary_tolerance=args.boundary_tolerance)
    psi = initial_wave(grid)
    positions, proposals = sample_gaussian(grid, args.Np, args.seed)
    derivatives = grid.derivatives(psi)
    evolution.inspect_wave(psi)
    arrays, measurements = {}, {}

    def record(step, wave, points):
        if step not in requested:
            return
        time = requested[step]
        wrapped = grid.wrap(points)
        counts = np.histogram2d(wrapped[:, 0], wrapped[:, 1],
                               bins=(grid.x_edges, grid.y_edges))[0].T
        empirical = counts / args.Np
        reference = abs(wave)**2 * grid.dx * grid.dy
        empirical_smoothed = smooth_mass(empirical, grid, args.sigma_kde)
        reference_smoothed = smooth_mass(reference, grid, args.sigma_kde)
        key = f"step_{step}"
        arrays[key+"_empirical_mass"] = empirical
        arrays[key+"_reference_mass"] = reference
        arrays[key+"_empirical_smoothed_mass"] = empirical_smoothed
        arrays[key+"_reference_smoothed_mass"] = reference_smoothed
        density = abs(wave)**2
        edge_density = density[evolution.edge_mask]
        measurements[str(time)] = {
            "requested_time": time, "actual_time": step*args.dt,
            "array_prefix": key, "raw_rel_l2": relative_l2(empirical, reference),
            "equally_smoothed_rel_l2": relative_l2(empirical_smoothed, reference_smoothed),
            "empirical_total_mass": float(np.sum(empirical)),
            "reference_total_mass": float(np.sum(reference)),
            "outside_unwrapped_particles": int(np.count_nonzero(grid.outside(points))),
            "raw_count": int(np.sum(counts)),
            "edge_amplitude": float(np.sqrt(np.max(edge_density))),
            "edge_mass": float(np.sum(edge_density) * grid.dx * grid.dy)}

    record(0, psi, positions)
    for completed in range(1, total+1):
        psi, derivatives, positions = evolution.step(psi, derivatives, positions)
        record(completed, psi, positions)
    issues = evolution.issues()
    if len(measurements) != len(requested):
        issues.append("required_measurement_missing")
    arrays.update(x_edges=grid.x_edges, y_edges=grid.y_edges)
    raw = [item["raw_rel_l2"] for item in measurements.values()]
    smoothed = [item["equally_smoothed_rel_l2"] for item in measurements.values()]
    result = {"kind": "equilibrium", **metadata(args, __file__),
        "observable": "fixed_time_probability_mass_rel_l2_error", "reference": REFERENCE,
        "metrics": {"rel_l2_error": max(raw) if raw and None not in raw else None,
                    "equally_smoothed_rel_l2_error": max(smoothed) if smoothed and None not in smoothed else None,
                    "grid_change": None, "monte_carlo_std": None},
        "measurements": measurements, "validity": validity(issues),
        "diagnostics": {**evolution.diag, "sampling_proposals": proposals,
                        "norm_end": float(np.sum(abs(psi)**2)*grid.dx*grid.dy),
                        "clipping_counter_unit": "particle-stage evaluations before modification",
                        "raw_metric_note": "Unsmoothed histogram sampling noise is included; fixed-bandwidth error is reported separately."}}
    write_result(args.out, result, arrays)


if __name__ == "__main__":
    main()
