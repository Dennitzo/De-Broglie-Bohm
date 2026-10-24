"""Free spreading Gaussian: center AND off-center Bohm trajectories.

This is a future control calculation, not supplied numerical evidence. Common
time trajectory arrays, rather than differences of scalar error norms, measure
time refinement. Spatial and domain-size convergence remain separate checks.
"""
import argparse
from pathlib import Path
import numpy as np
from bohm_numerics import (Grid, Evolution, initial_wave, output_paths, source_hashes,
                          step_count, validity, write_result)

REFERENCE = "Analytische Bohmsche Trajektorien eines freien, sich ausbreitenden Gauss-Wellenpakets: X(t)=q0+(p0/m)t+(X0-q0)*sqrt(1+(hbar*t/(m*sigma0^2))^2), fuer psi0 proportional exp(-(x-q0)^2/(2*sigma0^2)); Zentrum und nichtzentrale Anfangspositionen."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Lx", type=float, default=60.0)
    parser.add_argument("--Ly", type=float, default=30.0)
    parser.add_argument("--Nx", type=int, default=768)
    parser.add_argument("--Ny", type=int, default=384)
    parser.add_argument("--T", type=float, default=2.0)
    parser.add_argument("--dts", type=float, nargs="+", default=[0.02, 0.01, 0.005])
    parser.add_argument("--out", default="control_free_packet_corrected_metrics.json")
    parser.add_argument("--center_out", default="control_free_packet_center_corrected_metrics.json")
    args = parser.parse_args()
    output_paths(args.out)
    output_paths(args.center_out)
    if Path(args.out).resolve() == Path(args.center_out).resolve():
        raise ValueError("Center and general control need distinct output paths")
    dts = sorted(set(args.dts), reverse=True)
    if len(dts) < 3:
        raise ValueError("Supply at least three distinct time steps")
    for dt in dts:
        step_count(args.T, dt)
    for coarse, fine in zip(dts, dts[1:]):
        step_count(coarse, fine, "refinement ratio")
    grid = Grid(args.Lx, args.Ly, args.Nx, args.Ny)
    center, sigma, momentum = np.array([-12.0, 0.0]), 2.1, 5.0
    offsets = sigma * np.array([[0,0], [-2,0], [-1,0], [1,0], [2,0],
                               [0,-2], [0,-1], [0,1], [0,2], [1,1]], dtype=float)
    initial = center + offsets
    if np.any(grid.outside(initial)):
        raise ValueError("Domain does not contain all control trajectories at t=0")
    arrays = {"initial_positions": initial}
    errors, center_errors, differences, center_differences = [], [], [], []
    diagnostics, histories, issues = [], [], []
    for run_index, dt in enumerate(dts):
        count = step_count(args.T, dt)
        psi = initial_wave(grid)
        derivatives = grid.derivatives(psi)
        evolution = Evolution(grid, np.zeros((grid.Ny, grid.Nx)), dt, vmax=0.0)
        evolution.inspect_wave(psi)
        positions = initial.copy()
        history = np.empty((count+1, len(initial), 2))
        history[0] = positions
        for step in range(1, count+1):
            psi, derivatives, positions = evolution.step(psi, derivatives, positions)
            history[step] = positions
        times = np.arange(count+1) * dt
        dilation = np.sqrt(1 + (times/sigma**2)**2)  # hbar=m=1; amplitude width sigma
        analytic = center[None, None, :] + offsets[None, :, :] * dilation[:, None, None]
        analytic[:, :, 0] += momentum * times[:, None]
        errors.append(float(np.max(np.linalg.norm(history-analytic, axis=2))))
        center_errors.append(float(np.max(np.linalg.norm(history[:,0,:]-analytic[:,0,:],axis=1))))
        histories.append(history)
        arrays[f"run_{run_index}_times"] = times
        arrays[f"run_{run_index}_trajectories"] = history
        arrays[f"run_{run_index}_analytic"] = analytic
        diagnostics.append({"dt": dt, **evolution.diag})
        issues.extend(f"dt={dt}:{issue}" for issue in evolution.issues())
    for index in range(len(dts)-1):
        stride = step_count(dts[index], dts[index+1], "refinement ratio")
        aligned_fine = histories[index+1][::stride]
        if aligned_fine.shape != histories[index].shape:
            raise ValueError("Refinements do not share identical physical measurement times")
        differences.append(float(np.max(np.linalg.norm(histories[index]-aligned_fine, axis=2))))
        center_differences.append(float(np.max(np.linalg.norm(histories[index][:,0,:]-aligned_fine[:,0,:],axis=1))))
    result = {"kind": "control", "model": "free_gaussian_infinite_domain_analytic_reference",
        "source_sha256": source_hashes(__file__), "config": {k:v for k,v in vars(args).items() if k != "out"},
        "reference": REFERENCE, "observable": "trajectory_max_abs_error",
        "metrics": {"max_abs_error": errors[-1], "traj_grid_change": differences[-1]},
        "diagnostics": {"dts": dts, "max_error_by_dt": errors,
                        "successive_trajectory_distances": differences, "runs": diagnostics,
                        "coupled_time_order": 2, "units": "hbar=m=1; lengths and times nondimensional",
                        "pending": "Spatial and domain-size checks are not replaced by this time-step series."},
        "validity": validity(issues)}
    write_result(args.out, result, arrays)
    center_result = {**result,
        "reference": "analytische Bohmsche Trajektorie X(t)=x0+(p0/m)*t fuer freies Gauß-Wellenpaket (V=0)",
        "metrics": {"max_abs_error":center_errors[-1], "traj_grid_change":center_differences[-1]},
        "diagnostics": {**result["diagnostics"], "scope":"Packet center only; general spreading is a separate required report.",
                        "max_error_by_dt":center_errors, "successive_trajectory_distances":center_differences}}
    center_arrays = {key:(value[:,0,:] if value.ndim == 3 else value)
                     for key,value in arrays.items() if key != "initial_positions"}
    write_result(args.center_out, center_result, center_arrays)


if __name__ == "__main__":
    main()
