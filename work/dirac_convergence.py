"""Dirac discretization audit: vary time step and spatial grid independently."""
import argparse
from pathlib import Path
from bohm_numerics import output_paths, sha256, source_hashes, validity, write_result
from dirac_numerics import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Ns", type=int, nargs="+", default=[256,512,1024])
    parser.add_argument("--dts", type=float, nargs="+", default=[0.01,0.005,0.0025])
    parser.add_argument("--T", type=float, default=5.0)
    parser.add_argument("--L", type=float, default=40.0)
    parser.add_argument("--out", default="dirac_convergence_corrected_metrics.json")
    args = parser.parse_args()
    output_paths(args.out)
    if len(set(args.Ns)) < 3 or len(set(args.dts)) < 3:
        raise ValueError("At least three distinct Ns and dts are required")
    finest_N, finest_dt = max(args.Ns), min(args.dts)
    configs = sorted({(N,finest_dt) for N in args.Ns} | {(finest_N,dt) for dt in args.dts})
    records, arrays = [], {}
    for index, (N,dt) in enumerate(configs):
        record, fields = run(N,dt,L=args.L,T=args.T)
        record["array_prefix"] = f"run_{index}"
        records.append(record)
        arrays.update({f"run_{index}_{key}":value for key,value in fields.items()})
    hashes = {**source_hashes(__file__), "dirac_numerics.py": sha256(Path(__file__).with_name("dirac_numerics.py"))}
    write_result(args.out, {"kind": "dirac_convergence",
        "model": "central_difference_Dirac_CN_positive_energy_lattice_projection_v2",
        "source_sha256": hashes, "config": {k:v for k,v in vars(args).items() if k != "out"},
        "reference": "Separate fixed-grid time refinement and fixed-time-step space refinement",
        "metrics": {"validated_convergence_order": None}, "runs": records,
        "diagnostics": {"time_series_fixed_N": finest_N, "space_series_fixed_dt": finest_dt,
                        "interpretation": "Exact discrete link continuity does not establish continuum-solution convergence."},
        "validity": validity([])}, arrays)


if __name__ == "__main__":
    main()
