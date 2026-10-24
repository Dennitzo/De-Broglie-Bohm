"""Shared, explicitly versioned numerical contract for the repaired project.

Formal repair only: no new numerical result is supplied with this source.
The coupled method is second order (Strang propagation), even though particle
quadrature uses RK4 stages. Periodic boundaries approximate an open experiment
only after boundary and domain-size checks. No function runs on import.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import numpy as np

SCHEMA = "missum.bohmian.result.v2"
VERSION = "2026-10-13-formal-1"
MODEL = "gaussian-soft-double-slit-fixed-potential-v2"


def positive(value, name, allow_zero=False):
    if not math.isfinite(float(value)) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")


def step_count(time, dt, name="T"):
    positive(dt, "dt")
    positive(time, name, allow_zero=True)
    count = round(time / dt)
    if not math.isclose(count * dt, time, rel_tol=1e-11, abs_tol=1e-12):
        raise ValueError(f"{name}={time} is not exactly represented by dt={dt}; do not relabel a nearby time")
    return count


@dataclass
class Grid:
    Lx: float
    Ly: float
    Nx: int
    Ny: int

    def __post_init__(self):
        for key in ("Lx", "Ly"):
            positive(getattr(self, key), key)
        if self.Nx < 8 or self.Ny < 8:
            raise ValueError("Nx and Ny must be at least 8")
        self.dx, self.dy = self.Lx / self.Nx, self.Ly / self.Ny
        # FFT samples are cell centres. Histogram cells cover precisely Lx,Ly.
        self.x_edges = np.linspace(-self.Lx / 2, self.Lx / 2, self.Nx + 1)
        self.y_edges = np.linspace(-self.Ly / 2, self.Ly / 2, self.Ny + 1)
        self.x = (self.x_edges[:-1] + self.x_edges[1:]) / 2
        self.y = (self.y_edges[:-1] + self.y_edges[1:]) / 2
        self.kx = 2 * np.pi * np.fft.fftfreq(self.Nx, d=self.dx)
        self.ky = 2 * np.pi * np.fft.fftfreq(self.Ny, d=self.dy)
        self.k2 = self.kx[None, :] ** 2 + self.ky[:, None] ** 2

    def interpolate(self, field, positions):
        # Both wave and particles use the same periodic cell-centred lattice.
        fi = (positions[:, 0] - self.x[0]) / self.dx
        fj = (positions[:, 1] - self.y[0]) / self.dy
        bi, bj = np.floor(fi), np.floor(fj)
        ai, aj = fi - bi, fj - bj
        i, j = bi.astype(np.int64) % self.Nx, bj.astype(np.int64) % self.Ny
        ip, jp = (i + 1) % self.Nx, (j + 1) % self.Ny
        return (field[j, i] * (1-ai) * (1-aj) + field[j, ip] * ai * (1-aj)
                + field[jp, i] * (1-ai) * aj + field[jp, ip] * ai * aj)

    def derivatives(self, psi):
        transformed = np.fft.fft2(psi)
        return (np.fft.ifft2(1j * self.kx[None, :] * transformed),
                np.fft.ifft2(1j * self.ky[:, None] * transformed))

    def wrap(self, positions):
        length = np.array([self.Lx, self.Ly])
        return (positions + length / 2) % length - length / 2

    def outside(self, positions):
        return np.any((positions < -np.array([self.Lx, self.Ly]) / 2)
                      | (positions >= np.array([self.Lx, self.Ly]) / 2), axis=1)


def fixed_wall(grid, sigma_x, sigma_y, height=120.0, separation=3.0):
    for value, name in ((sigma_x, "wall_sigma_x"), (sigma_y, "wall_sigma_y"),
                        (height, "wall_height"), (separation, "slit_separation")):
        positive(value, name)
    if sigma_x < 3 * grid.dx * (1 - 1e-12) or sigma_y < 3 * grid.dy * (1 - 1e-12):
        raise ValueError("Underresolved fixed wall: require sigma_x >= 3*dx and sigma_y >= 3*dy. "
                         "Refine the grid; do not widen the physical wall to pass this check.")
    gx = np.exp(-grid.x[None, :] ** 2 / (2 * sigma_x ** 2))
    g1 = np.exp(-(grid.y[:, None] - separation / 2) ** 2 / (2 * sigma_y ** 2))
    g2 = np.exp(-(grid.y[:, None] + separation / 2) ** 2 / (2 * sigma_y ** 2))
    return height * gx * (1 - g1) * (1 - g2)


def initial_wave(grid, center=(-12.0, 0.0), sigma=2.1, momentum=5.0, hbar=1.0):
    psi = np.exp(-((grid.x[None, :] - center[0]) ** 2
                   + (grid.y[:, None] - center[1]) ** 2) / (2 * sigma ** 2))
    psi = psi * np.exp(1j * momentum * grid.x[None, :] / hbar)
    return psi / np.sqrt(np.sum(abs(psi) ** 2) * grid.dx * grid.dy)


def sample_gaussian(grid, count, seed, center=(-12.0, 0.0), sigma=2.1):
    """Exact continuous initial Gaussian, conditioned on the stated box.

    This avoids uniform rejection over a mostly empty large rectangle. The
    finite-grid quadrature approximation is recorded and still needs refinement.
    """
    if count < 1:
        raise ValueError("Np must be positive")
    rng = np.random.default_rng(seed)
    positions = np.empty((count, 2))
    filled = proposals = 0
    while filled < count:
        batch = rng.normal(center, sigma / np.sqrt(2), size=(count-filled, 2))
        proposals += len(batch)
        batch = batch[~grid.outside(batch)]
        positions[filled:filled+len(batch)] = batch
        filled += len(batch)
    return positions, proposals


class Evolution:
    """Physical half-time wave stages and auditable particle interventions."""
    def __init__(self, grid, potential, dt, vmax=0.0, hbar=1.0, mass=1.0,
                 edge_width=1.0, boundary_tolerance=1e-8):
        for value, name in ((dt, "dt"), (hbar, "hbar"), (mass, "mass"),
                            (edge_width, "edge_width")):
            positive(value, name)
        positive(vmax, "vmax", allow_zero=True)
        positive(boundary_tolerance, "boundary_tolerance", allow_zero=True)
        if edge_width >= min(grid.Lx, grid.Ly) / 2:
            raise ValueError("edge_width must be smaller than both domain half-widths")
        self.grid, self.dt, self.vmax = grid, dt, vmax
        self.hbar, self.mass = hbar, mass
        self.boundary_tolerance = boundary_tolerance
        half = dt / 2
        self.phase_v = np.exp(-1j * potential * half / (2 * hbar))
        self.phase_k = np.exp(-1j * hbar * grid.k2 * half / (2 * mass))
        self.edge_mask = ((abs(grid.x)[None, :] >= grid.Lx/2-edge_width)
                          | (abs(grid.y)[:, None] >= grid.Ly/2-edge_width))
        self.diag = {"raw_speed_max": 0.0, "clipped_stage_evaluations": 0,
                     "outside_stage_evaluations": 0, "outside_completed_positions": 0, "edge_mass_max": 0.0,
                     "edge_amplitude_max": 0.0, "norm_error_max": 0.0}

    def half_step(self, psi):
        return self.phase_v * np.fft.ifft2(self.phase_k * np.fft.fft2(self.phase_v * psi))

    def inspect_wave(self, psi):
        density = abs(psi) ** 2
        cell_area = self.grid.dx * self.grid.dy
        self.diag["edge_mass_max"] = max(self.diag["edge_mass_max"], float(np.sum(density[self.edge_mask]) * cell_area))
        self.diag["edge_amplitude_max"] = max(self.diag["edge_amplitude_max"], float(np.sqrt(np.max(density[self.edge_mask]))))
        self.diag["norm_error_max"] = max(self.diag["norm_error_max"], abs(float(np.sum(density) * cell_area)-1))

    def velocity(self, positions, fields):
        if not np.all(np.isfinite(positions)):
            raise FloatingPointError("Nonfinite trajectory stage; no fabricated velocity substituted")
        self.diag["outside_stage_evaluations"] += int(np.count_nonzero(self.grid.outside(positions)))
        psi, px, py = [self.grid.interpolate(field, positions) for field in fields]
        rho = abs(psi) ** 2
        if np.any(rho <= np.finfo(float).tiny):
            raise FloatingPointError("Trajectory stage reached an unresolved node; refine integration instead of adding to complex psi")
        velocity = (self.hbar/self.mass) * np.stack(
            (np.imag(np.conj(psi)*px)/rho, np.imag(np.conj(psi)*py)/rho), axis=1)
        if not np.all(np.isfinite(velocity)):
            raise FloatingPointError("Nonfinite guidance velocity")
        speed = np.linalg.norm(velocity, axis=1)
        self.diag["raw_speed_max"] = max(self.diag["raw_speed_max"], float(np.max(speed)))
        if self.vmax > 0:
            mask = speed > self.vmax
            self.diag["clipped_stage_evaluations"] += int(np.count_nonzero(mask))
            velocity[mask] *= (self.vmax / speed[mask])[:, None]
        return velocity

    def step(self, psi, derivatives, positions):
        mid = self.half_step(psi)
        new = self.half_step(mid)
        mid_derivatives = self.grid.derivatives(mid)
        new_derivatives = self.grid.derivatives(new)
        k1 = self.velocity(positions, (psi, *derivatives))
        k2 = self.velocity(positions+self.dt*k1/2, (mid, *mid_derivatives))
        k3 = self.velocity(positions+self.dt*k2/2, (mid, *mid_derivatives))
        k4 = self.velocity(positions+self.dt*k3, (new, *new_derivatives))
        updated = positions + self.dt * (k1+2*k2+2*k3+k4) / 6
        if not np.all(np.isfinite(updated)):
            raise FloatingPointError("Nonfinite completed trajectory step")
        self.diag["outside_completed_positions"] += int(np.count_nonzero(self.grid.outside(updated)))
        self.inspect_wave(mid)
        self.inspect_wave(new)
        return new, new_derivatives, updated

    def issues(self):
        result = []
        if self.diag["clipped_stage_evaluations"]:
            result.append("guidance_velocity_modified_by_clipping")
        if self.diag["outside_stage_evaluations"] or self.diag["outside_completed_positions"]:
            result.append("particle_stages_used_periodic_wrap")
        if self.diag["edge_mass_max"] > self.boundary_tolerance:
            result.append("boundary_layer_probability_exceeds_declared_tolerance")
        if self.diag["norm_error_max"] > 1e-9:
            result.append("wave_norm_drift")
        return result


def relative_l2(actual, reference):
    actual, reference = np.asarray(actual), np.asarray(reference)
    if actual.shape != reference.shape or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(reference)):
        raise ValueError("Distribution shapes or finite values do not match")
    denominator = np.linalg.norm(reference.ravel())
    return None if denominator == 0 else float(np.linalg.norm((actual-reference).ravel()) / denominator)


def smooth_mass(mass, grid, sigma):
    positive(sigma, "sigma_kde", allow_zero=True)
    if sigma == 0:
        return mass.copy()
    return np.fft.ifft2(np.fft.fft2(mass) * np.exp(-0.5*grid.k2*sigma**2)).real


def overlap_weights(cell_edges, bin_edges):
    """Integrate a cell-constant field over arbitrary bins, without dropping edges."""
    return np.maximum(0.0, np.minimum(bin_edges[1:, None], cell_edges[None, 1:])
                      - np.maximum(bin_edges[:-1, None], cell_edges[None, :-1]))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_hashes(entrypoint):
    return {path.name: sha256(path) for path in (Path(__file__), Path(entrypoint))}


def output_paths(output):
    path = Path(output).resolve()
    if path.suffix.lower() != ".json":
        raise ValueError("--out must be a JSON path")
    array_path = path.with_suffix(".arrays.npz")
    if path.exists() or array_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing evidence: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path, array_path


def write_result(output, result, arrays=None):
    path, array_path = output_paths(output)
    payload = {"schema": SCHEMA, "code_version": VERSION, **result}
    if arrays:
        with open(array_path, "xb") as stream:
            np.savez_compressed(stream, **arrays)
        payload["array_artifact"] = {"path": array_path.name, "sha256": sha256(array_path),
                                     "keys": sorted(arrays)}
    with open(path, "x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
    print(json.dumps({"output": str(path), "metrics": payload.get("metrics"),
                      "validity": payload.get("validity")}, ensure_ascii=False, allow_nan=False))


def common_parser(parser, particle_default):
    parser.add_argument("--Lx", type=float, default=60.0)
    parser.add_argument("--Ly", type=float, default=20.0)
    parser.add_argument("--Nx", type=int, default=1201)
    parser.add_argument("--Ny", type=int, default=401)
    parser.add_argument("--T", type=float, default=6.0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--Np", type=int, default=particle_default)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--vmax", type=float, default=0.0,
                        help="0 disables artificial speed clipping; nonzero interventions invalidate acceptance")
    parser.add_argument("--edge_width", type=float, default=1.0)
    parser.add_argument("--boundary_tolerance", type=float, default=1e-8)


def metadata(args, entrypoint):
    return {"model": MODEL, "config": {k: v for k, v in vars(args).items() if k != "out"},
            "source_sha256": source_hashes(entrypoint),
            "method": {"grid": "periodic_cell_centres", "wave": "two_Strang_half_steps",
                       "particles": "RK4_physical_half_time_wave_stages",
                       "coupled_time_order": 2, "spatial_interpolation": "bilinear_complex_wave_and_spectral_derivatives",
                       "initial_sampling": "continuous_Gaussian_conditioned_on_box"}}


def validity(issues):
    return {"eligible_for_convergence_analysis": not issues, "issues": issues,
            "quantitative_acceptance": "pending_independent_convergence_domain_and_seed_checks",
            "warning": "A single run is not proof of quantitative correctness."}
