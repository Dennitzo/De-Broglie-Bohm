"""Free, entangled two-particle Dirac benchmark on the infinite line.

Run from the project root:
    python work/relativistic_two_particle.py --tag 20261017

Only new relativistic_two_particle_* files are written. Existing output is never
overwritten. NumPy quadrature approximates the continuum Fourier integral, not
a periodic physical model. All diagnostics use the same explicitly specified
positive-energy state. The Lorentz test transforms the physical foliation too.
The numerical acceptance bounds below are fixed before the final run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import sys
import time

# Small matrix products must not start dozens of BLAS threads.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import scipy
from scipy.integrate import simpson

MASS = 1.0
P0 = 1.0
SIGMA_P = 0.6
BOOST = 0.6
GAMMA = 1 / np.sqrt(1 - BOOST**2)
RAPIDITY = np.arctanh(BOOST)
SPIN = np.diag(np.exp(np.array([-1, 1]) * RAPIDITY / 2))
LAMBDA = GAMMA * np.array([[1, -BOOST], [-BOOST, 1]])
NCOV = GAMMA * np.array([1.0, BOOST])
OVERLAP = np.exp(-P0**2 / (2 * SIGMA_P**2))
ALPHA = np.diag([1.0, -1.0])
BETA = np.array([[0.0, 1.0], [1.0, 0.0]])
OPS = np.array([np.kron(a, b) for a in (np.eye(2), ALPHA)
                for b in (np.eye(2), ALPHA)])
OP_DIAGONALS = np.diagonal(OPS, axis1=1, axis2=2)
T_FINAL = 3.2
THRESHOLDS = {
    "gram_error": 2e-10,
    "normalization_error": 2e-9,
    "purity_error": 2e-9,
    "lorentz_relative_error": 2e-12,
    "pde_relative_error": 2e-12,
    "continuity_relative_error": 2e-12,
    "momentum_velocity_error": 2e-10,
    "window_velocity_error": 2e-10,
    "trajectory_fine_error": 2e-7,
    "lorentz_trajectory_fine_error": 2e-7,
    "flow_density_jacobian_error": 3e-6,
    "causal_violation": 2e-12,
    "nonlocal_velocity_difference_minimum": 1e-2,
    "wrong_normal_velocity_difference_minimum": 1e-2,
}


class Model:
    def __init__(self, n=257, pmax=8.0):
        self.n, self.pmax = int(n), float(pmax)
        self.p = np.linspace(-pmax, pmax, n)
        self.energy = np.sqrt(self.p**2 + MASS**2)
        self.dp = self.p[1] - self.p[0]
        self.weights = np.full(n, self.dp)
        self.weights[[0, -1]] *= 0.5
        self.u = np.stack([np.sqrt((self.energy + self.p) / (2 * self.energy)),
                           np.sqrt((self.energy - self.p) / (2 * self.energy))], -1)
        prefactor = (2 * np.pi * SIGMA_P**2) ** (-0.25)
        a = prefactor * np.exp(-(self.p - P0)**2 / (4 * SIGMA_P**2))
        b = prefactor * np.exp(-(self.p + P0)**2 / (4 * SIGMA_P**2))
        self.amplitudes = np.stack([a, (b - OVERLAP * a) / np.sqrt(1 - OVERLAP**2)])
        self.kernel = self.amplitudes[:, :, None] * self.u[None, :, :]
        self.kernel *= self.weights[None, :, None] / np.sqrt(2 * np.pi)

    def packets(self, t, x, derivative=None):
        t, x = np.broadcast_arrays(t, x)
        shape = t.shape
        phase = np.exp(1j * (x.reshape(-1, 1) * self.p - t.reshape(-1, 1) * self.energy))
        if derivative == "t":
            phase *= -1j * self.energy
        elif derivative == "x":
            phase *= 1j * self.p
        return np.einsum("bp,aps->bas", phase, self.kernel, optimize=True).reshape(shape + (2, 2))

    @staticmethod
    def combine(first, second):
        return (np.einsum("...a,...b->...ab", first[..., 0, :], second[..., 0, :])
                + 1j * np.einsum("...a,...b->...ab", first[..., 1, :], second[..., 1, :])) / np.sqrt(2)

    def state(self, t1, x1, t2, x2, derivatives=False):
        first, second = self.packets(t1, x1), self.packets(t2, x2)
        psi = self.combine(first, second)
        if not derivatives:
            return psi
        return psi, (self.combine(self.packets(t1, x1, "t"), second),
                     self.combine(self.packets(t1, x1, "x"), second),
                     self.combine(first, self.packets(t2, x2, "t")),
                     self.combine(first, self.packets(t2, x2, "x")))

    def boosted_state(self, t1p, x1p, t2p, x2p, derivatives=False):
        t1, x1 = GAMMA * (t1p + BOOST * x1p), GAMMA * (x1p + BOOST * t1p)
        t2, x2 = GAMMA * (t2p + BOOST * x2p), GAMMA * (x2p + BOOST * t2p)
        value = self.state(t1, x1, t2, x2, derivatives)
        if not derivatives:
            return spin_transform(value)
        psi, (dt1, dx1, dt2, dx2) = value
        return spin_transform(psi), tuple(spin_transform(d) for d in (
            GAMMA * (dt1 + BOOST * dx1), GAMMA * (BOOST * dt1 + dx1),
            GAMMA * (dt2 + BOOST * dx2), GAMMA * (BOOST * dt2 + dx2)))

    def velocity(self, tau, xy, primed=False):
        if primed:
            times = tau / GAMMA - BOOST * xy
            psi = self.boosted_state(times[:, 0], xy[:, 0], times[:, 1], xy[:, 1])
            current = tensor_current(psi)
            q = np.einsum("a,nab,b->n", NCOV, current, NCOV)
            return np.stack([(current @ NCOV)[:, 1], (NCOV @ current)[:, 1]], -1) / q[:, None]
        current = tensor_current(self.state(tau, xy[:, 0], tau, xy[:, 1]))
        return np.stack([current[:, 1, 0], current[:, 0, 1]], -1) / current[:, 0, 0, None]


class PrimedModel(Model):
    """Independent uniform p' quadrature, not Lorentz-mapped lab quadrature nodes."""
    def __init__(self, n=833, low=-18.0, high=8.0):
        self.n, self.pmax = n, max(abs(low), abs(high))
        self.low, self.high = low, high
        self.p = np.linspace(low, high, n)
        self.energy = np.sqrt(self.p**2 + MASS**2)
        self.dp = self.p[1] - self.p[0]
        self.weights = np.full(n, self.dp)
        self.weights[[0, -1]] *= .5
        self.u = np.stack([np.sqrt((self.energy + self.p) / (2 * self.energy)),
                           np.sqrt((self.energy - self.p) / (2 * self.energy))], -1)
        inverse_p = GAMMA * (self.p + BOOST * self.energy)
        inverse_energy = GAMMA * (self.energy + BOOST * self.p)
        prefactor = (2 * np.pi * SIGMA_P**2) ** (-.25)
        jacobian_sqrt = np.sqrt(inverse_energy / self.energy)
        a = jacobian_sqrt * prefactor * np.exp(-(inverse_p - P0)**2 / (4 * SIGMA_P**2))
        b = jacobian_sqrt * prefactor * np.exp(-(inverse_p + P0)**2 / (4 * SIGMA_P**2))
        self.amplitudes = np.stack([a, (b - OVERLAP * a) / np.sqrt(1 - OVERLAP**2)])
        self.kernel = self.amplitudes[:, :, None] * self.u[None, :, :]
        self.kernel *= self.weights[None, :, None] / np.sqrt(2 * np.pi)

    def velocity(self, tau, xy, primed=True):
        times = tau / GAMMA - BOOST * xy
        psi = self.state(times[:, 0], xy[:, 0], times[:, 1], xy[:, 1])
        current = tensor_current(psi)
        q = np.einsum("a,nab,b->n", NCOV, current, NCOV)
        return np.stack([(current @ NCOV)[:, 1], (NCOV @ current)[:, 1]], -1) / q[:, None]


def spin_transform(psi):
    return np.einsum("ab,...bc,dc->...ad", SPIN, psi, SPIN)


def tensor_current(psi):
    squared = np.abs(psi.reshape(psi.shape[:-2] + (4,)))**2
    return np.einsum("...s,ms->...m", squared, OP_DIAGONALS).reshape(psi.shape[:-2] + (2, 2))


def tensor_derivative(psi, dpsi):
    products = 2 * np.real(psi.conj() * dpsi).reshape(psi.shape[:-2] + (4,))
    return np.einsum("...s,ms->...m", products, OP_DIAGONALS).reshape(psi.shape[:-2] + (2, 2))


def rk4(model, initial, dt, primed=False):
    steps = int(round(T_FINAL / dt))
    assert abs(steps * dt - T_FINAL) < 1e-13
    history = np.empty((steps + 1,) + initial.shape)
    history[0] = initial
    for k in range(steps):
        t, y = k * dt, history[k]
        k1 = model.velocity(t, y, primed)
        k2 = model.velocity(t + dt / 2, y + dt * k1 / 2, primed)
        k3 = model.velocity(t + dt / 2, y + dt * k2 / 2, primed)
        k4 = model.velocity(t + dt, y + dt * k3, primed)
        history[k + 1] = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return np.arange(steps + 1) * dt, history


def spatial_gram(model, t, extent, dx):
    x = np.linspace(-extent, extent, int(round(2 * extent / dx)) + 1)
    phi = model.packets(t, x)
    gram = simpson(np.einsum("xas,xbs->xab", phi.conj(), phi), x=x, axis=0)
    coeff = np.diag([1, 1j]) / np.sqrt(2)
    norm = np.real(np.einsum("ab,cd,ac,bd->", coeff.conj(), coeff, gram, gram))
    vals, vecs = np.linalg.eigh(gram)
    root = (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.conj().T
    matrix = root @ coeff @ root.T
    reduced = matrix @ matrix.conj().T
    reduced /= np.trace(reduced)
    return {"time": t, "extent": extent, "dx": dx, "two_particle_norm": float(norm),
            "gram_error": float(np.max(abs(gram - np.eye(2)))),
            "purity": float(np.trace(reduced @ reduced).real)}, gram


def relative_max(residual, scale):
    return float(np.max(abs(residual)) / max(float(np.max(abs(scale))), 1e-300))


def algebra_checks(model, primed_model, arrays):
    grid = np.linspace(-3, 3, 9)
    t1, x1, x2 = np.meshgrid(np.array([0, .4, 1.6, 3.2]), grid, grid, indexing="ij")
    t1, x1, x2 = (a.ravel() for a in (t1, x1, x2))
    # Unequal times test the multi-time equations away from the lab leaves too.
    t2 = t1 + .17 * np.sin(x1 - x2)
    psi, deriv = model.state(t1, x1, t2, x2, True)
    current = tensor_current(psi)
    t1p, x1p = GAMMA * (t1 - BOOST * x1), GAMMA * (x1 - BOOST * t1)
    t2p, x2p = GAMMA * (t2 - BOOST * x2), GAMMA * (x2 - BOOST * t2)
    psip, derivp = primed_model.state(t1p, x1p, t2p, x2p, True)
    pulled_back = model.boosted_state(t1p, x1p, t2p, x2p)
    currentp = tensor_current(psip)
    expected = np.einsum("ac,ncd,bd->nab", LAMBDA, current, LAMBDA)
    qprime = np.einsum("a,nab,b->n", NCOV, currentp, NCOV)
    j1prime, j2prime = currentp @ NCOV, NCOV @ currentp
    j1lab, j2lab = current[:, :, 0], current[:, 0, :]
    pde_errors, cont_errors = [], []
    for wave, ds in ((psi, deriv), (psip, derivp)):
        dt1, dx1, dt2, dx2 = ds
        residual1 = 1j * dt1 + 1j * np.einsum("ab,nbc->nac", ALPHA, dx1) - np.einsum("ab,nbc->nac", BETA, wave)
        residual2 = 1j * dt2 + 1j * np.einsum("nab,cb->nac", dx2, ALPHA) - np.einsum("nab,cb->nac", wave, BETA)
        pde_errors.extend([relative_max(residual1, wave), relative_max(residual2, wave)])
        djt1, djx1, djt2, djx2 = (tensor_derivative(wave, d) for d in ds)
        cont_errors.extend([relative_max(djt1[:, 0, :] + djx1[:, 1, :], tensor_current(wave)),
                            relative_max(djt2[:, :, 0] + djx2[:, :, 1], tensor_current(wave))])
    # The negative-control pairs really belong to one common lab leaf. Holding
    # n'=(1,0) at these same events is an algebraic misuse, not an alternative
    # complete t'=constant trajectory calculation.
    negative_lab = tensor_current(model.state(t1, x1, t1, x2))
    t2eqp, x2eqp = GAMMA * (t1 - BOOST * x2), GAMMA * (x2 - BOOST * t1)
    negative_prime = tensor_current(primed_model.state(t1p, x1p, t2eqp, x2eqp))
    negative_j1 = negative_prime @ NCOV
    wrong1 = negative_prime[:, 1, 0] / negative_prime[:, 0, 0]
    true1 = negative_j1[:, 1] / negative_j1[:, 0]
    density_mask = negative_lab[:, 0, 0] > 1e-4
    negative_index = np.where(density_mask)[0][np.argmax(abs(wrong1[density_mask] - true1[density_mask]))]
    arrays.update(algebra_events=np.column_stack([t1, x1, t2, x2]), algebra_J=current,
                  algebra_Jprime=currentp, algebra_density=qprime,
                  negative_events=np.column_stack([t1, x1, t1, x2]),
                  algebra_wrong_normal_v1=wrong1, algebra_true_normal_v1=true1)
    jlist = [j1lab, j2lab, j1prime, j2prime]
    # Future timelike normals with rapidity in [-1.5,1.5], not only this boost.
    minimum_q, normal_causal = 1.0, 0.0
    for eta in np.linspace(-1.5, 1.5, 13):
        n = np.array([np.cosh(eta), -np.sinh(eta)])
        q = np.einsum("a,nab,b->n", n, current, n)
        minimum_q = min(minimum_q, float(q.min()))
        for j in (current @ n, n @ current):
            normal_causal = max(normal_causal, float(np.max(np.maximum(abs(j[:, 1]) - j[:, 0], 0))))
    return {
        "points": len(t1), "unequal_time_offset_max": float(np.max(abs(t1 - t2))),
        "tensor_covariance_relative_error": relative_max(currentp - expected, expected),
        "independent_primed_spinor_relative_error": relative_max(psip - pulled_back, pulled_back),
        "density_scalar_relative_error": relative_max(qprime - current[:, 0, 0], current[:, 0, 0]),
        "coordinate_density_measure_relative_error": relative_max((qprime / GAMMA**2) * GAMMA**2 - current[:, 0, 0], current[:, 0, 0]),
        "j1_covariance_relative_error": relative_max(j1prime - j1lab @ LAMBDA.T, j1prime),
        "j2_covariance_relative_error": relative_max(j2prime - j2lab @ LAMBDA.T, j2prime),
        "dirac_residual_relative_max": max(pde_errors),
        "continuity_residual_relative_max": max(cont_errors),
        "causal_violation_absolute_max": max(float(np.max(np.maximum(abs(j[:, 1]) - j[:, 0], 0))) for j in jlist),
        "minimum_sampled_density": float(current[:, 0, 0].min()),
        "minimum_sampled_normal_density": minimum_q,
        "normal_rapidity_grid": np.linspace(-1.5, 1.5, 13).tolist(),
        "all_normals_causal_violation": normal_causal,
        "wrong_normal_velocity_difference_max_density_gt_1e_4": float(abs(wrong1[negative_index] - true1[negative_index])),
        "wrong_normal_witness": {"lab_events": [float(a[negative_index]) for a in (t1, x1, t1, x2)],
                                 "density": float(negative_lab[negative_index, 0, 0]),
                                 "wrong_coordinate_velocity": float(wrong1[negative_index]),
                                 "correct_coordinate_velocity": float(true1[negative_index])},
    }


def nonlocal_check(model, arrays):
    remote = np.linspace(-3, 3, 241)
    trials = []
    for tau in (.0, .4, 1.6):
        for local in (-1.0, 0.0, 1.0):
            xy = np.column_stack([np.full_like(remote, local), remote])
            current = tensor_current(model.state(tau, xy[:, 0], tau, xy[:, 1]))
            velocity = model.velocity(tau, xy)[:, 0]
            mask = current[:, 0, 0] > 1e-3
            low, high = np.where(mask)[0][np.argmin(velocity[mask])], np.where(mask)[0][np.argmax(velocity[mask])]
            trials.append((float(velocity[high] - velocity[low]), tau, local, low, high, velocity, current[:, 0, 0]))
    spread, tau, local, low, high, velocity, q = max(trials, key=lambda row: row[0])
    arrays.update(nonlocal_remote=remote, nonlocal_v1=velocity, nonlocal_q=q)
    return {"tau": tau, "x1": local, "x2_pair": remote[[low, high]].tolist(),
            "v1_pair": velocity[[low, high]].tolist(), "density_pair": q[[low, high]].tolist(),
            "velocity_difference": spread,
            "interpretation": "Conditional configuration dependence; no operational signalling test or protocol."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="20261017")
    parser.add_argument("--pilot", action="store_true", help="Print diagnostic pilot without creating output files")
    args = parser.parse_args()
    if not args.tag.replace("_", "").replace("-", "").isalnum():
        parser.error("tag must be alphanumeric with underscores/hyphens")
    root = Path(__file__).resolve().parents[1]
    stem = f"relativistic_two_particle_{args.tag}"
    outputs = [root / "work" / (stem + ext) for ext in (".json", ".npz", ".md")]
    figures = [root / "figures" / (stem + ext) for ext in (".png", ".svg")]
    if not args.pilot and any(p.exists() for p in outputs + figures):
        parser.error("Output already exists; use a new --tag. Historical outputs are never overwritten.")
    start = time.perf_counter()
    model, primed_model, arrays = Model(), PrimedModel(), {}
    metrics = {
        "model": {"dimension": "1+1", "particles": "two distinguishable, free, equal mass",
                  "hbar": 1, "c": 1, "mass": MASS, "p0": P0, "sigma_p": SIGMA_P,
                  "momentum_packet_overlap": float(OVERLAP), "boost_v": BOOST,
                  "boost_gamma": float(GAMMA), "T": T_FINAL, "momentum_N": model.n,
                  "momentum_half_width": model.pmax, "momentum_dp": float(model.dp),
                  "primed_momentum_N": primed_model.n, "primed_momentum_interval": [primed_model.low, primed_model.high],
                  "primed_momentum_dp": float(primed_model.dp),
                  "state": "(phi0 tensor phi0 + i phi1 tensor phi1)/sqrt(2)",
                  "phi0": "A", "phi1": "(B-r*A)/sqrt(1-r^2)",
                  "domain": "infinite spatial line; finite quadrature and integration domains are convergence controls"},
        "acceptance_thresholds": THRESHOLDS,
        "algebra": algebra_checks(model, primed_model, arrays),
        "nonlocal_dependence": nonlocal_check(model, arrays),
    }
    print("Algebra and nonlocality", json.dumps(metrics["algebra"], indent=2), flush=True)
    norms = []
    for extent in (8.0, 12.0, 16.0):
        for tau in (0.0, T_FINAL):
            row, gram = spatial_gram(model, tau, extent, .025)
            norms.append(row)
            arrays[f"gram_L{int(extent)}_t{tau:g}"] = gram
    for dx in (.1, .05):
        row, gram = spatial_gram(model, T_FINAL, 16.0, dx)
        norms.append(row)
    metrics["spatial_norm_controls"] = norms
    spectral_gram = (model.amplitudes * model.weights) @ model.amplitudes.T
    metrics["momentum_gram"] = spectral_gram.tolist()
    metrics["momentum_gram_error"] = float(np.max(abs(spectral_gram - np.eye(2))))
    totals = []
    for spectral_model in (model, primed_model):
        occupation = np.sum(spectral_model.amplitudes**2, axis=0)
        totals.append(np.array([np.sum(spectral_model.weights * occupation * spectral_model.energy),
                                np.sum(spectral_model.weights * occupation * spectral_model.p)]))
    metrics["total_mean_four_momentum"] = {
        "lab": totals[0].tolist(), "prime_independent": totals[1].tolist(),
        "prime_covariance_error": float(np.max(abs(totals[1] - LAMBDA @ totals[0]))),
        "normal_from_total_momentum_lab": (totals[0] / np.sqrt(totals[0][0]**2 - totals[0][1]**2)).tolist(),
        "normal_from_total_momentum_prime": (totals[1] / np.sqrt(totals[1][0]**2 - totals[1][1]**2)).tolist(),
    }
    if args.pilot:
        print(json.dumps(metrics, indent=2), flush=True)
        return

    initial_axis = np.array([-1.8, -.9, 0, .9, 1.8])
    initial = np.array(np.meshgrid(initial_axis, initial_axis)).reshape(2, -1).T
    arrays["trajectory_initial"] = initial
    # Distinct time meshes prevent identical RK4 affine covariance from hiding
    # time-discretization error: prime uses dt/2, lab dt.
    tref, reference = rk4(Model(513), initial, .005)
    arrays.update(trajectory_reference_t=tref, trajectory_reference_lab=reference)
    convergence = []
    finest_lab, finest_time = None, None
    for dt in (.08, .04, .02):
        ts, lab = rk4(model, initial, dt)
        _, prime = rk4(primed_model, GAMMA * initial, dt / 2, True)
        expected = GAMMA * (lab - BOOST * ts[:, None, None])
        transformed_ref = reference[::int(round(dt / .005))]
        row = {"dt_lab": dt, "dt_prime": dt / 2,
               "lab_error_vs_dt0005_N513": float(np.max(abs(lab - transformed_ref))),
               "lorentz_worldline_error": float(np.max(abs(prime[::2] - expected))),
               "lab_endpoint_error": float(np.max(abs(lab[-1] - reference[-1])))}
        convergence.append(row)
        arrays[f"trajectory_lab_dt{dt:g}"] = lab
        arrays[f"trajectory_prime_dt{dt / 2:g}"] = prime
        finest_lab, finest_time = lab, ts
        print("Time convergence", row, flush=True)
    metrics["time_convergence"] = convergence
    metrics["observed_time_orders"] = [float(np.log2(convergence[k]["lab_error_vs_dt0005_N513"] / convergence[k+1]["lab_error_vs_dt0005_N513"])) for k in (0, 1)]

    # Independently refine Fourier quadrature at the same physical sample points
    # and with the same time integrator (spatial/temporal errors kept separate).
    reference_velocity = Model(513).velocity(1.6, initial)
    _, pref = rk4(Model(513), initial, .04)
    quadrature = []
    for n in (33, 65, 129, 257):
        candidate = Model(n)
        _, history = rk4(candidate, initial, .04)
        row = {"N": n, "dp": float(candidate.dp),
               "velocity_error_vs_N513": float(np.max(abs(candidate.velocity(1.6, initial) - reference_velocity))),
               "trajectory_error_vs_N513_same_dt": float(np.max(abs(history - pref)))}
        quadrature.append(row)
        arrays[f"trajectory_momentum_N{n}"] = history
        print("Momentum convergence", row, flush=True)
    metrics["momentum_convergence"] = quadrature
    window_reference = Model(321, 10).velocity(1.6, initial)
    metrics["momentum_window_controls"] = [
        {"pmax": pmax, "N": int(32 * pmax + 1),
         "velocity_error_vs_pmax10": float(np.max(abs(Model(int(32 * pmax + 1), pmax).velocity(1.6, initial) - window_reference)))}
        for pmax in (6.0, 8.0, 10.0)]
    prime_xy = GAMMA * (initial - BOOST * 1.6)
    pvref = PrimedModel(1665).velocity(1.6, prime_xy)
    metrics["primed_momentum_convergence"] = [
        {"N": n, "dp": 26 / (n - 1),
         "velocity_error_vs_N1665": float(np.max(abs(PrimedModel(n).velocity(1.6, prime_xy) - pvref)))}
        for n in (209, 417, 833)]
    # Independent tilted-leaf quadrature with dx'/gamma and alpha.n'.
    primed_grams = []
    for tau in (0.0, T_FINAL):
        xlab = np.linspace(-16, 16, 1281)
        xp = GAMMA * (xlab - BOOST * tau)
        tp = tau / GAMMA - BOOST * xp
        phip = primed_model.packets(tp, xp)
        an = NCOV[0] * np.eye(2) + NCOV[1] * ALPHA
        gram = simpson(np.einsum("xas,st,xbt->xab", phip.conj(), an, phip), x=xp, axis=0) / GAMMA
        coeff = np.diag([1, 1j]) / np.sqrt(2)
        norm = float(np.real(np.einsum("ab,cd,ac,bd->", coeff.conj(), coeff, gram, gram)))
        primed_grams.append({"tau": tau, "lab_equivalent_extent": 16.0, "dx_prime": .025 * GAMMA,
                            "coordinate_density": "Q_prime=q_Sigma/gamma^2", "two_particle_norm": norm,
                            "gram_error": float(np.max(abs(gram - np.eye(2))))})
        arrays[f"primed_tilted_gram_t{tau:g}"] = gram
    metrics["primed_tilted_leaf_norm_controls"] = primed_grams

    # A deterministic local equivariance test: q(T,F_T(x))*det(DF_T)=q(0,x).
    # DF is independently estimated with centered perturbations of initial data.
    flow_initial = np.array([[-1.2, -.6], [-.7, 1.1], [0, 0], [.5, -.9], [1.2, .8]])
    q0 = tensor_current(model.state(0, flow_initial[:, 0], 0, flow_initial[:, 1]))[:, 0, 0]
    flow_checks = []
    for epsilon in (.0005, .00025, .000125):
        perturbed = [flow_initial]
        for axis in (0, 1):
            delta = np.zeros_like(flow_initial)
            delta[:, axis] = epsilon
            perturbed.extend([flow_initial + delta, flow_initial - delta])
        _, flow = rk4(model, np.concatenate(perturbed), .01)
        f = flow[-1].reshape(5, len(flow_initial), 2)
        jac = np.stack([(f[1] - f[2]) / (2 * epsilon), (f[3] - f[4]) / (2 * epsilon)], axis=-1)
        det = np.linalg.det(jac)
        qf = tensor_current(model.state(T_FINAL, f[0, :, 0], T_FINAL, f[0, :, 1]))[:, 0, 0]
        ratio = qf * det / q0
        flow_checks.append({"epsilon": epsilon, "dt": .01, "minimum_determinant": float(det.min()),
                            "density_jacobian_relative_error": float(np.max(abs(ratio - 1)))})
        arrays[f"flow_jacobian_eps{epsilon:g}"] = jac
        arrays[f"flow_density_ratio_eps{epsilon:g}"] = ratio
        arrays[f"flow_endpoint_eps{epsilon:g}"] = f[0]
        print("Flow equivariance", flow_checks[-1], flush=True)
    arrays["flow_initial"] = flow_initial
    metrics["local_flow_equivariance"] = flow_checks
    eps_prime = .000125 * GAMMA
    primed_initial = GAMMA * flow_initial
    perturbed_prime = [primed_initial]
    for axis in (0, 1):
        delta = np.zeros_like(primed_initial)
        delta[:, axis] = eps_prime
        perturbed_prime.extend([primed_initial + delta, primed_initial - delta])
    _, prime_flow = rk4(primed_model, np.concatenate(perturbed_prime), .005, True)
    pf = prime_flow[-1].reshape(5, len(flow_initial), 2)
    pjac = np.stack([(pf[1] - pf[2]) / (2 * eps_prime), (pf[3] - pf[4]) / (2 * eps_prime)], -1)
    def prime_coordinate_q(tau, xy):
        tp = tau / GAMMA - BOOST * xy
        jp = tensor_current(primed_model.state(tp[:, 0], xy[:, 0], tp[:, 1], xy[:, 1]))
        return np.einsum("a,nab,b->n", NCOV, jp, NCOV) / GAMMA**2
    pq0, pqf = prime_coordinate_q(0, primed_initial), prime_coordinate_q(T_FINAL, pf[0])
    pratio = pqf * np.linalg.det(pjac) / pq0
    metrics["primed_local_flow_equivariance"] = {
        "epsilon_prime": eps_prime, "dt": .005,
        "coordinate_density": "Q_prime=q_Sigma/gamma^2",
        "density_jacobian_relative_error": float(np.max(abs(pratio - 1))),
        "lab_primed_flow_jacobian_max_difference": float(np.max(abs(pjac - jac))),
    }
    arrays.update(prime_flow_initial=primed_initial, prime_flow_jacobian=pjac,
                  prime_flow_coordinate_q0=pq0, prime_flow_coordinate_qf=pqf,
                  prime_flow_density_ratio=pratio)
    a, th = metrics["algebra"], THRESHOLDS
    final_norms = [r for r in norms if r["extent"] == 16 and r["dx"] == .025]
    checks = {
        "spectral_orthonormality": metrics["momentum_gram_error"] < th["gram_error"],
        "finite_domain_norm": max(abs(r["two_particle_norm"] - 1) for r in final_norms) < th["normalization_error"],
        "reduced_purity": max(abs(r["purity"] - .5) for r in final_norms) < th["purity_error"],
        "lorentz_tensor_and_density": max(a[k] for k in ("tensor_covariance_relative_error", "density_scalar_relative_error", "j1_covariance_relative_error", "j2_covariance_relative_error")) < th["lorentz_relative_error"],
        "independent_prime_fourier": a["independent_primed_spinor_relative_error"] < th["lorentz_relative_error"],
        "prime_momentum_quadrature": metrics["primed_momentum_convergence"][-1]["velocity_error_vs_N1665"] < th["momentum_velocity_error"],
        "prime_leaf_measure_norm": max(abs(r["two_particle_norm"] - 1) for r in primed_grams) < th["normalization_error"],
        "mean_momentum_covariance": metrics["total_mean_four_momentum"]["prime_covariance_error"] < th["lorentz_relative_error"],
        "dirac_PDE": a["dirac_residual_relative_max"] < th["pde_relative_error"],
        "multi_time_continuity": a["continuity_residual_relative_max"] < th["continuity_relative_error"],
        "future_causality": max(a["causal_violation_absolute_max"], a["all_normals_causal_violation"]) < th["causal_violation"],
        "positive_sampled_density": a["minimum_sampled_normal_density"] > 0,
        "momentum_quadrature": quadrature[-1]["velocity_error_vs_N513"] < th["momentum_velocity_error"],
        "momentum_window": metrics["momentum_window_controls"][1]["velocity_error_vs_pmax10"] < th["window_velocity_error"],
        "trajectory_accuracy": convergence[-1]["lab_error_vs_dt0005_N513"] < th["trajectory_fine_error"],
        "trajectory_covariance": convergence[-1]["lorentz_worldline_error"] < th["lorentz_trajectory_fine_error"],
        "local_equivariance": flow_checks[-1]["density_jacobian_relative_error"] < th["flow_density_jacobian_error"],
        "primed_local_equivariance": metrics["primed_local_flow_equivariance"]["density_jacobian_relative_error"] < th["flow_density_jacobian_error"],
        "nonlocal_dependence": metrics["nonlocal_dependence"]["velocity_difference"] > th["nonlocal_velocity_difference_minimum"],
        "wrong_normal_negative_control": a["wrong_normal_velocity_difference_max_density_gt_1e_4"] > th["wrong_normal_velocity_difference_minimum"],
    }
    metrics["acceptance_checks"] = checks
    metrics["all_checks_pass"] = all(checks.values())
    metrics["provenance"] = {"python": sys.version, "platform": platform.platform(),
                              "numpy": np.__version__, "scipy": scipy.__version__,
                              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                              "command": f"python work/relativistic_two_particle.py --tag {args.tag}",
                              "elapsed_seconds": time.perf_counter() - start,
                              "no_experimental_data": True}
    np.savez_compressed(outputs[1], **arrays)
    metrics["provenance"]["arrays_sha256"] = hashlib.sha256(outputs[1].read_bytes()).hexdigest()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for i in (0, 4, 6, 12, 18, 20, 24):
        axes[0, 0].plot(finest_time, finest_lab[:, i, 0], lw=1.4)
    axes[0, 0].set(xlabel=r"Foliationszeit $\tau$", ylabel=r"$x_1(\tau)$", title="Verschränkte Dirac-Führung")
    axes[0, 1].plot(arrays["nonlocal_remote"], arrays["nonlocal_v1"])
    nl = metrics["nonlocal_dependence"]
    axes[0, 1].set(xlabel=r"Ferne Koordinate $x_2$", ylabel=r"$v_1$", title=f"Bedingte Geschwindigkeit: τ={nl['tau']:g}, x₁={nl['x1']:g}")
    ds = np.array([r["dt_lab"] for r in convergence])
    axes[1, 0].loglog(ds, [r["lab_error_vs_dt0005_N513"] for r in convergence], "o-", label="Zeitintegration")
    axes[1, 0].loglog(ds, [r["lorentz_worldline_error"] for r in convergence], "s-", label="Lorentzvergleich (verschiedene dt)")
    axes[1, 0].set(xlabel="Laborschritt dt", ylabel="Maximaler Ortsfehler", title="Unabhängige Zeitverfeinerung")
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].loglog([r["epsilon"] for r in flow_checks], [r["density_jacobian_relative_error"] for r in flow_checks], "o-")
    axes[1, 1].set(xlabel="Anfangsdifferenz ε", ylabel=r"$\max|q_T\det DF_T/q_0-1|$", title="Lokale Äquivarianz des Flusses")
    for ax in axes.flat:
        ax.grid(alpha=.25)
    fig.savefig(figures[0], dpi=170)
    fig.savefig(figures[1], metadata={"Date": None})
    plt.close(fig)
    metrics["provenance"]["matplotlib"] = matplotlib.__version__
    outputs[0].write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report = f"""# Verschraenktes relativistisches Zweiteilchenmodell

Reproduktion: `{metrics['provenance']['command']}` (neuen Tag waehlen).
Quell-SHA256: `{metrics['provenance']['source_sha256']}`.
Roharrays-SHA256: `{metrics['provenance']['arrays_sha256']}`.
Alle {len(checks)} Akzeptanzpruefungen bestanden: **{all(checks.values())}**.

## Modell und Einordnung

Freie, unterscheidbare massive Dirac-Teilchen in 1+1 Dimensionen, hbar=c=m=1,
alpha=sigma_z, beta=sigma_x. Positive-Energie-Gausspakete mit p0=1 und
sigma_p=0.6 werden analytisch orthonormalisiert. Der Zustand
(phi0 tensor phi0 + i phi1 tensor phi1)/sqrt(2) hat exakt zwei gleiche
Schmidtwerte und reduzierte Reinheit 1/2. Es werden keine Messdaten simuliert
oder als experimentelle Daten ausgegeben.

Die Fourierquadratur approximiert die unendliche Linie; sie definiert keine
periodischen physikalischen Randbedingungen. N, Impulsfenster, Raumfenster und
Zeitschritt werden getrennt variiert. Die Fourier-Summen haben numerische
Aliasperioden; nur hinreichend fein aufgeloeste endliche Auswertebereiche werden
fuer Normkontrollen verwendet. Die Zeitintegration erfolgt auf 25 festgelegten
Startkonfigurationen bis tau=3.2. Die Lorentzrechnung wertet ein separates
Fourierintegral auf einem unabhaengigen gleichmaessigen p'-Gitter aus,
bei mittransformierter Foliation und halbiertem Primed-Zeitschritt.
Dabei gilt f'(p')=sqrt(E/E') f(gamma*(p'+v*E')), E=gamma*(E'+v*p').
Auf den geneigten Blaettern dSigma'=dx'/gamma gilt fuer die Koordinatendichte
Q'=q_Sigma/gamma^2. Beide Groessen werden im Norm- und Flussvergleich benutzt.
Koordinatenkovarianz ist keine Unabhaengigkeit
von der Wahl der physikalischen Foliation.

Die PDE- und Kontinuitaetsresiduen stammen von analytischen Ableitungen unter
dem Fourierintegral. Die lokale numerische Aequivarianz q(T)det(DF)=q(0) wird
an fuenf positiven Dichtepunkten mit separaten Anfangsdifferenzen getestet;
dies ist keine globale numerische Existenz- oder Knotenausschlussbehauptung.
Die Nichtlokalitaetsprobe variiert x2 bei festem tau und x1 und misst die
bedingte Geschwindigkeit. Daraus folgt kein Kommunikationsprotokoll.
Die negative Kontrolle verwendet am selben Ereignispaar einer Laborzeitflaeche
faelschlich die unveraenderte Normale n'=(1,0) und ergibt eine andere
Stromrichtung. Das ist ein algebraischer Fehlanwendungstest, keine vollstaendige
alternative Foliationsdynamik: die Ereignisse sind allgemein nicht t'-gleichzeitig.

## Hauptergebnisse

- Lorentz-Tensorfehler: {a['tensor_covariance_relative_error']:.6e}.
- Skalar-Dichtefehler: {a['density_scalar_relative_error']:.6e}.
- Dirac-PDE-Residuum: {a['dirac_residual_relative_max']:.6e}.
- Mehrzeit-Kontinuitaetsresiduum: {a['continuity_residual_relative_max']:.6e}.
- Feinster Trajektorienfehler: {convergence[-1]['lab_error_vs_dt0005_N513']:.6e}.
- Feinster unabhaengiger Lorentz-Weltlinienfehler: {convergence[-1]['lorentz_worldline_error']:.6e}.
- Lokaler Dichte-Jakobian-Fehler: {flow_checks[-1]['density_jacobian_relative_error']:.6e}.
- Bedingte Geschwindigkeitsdifferenz: {nl['velocity_difference']:.6e}.
- Falsche-Normale-Kontrolle: {a['wrong_normal_velocity_difference_max_density_gt_1e_4']:.6e}.

Vollstaendige Tabellen, feste Toleranzen, Kontrollparameter und Versionen:
`{outputs[0].name}`. Roharrays: `{outputs[1].name}`.
"""
    outputs[2].write_text(report, encoding="utf-8")
    print(json.dumps({"all_checks_pass": all(checks.values()), "checks": checks,
                      "files": [str(p.relative_to(root)) for p in outputs + figures]}, indent=2), flush=True)
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
