"""Massive scalar on a static cubic three-torus: two continuum calculations.

The subtraction is the local Minkowski vacuum, not the entire torus vacuum.
No instantaneous time-dependent Casimir source or gravity solution is assumed.
"""
from __future__ import annotations

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import kv

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(maxsize=16)
def degeneracies(cutoff: int) -> tuple[np.ndarray, np.ndarray]:
    """Exact integer multiplicities for |n_i| <= cutoff, excluding n=0."""
    n = np.arange(cutoff + 1, dtype=np.int64)
    n2 = n * n
    weight = np.where(n == 0, 1, 2)
    squared = (n2[:, None, None] + n2[None, :, None] + n2[None, None, :]).ravel()
    counts = (weight[:, None, None] * weight[None, :, None] * weight[None, None, :]).ravel()
    multiplicity = np.bincount(squared, weights=counts).astype(np.int64)
    present = np.flatnonzero(multiplicity)[1:]
    assert int(multiplicity.sum()) == (2 * cutoff + 1) ** 3
    return np.sqrt(present.astype(float)), multiplicity[present]


def bessel_sum(z: float, cutoff: int) -> dict[str, float]:
    """Return L^4 rho and L^4 p, where z=mL>0, with absolute tail bounds."""
    if z <= 0 or cutoff < 1:
        raise ValueError("Require z>0 and cutoff>=1")
    r, multiplicity = degeneracies(cutoff)
    s2 = float(np.dot(multiplicity, kv(2, z * r) / r**2))
    s1 = float(np.dot(multiplicity, kv(1, z * r) / r))
    rho = -z * z * s2 / (4 * math.pi**2)
    pressure = -z * z * (s2 + z * s1) / (12 * math.pi**2)
    # |n|_infty=s has 24s^2+2 points. K1 <= K2, and
    # K2(x) <= sqrt(2pi/x)*exp(-x+2/x) from its positive integral.
    start = cutoff + 1
    q = math.exp(-z)
    pref = 26 * math.sqrt(2 * math.pi / (z * start)) * math.exp(2 / (z * start) - z * start)
    s2_tail = pref / (1 - q)
    pressure_sum_tail = pref * ((1 + z * start) / (1 - q) + z * q / (1 - q)**2)
    return {
        "rho_L4": rho,
        "pressure_L4": pressure,
        "nec_L4": rho + pressure,
        "rho_tail_bound_L4": z*z*s2_tail/(4*math.pi**2),
        "pressure_tail_bound_L4": z*z*pressure_sum_tail/(12*math.pi**2),
        "cutoff": cutoff,
        "mode_count": (2*cutoff+1)**3-1,
    }


def theta_images(t: float) -> tuple[float, float]:
    """Sum exp(-|n|²/4t), n!=0, and sum |n|² exp(-|n|²/4t).

    Direct and modular representations avoid cancellation and excessive sums.
    This calculation never calls a Bessel function.
    """
    eps = 1e-18
    if t <= 0.2:
        nmax = max(1, math.ceil(math.sqrt(-4*t*math.log(eps))))
        n = np.arange(1, nmax+1, dtype=float)
        terms = np.exp(-n*n/(4*t))
        part = float(terms.sum())
        theta = 1 + 2*part
        images = 6*part + 12*part*part + 8*part**3
        moment = 2*float(np.dot(n*n, terms))
    else:
        nmax = max(1, math.ceil(math.sqrt(-math.log(eps)/(4*math.pi**2*t))))
        n = np.arange(1, nmax+1, dtype=float)
        terms = np.exp(-4*math.pi**2*t*n*n)
        root = math.sqrt(4*math.pi*t)
        theta = root * (1 + 2*float(terms.sum()))
        images = theta**3 - 1
        moment = 2*t*theta - 32*math.pi**2*t*t*root*float(np.dot(n*n, terms))
    return images, 3*theta*theta*moment


def heat_kernel(z: float) -> dict[str, float]:
    """Independent integral, dimensionless t=s/L² and logarithmic quadrature.

    Integrating to +/-infinity would under/overflow; the chosen endpoints have
    suppression e^-5506 or e^-120. Reported quadrature errors are numerical,
    distinct from the analytic lattice-tail bounds in bessel_sum.
    """
    lo = -10.0
    hi = math.log(120/(z*z))

    def integrand(log_t: float, moment: bool) -> float:
        t = math.exp(log_t)
        images, square_sum = theta_images(t)
        if moment:
            return math.exp(-z*z*t) * square_sum / t**3
        return math.exp(-z*z*t) * images / t**2

    density, err_density = quad(lambda x: integrand(x, False), lo, hi, epsabs=2e-11, epsrel=2e-12, limit=300)
    nec_int, err_nec = quad(lambda x: integrand(x, True), lo, hi, epsabs=2e-10, epsrel=2e-12, limit=300)
    rho = -density/(32*math.pi**2)
    nec = -nec_int/(192*math.pi**2)
    return {
        "rho_L4": rho,
        "pressure_L4": nec-rho,
        "nec_L4": nec,
        "rho_quadrature_error_L4": err_density/(32*math.pi**2),
        "nec_quadrature_error_L4": err_nec/(192*math.pi**2),
        "log_t_interval": [lo, hi],
    }


def main() -> None:
    checks = []
    points = []
    for z in [0.5, 1.0, 2.0, 4.0, 8.0, 12.0]:
        cutoff = max(8, math.ceil(42/z))
        series = bessel_sum(z, cutoff)
        integral = heat_kernel(z)
        delta_rho = abs(series["rho_L4"]-integral["rho_L4"])
        delta_p = abs(series["pressure_L4"]-integral["pressure_L4"])
        checks.extend([
            {"name": f"rho_independent_z{z}", "passed": delta_rho < 3e-11},
            {"name": f"pressure_independent_z{z}", "passed": delta_p < 3e-11},
            {"name": f"negative_nec_z{z}", "passed": series["nec_L4"] < 0},
            {"name": f"tails_z{z}", "passed": max(series["rho_tail_bound_L4"], series["pressure_tail_bound_L4"]) < 1e-13},
        ])
        # Independent thermodynamic differentiation at fixed m=1: L=z.
        step = z*1e-4
        def rho_at_L(L: float) -> float:
            return bessel_sum(L, cutoff)["rho_L4"]/L**4
        derivative = (rho_at_L(z-2*step)-8*rho_at_L(z-step)+8*rho_at_L(z+step)-rho_at_L(z+2*step))/(12*step)
        pressure_fd = (-rho_at_L(z)-z*derivative/3)*z**4
        checks.append({"name": f"thermodynamic_pressure_z{z}", "passed": abs(pressure_fd-series["pressure_L4"]) < 1e-10})
        points.append({"mL":z, "bessel":series, "heat_kernel":integral,
                       "rho_method_difference_L4":delta_rho, "pressure_method_difference_L4":delta_p,
                       "pressure_five_point_difference_L4":abs(pressure_fd-series["pressure_L4"])})

    zgrid = np.geomspace(0.5, 12, 90)
    grid = np.array([[row[key] for key in ["rho_L4","pressure_L4","nec_L4","rho_tail_bound_L4","pressure_tail_bound_L4"]]
                    for z in zgrid for row in [bessel_sum(float(z),max(8,math.ceil(42/z)))]] )
    cutoff_grid = np.arange(1, 65)
    variance = (2*cutoff_grid+1)**3 * (1+cutoff_grid*(cutoff_grid+1)) / 8
    direct_variance = []
    for cutoff in [1,2,4,8]:
        n = np.arange(-cutoff,cutoff+1)
        direct = float(np.sum(1+n[:,None,None]**2+n[None,:,None]**2+n[None,None,:]**2)/8)
        exact = (2*cutoff+1)**3*(1+cutoff*(cutoff+1))/8
        direct_variance.append({"cutoff":cutoff,"direct_variance":direct,"closed_variance":exact})
        checks.append({"name":f"gaussian_variance_mode_sum_K{cutoff}","passed":direct==exact})
    data_path = ROOT/'inputs/inhom_casimir_data.npz'
    np.savez(data_path, z=zgrid, rho_dimless=grid[:,0], p_dimless=grid[:,1],
             nec_dimless=grid[:,2], rho_tail_dimless=grid[:,3], p_tail_dimless=grid[:,4],
             varcutoff=cutoff_grid, centered_potential_variance_m1_L2pi=variance)
    report = {
        "model":"real free massive scalar, static R x T_L^3, local Minkowski vacuum subtraction",
        "units":"hbar=c=1; dimensionless outputs L^4 rho and L^4 p",
        "continuum_scope":"Infinite absolutely convergent image sums for renormalized expectation; no configuration-local Q-density continuum construction",
        "renormalization":"Minkowski vacuum expectation set to zero; residual torus Casimir energy retained",
        "points":points,"gaussian_variance_checks":direct_variance,
        "checks":checks,"passed":all(c['passed'] for c in checks),
        "data_sha256":hashlib.sha256(data_path.read_bytes()).hexdigest(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "heat_kernel_quadrature_is_not_interval_certified":True,
        "lattice_tail_bounds_are_analytic":True,
    }
    report_path = ROOT/'work/inhom_casimir_verification.json'
    report_path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({"passed":report['passed'],"checks":len(checks),"points":[{'mL':r['mL'],'rho':r['bessel']['rho_L4'],'p':r['bessel']['pressure_L4'],'max_difference':max(r['rho_method_difference_L4'],r['pressure_method_difference_L4'])}for r in points]},indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
