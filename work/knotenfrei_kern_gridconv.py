"""Direkte Loesungsdistanz-Konvergenz fuer den knotenfreien Kernkontrollfall.

Konvergenz wird hier ECHT als Loesungsdistanz gemessen: Dieselben determini-
istischen Startpunkte X0 werden mit dt_fine und dt_coarse gefuehrt, und die
Abweichung der Trajektorien-Endpunkte wird direkt gegeneinander gerechnet
(relativer Abstand der Loesungen auf gemeinsamer Vergleichsbasis), NICHT als
Differenz zweier skalarer Fehlernormen.

Zusaetzlich: Abweichung der feinen Loesung von der geschlossenen analytischen
Trajektorie (accuracy der Konvergenz).
"""
import argparse
import json
import numpy as np
from scipy.special import erfinv


def psi_initial(x, sigma0):
    return (np.pi * sigma0 ** 2) ** (-0.25) * np.exp(-x * x / (2.0 * sigma0 ** 2))


def traj_analytic(x0, t, sigma0):
    return x0 * np.sqrt(1.0 + (t / (sigma0 * sigma0)) ** 2)


def interp1d_uniform(grid_x, f, xq):
    dx = grid_x[1] - grid_x[0]
    pos = (xq - grid_x[0]) / dx
    i = np.clip(np.floor(pos).astype(int), 0, len(grid_x) - 2)
    w = pos - i
    return f[i] * (1.0 - w) + f[i + 1] * w


def integrate(N, Lx, sigma0, t_end, dt, X0):
    m = 1.0
    dx = Lx / N
    x = -Lx / 2.0 + dx * np.arange(N)
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    nt = int(round(t_end / dt))
    dt = t_end / nt
    psi_hat_t = np.fft.fft(psi_initial(x, sigma0))
    factor = np.exp(-1.0j * (k * k) * dt / (2.0 * m))
    xq = X0.copy()
    for _ in range(nt):
        psi_t = np.fft.ifft(psi_hat_t)
        dpsi_t = np.fft.ifft(1.0j * k * psi_t)
        rho_t = np.abs(psi_t) ** 2
        v_t = np.imag(np.conj(psi_t) * dpsi_t) / np.maximum(rho_t, 1e-300)

        def f(p):
            return interp1d_uniform(x, v_t, p)

        k1 = f(xq)
        k2 = f(xq + 0.5 * dt * k1)
        k3 = f(xq + 0.5 * dt * k2)
        k4 = f(xq + dt * k3)
        xq = xq + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        psi_hat_t = psi_hat_t * factor
    return xq, nt, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=4096)
    ap.add_argument("--Lx", type=float, default=40.0)
    ap.add_argument("--sigma0", type=float, default=1.0)
    ap.add_argument("--t_end", type=float, default=1.0)
    ap.add_argument("--Ntraj", type=int, default=20000)
    ap.add_argument("--dts", type=str, default="0.04,0.02,0.01")
    ap.add_argument("--out", type=str, default="knotenfrei_gridconv.json")
    args = ap.parse_args()

    q = (np.arange(args.Ntraj) + 0.5) / args.Ntraj
    std0 = args.sigma0 / np.sqrt(2.0)
    X0 = std0 * (np.sqrt(2.0) * erfinv(2.0 * q - 1.0))

    dts = [float(s) for s in args.dts.split(",") if s.strip() != ""]
    solutions = {}
    info = {}
    for dt in dts:
        xq, nt, dt_used = integrate(args.N, args.Lx, args.sigma0, args.t_end, dt, X0)
        solutions[dt] = xq
        info[dt] = {"dt_used": dt_used, "nt": nt}

    # Loesungen auf gemeinsamer Basis (identische Startpunkte X0) vergleichen.
    x_an = traj_analytic(X0, args.t_end, args.sigma0)
    ref_norm = float(np.linalg.norm(x_an))

    metrics = []
    for dt in dts:
        e = float(np.linalg.norm(solutions[dt] - x_an) / ref_norm)
        metrics.append({"dt": dt, "rel_l2_vs_analytic": e, "nt": info[dt]["nt"]})

    # Direkte Loesungsdistanz zwischen benachbarten Verfeinerungsstufen.
    grid_changes = []
    for i in range(len(dts) - 1):
        coarse, fine = dts[i], dts[i + 1]
        dist = float(np.linalg.norm(solutions[fine] - solutions[coarse]) / np.linalg.norm(solutions[fine]))
        grid_changes.append({"coarse_dt": coarse, "fine_dt": fine, "solution_rel_l2_distance": dist})

    # Konvergenzordnung aus Loesungsdistanzen (logarithmisch).
    order = None
    if len(grid_changes) >= 2:
        d1 = grid_changes[0]["solution_rel_l2_distance"]
        d2 = grid_changes[1]["solution_rel_l2_distance"]
        if d1 > 0 and d2 > 0:
            order = float(np.log(d1 / d2) / np.log(2.0))

    finest = dts[-1]
    out = {
        "observable": "bohm_trajektorien_loesungsdistanz_dt_verfeinerung",
        "reference": "geschlossene analytische Trajektorie x(t)=x0*sqrt(1+(t/sigma0^2)^2); Loesungen mit identischen Startpunkten direkt verglichen",
        "config": {"N": args.N, "Lx": args.Lx, "dx": args.Lx / args.N, "sigma0": args.sigma0,
                   "t_end": args.t_end, "dts": dts, "Ntraj": args.Ntraj},
        "per_dt": metrics,
        "grid_change_direct": grid_changes,
        "convergence_order": order,
        "finest_dt": finest,
        "finest_rel_l2_vs_analytic": float(np.linalg.norm(solutions[finest] - x_an) / ref_norm),
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
