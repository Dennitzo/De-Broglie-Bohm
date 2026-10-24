"""Loesungsdistanz-Konvergenz mit dem VERIFIZIERTEN Akzept-Fuehrungsalgorithmus.

Der Trajektorienloop ist woertlich aus knotenfrei_kern_accept.py uebernommen
(dort verifiziert: dt=0.002 -> rel_Fehler 5.0e-4, Ordnung 1). Dieselben
deterministischen Startpunkte X0 werden ueber mehrere dt gefuehrt und die
Loesungen DIREKT gegeneinander verglichen (Loesungsdistanz), nicht ueber
differenzen skalarer Fehlernormen.
"""
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


def run_traj(N, Lx, sigma0, t_end, dt, X0):
    m = 1.0
    dx = Lx / N
    x = -Lx / 2.0 + dx * np.arange(N)
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    nt = int(round(t_end / dt))
    dt = t_end / nt
    psi_hat0 = np.fft.fft(psi_initial(x, sigma0))
    factor = np.exp(-1.0j * (k * k) * dt / (2.0 * m))
    psi_hat_t = psi_hat0.copy()
    xq = X0.copy()
    for step in range(nt):
        psi_t = np.fft.ifft(psi_hat_t)
        dpsi_t = np.fft.ifft(1.0j * k * psi_hat_t)
        rho_t = np.abs(psi_t) ** 2
        v_t = np.imag(np.conj(psi_t) * dpsi_t) / np.maximum(rho_t, 1e-300)

        def f(vpos):
            return interp1d_uniform(x, v_t, vpos)

        k1 = f(xq)
        k2 = f(xq + 0.5 * dt * k1)
        k3 = f(xq + 0.5 * dt * k2)
        k4 = f(xq + dt * k3)
        xq = xq + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        psi_hat_t = psi_hat_t * factor
    return xq, dt, nt


def main():
    sigma0 = 1.0
    N = 4096
    Lx = 40.0
    t_end = 1.0
    Ntraj = 20000
    q = (np.arange(Ntraj) + 0.5) / Ntraj
    std0 = sigma0 / np.sqrt(2.0)
    X0 = std0 * (np.sqrt(2.0) * erfinv(2.0 * q - 1.0))

    dts = [0.04, 0.02, 0.01, 0.005]
    solutions = {}
    info = {}
    for dt in dts:
        xq, dt_used, nt = run_traj(N, Lx, sigma0, t_end, dt, X0)
        solutions[dt] = xq
        info[dt] = {"dt_used": dt_used, "nt": nt}

    x_an = traj_analytic(X0, t_end, sigma0)
    ref_norm = float(np.linalg.norm(x_an))

    per_dt = []
    for dt in dts:
        e = float(np.linalg.norm(solutions[dt] - x_an) / ref_norm)
        per_dt.append({"dt": info[dt]["dt_used"], "nt": info[dt]["nt"], "rel_l2_vs_analytic": e})

    # Direkte Loesungsdistanz auf IDENTISCHER Vergleichsbasis (X0 fix).
    grid_changes = []
    for i in range(len(dts) - 1):
        coarse, fine = dts[i], dts[i + 1]
        dist = float(np.linalg.norm(solutions[fine] - solutions[coarse]) / np.linalg.norm(solutions[fine]))
        grid_changes.append({"coarse_dt": dts[i], "fine_dt": dts[i + 1], "solution_rel_l2_distance": dist})

    order = None
    if len(grid_changes) >= 2:
        d1 = grid_changes[0]["solution_rel_l2_distance"]
        d2 = grid_changes[1]["solution_rel_l2_distance"]
        d3 = grid_changes[2]["solution_rel_l2_distance"]
        if d1 > 0 and d2 > 0 and d3 > 0:
            order = float(np.mean([np.log(d1 / d2) / np.log(2.0), np.log(d2 / d3) / np.log(2.0)]))

    finest = dts[-1]
    out = {
        "observable": "bohm_trajektorien_loesungsdistanz_dt_verfeinerung",
        "reference": "geschlossene analytische Trajektorie x(t)=x0*sqrt(1+(t/sigma0^2)^2); Loesungen mit identischen Startpunkten direkt verglichen",
        "config": {"N": N, "Lx": Lx, "dx": Lx / N, "sigma0": sigma0, "t_end": t_end,
                   "dts": dts, "Ntraj": Ntraj},
        "per_dt": per_dt,
        "grid_change_direct": grid_changes,
        "convergence_order": order,
        "finest_dt": info[finest]["dt_used"],
        "finest_rel_l2_vs_analytic": float(np.linalg.norm(solutions[finest] - x_an) / ref_norm),
    }
    with open("knotenfrei_conv.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
