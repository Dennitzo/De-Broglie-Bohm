"""Knotenfreier 1D-Kern-Kontrollfall: freies sich ausbreitendes Gauss-Paket.

Ziel: die komplette numerische Bohm-Pipeline (Fourier-Feld, spektrale Ableitung,
Interpolation, RK4-Trajektorien, Histogrammdichte) gegen drei UNABHAENGIGE
analytische Referenzen zu pruefen:
  (1) Feld      psi_num(x,t) vs geschlossene analytische Loesung
  (2) Fuehrung  v_num(x,t)  vs analytisch v = x*t*sigma0^2/(sigma0^4+t^2)
  (3) Dichte    Ensemble-Histogramm bei t_end vs exakte GauSSdichte s(t)

Das Paket hat p0=0 und keinen Rand in Kernnaehe -> keine Knoten, keine
Regularisierung noetig. Alle Einheiten: hbar=m=1.

Es wird genau eine Konfiguration pro Aufruf ausgefuehrt (wiederaufnehmbar).
Keine grossen .npz; nur schlanke JSON und optional ein kleiner Plot.
"""
import argparse
import json
import numpy as np


def gauss_density(x, sigma):
    return np.exp(-x * x / (sigma * sigma)) / (np.sqrt(np.pi) * sigma)


def psi_initial(x, sigma0):
    return (np.pi * sigma0 ** 2) ** (-0.25) * np.exp(-x * x / (2.0 * sigma0 ** 2))


def psi_analytic(x, t, sigma0):
    tau = t / (sigma0 * sigma0)
    denom = 1.0 + 1.0j * tau
    amp = (np.pi * sigma0 ** 2) ** (-0.25) / np.sqrt(denom)
    return amp * np.exp(-x * x / (2.0 * sigma0 ** 2 * denom))


def v_analytic(x, t, sigma0):
    # Bohm-Geschwindigkeit des freien propagierenden Gauss-Pakets (p0=0)
    return x * t * sigma0 ** 2 / (sigma0 ** 4 + t * t)


def interp1d_uniform(grid_x, f, xq):
    # lineare Interpolation auf gleichmaessigem Gitter
    dx = grid_x[1] - grid_x[0]
    pos = (xq - grid_x[0]) / dx
    i = np.clip(np.floor(pos).astype(int), 0, len(grid_x) - 2)
    w = pos - i
    return f[i] * (1.0 - w) + f[i + 1] * w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=4096)
    ap.add_argument("--Lx", type=float, default=40.0)
    ap.add_argument("--sigma0", type=float, default=1.0)
    ap.add_argument("--t_end", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--Np", type=int, default=50000)
    ap.add_argument("--seeds", type=str, default="1,2,3")
    ap.add_argument("--xlo", type=float, default=-8.0)
    ap.add_argument("--xhi", type=float, default=8.0)
    ap.add_argument("--NBins", type=int, default=60)
    ap.add_argument("--out", type=str, default="knotenfrei_pilot.json")
    ap.add_argument("--plot", type=str, default="")
    args = ap.parse_args()

    m = 1.0
    sigma0 = args.sigma0
    N = args.N
    Lx = args.Lx
    dx = Lx / N
    x = -Lx / 2.0 + dx * np.arange(N)
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)

    # exakte Fourier-Zeitentwicklung des freien Feldes
    def propagate(psi_hat, dt, nsteps):
        factor = np.exp(-1.0j * (k * k) * dt / (2.0 * m))
        for _ in range(nsteps):
            psi_hat = psi_hat * factor
        return psi_hat

    nt = int(round(args.t_end / args.dt))
    dt = args.t_end / nt  # exakte Aufloesung des Zeitschritts

    # Anfangsfeld
    psi0 = psi_initial(x, sigma0)
    psi_hat0 = np.fft.fft(psi0)

    # Feld an t_end: numerisch vs analytisch
    psi_hat_end = propagate(psi_hat0.copy(), dt, nt)
    psi_num = np.fft.ifft(psi_hat_end)
    psi_an = psi_analytic(x, args.t_end, sigma0)
    field_err = np.linalg.norm(psi_num - psi_an) / np.linalg.norm(psi_an)

    # Randamplitude fuer Randfreiheit
    edge_mass = (np.sum(np.abs(psi_num[np.abs(x) > 18.0]) ** 2)) * dx
    rho_num = np.abs(psi_num) ** 2
    edge_amp = float(np.sqrt(np.max(rho_num[np.abs(x) > 18.0]))) if np.any(np.abs(x) > 18.0) else 0.0

    # Fuehrungsgeschwindigkeit aus numerischem Feld vs analytisch (Kernbereich)
    dpsi_num = np.fft.ifft(1.0j * k * psi_hat_end)
    v_num = np.imag(np.conj(psi_num) * dpsi_num) / np.maximum(rho_num, 1e-300)
    core = np.abs(x) <= 3.0
    v_err = np.max(np.abs(v_num[core] - v_analytic(x[core], args.t_end, sigma0)))
    v_rel = np.linalg.norm((v_num - v_analytic(x, args.t_end, sigma0))[core]) / np.linalg.norm(v_analytic(x[core], args.t_end, sigma0))

    # Trajektorien-Ensemble
    edges = np.linspace(args.xlo, args.xhi, args.NBins + 1)
    binw = edges[1] - edges[0]
    # analytische Bin-Dichte (unabhaengig, auf Gitter integriert)
    rho_an_end = np.abs(psi_an) ** 2
    a_hist, _ = np.histogram(x, bins=edges, weights=rho_an_end)
    a_density = a_hist * dx / binw
    a_norm = np.linalg.norm(a_density)

    seeds = [int(s) for s in args.seeds.split(",") if s.strip() != ""]
    density_metrics = []
    first_hist = None
    first_X = None
    for seed in seeds:
        rng = np.random.default_rng(seed)
        X = rng.normal(0.0, sigma0 / np.sqrt(2.0), size=args.Np)
        # RK4 mit numerischem Geschwindigkeitsfeld bei festem t_end (Feld zeitunabh. fuer Momentaufnahme
        # -> hier korrekte zeitabhaengige Fuehrung: wir propagieren Feld mit und integrieren simultan)
        psi_hat_t = psi_hat0.copy()
        xq = X.copy()
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
            psi_hat_t = psi_hat_t * np.exp(-1.0j * (k * k) * dt / (2.0 * m))

        h, _ = np.histogram(xq, bins=edges)
        h_density = h / (args.Np * binw)
        rel_l2 = float(np.linalg.norm(h_density - a_density) / a_norm)
        density_metrics.append(rel_l2)
        if first_hist is None:
            first_hist = h_density
            first_X = xq

    mean_d = float(np.mean(density_metrics))
    std_d = float(np.std(density_metrics, ddof=1)) if len(density_metrics) > 1 else None

    out = {
        "observable": "ensemble_density_vs_analytic_gauss_and_fuehrung",
        "reference": "geschlossene analytische Loesung des freien Gauss-Pakets (Feld, v, Dichte), p0=0",
        "config": {
            "N": N, "Lx": Lx, "sigma0": sigma0, "t_end": args.t_end, "dt": dt,
            "nt": nt, "Np": args.Np, "seeds": seeds, "NBins": args.NBins,
            "xlo": args.xlo, "xhi": args.xhi,
        },
        "field_rel_l2": float(field_err),
        "v_max_abs_err_core": float(v_err),
        "v_rel_l2_core": float(v_rel),
        "edge_mass": float(edge_mass),
        "edge_amplitude": float(edge_amp),
        "density_rel_l2_mean": mean_d,
        "density_rel_l2_std": std_d,
        "density_rel_l2_per_seed": density_metrics,
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            centers = 0.5 * (edges[:-1] + edges[1:])
            fig, ax = plt.subplots(1, 3, figsize=(15, 4))
            ax[0].plot(x, np.abs(psi_an) ** 2, "r-", label="analytic")
            ax[0].plot(x, np.abs(psi_num) ** 2, "b--", label="numeric")
            ax[0].set_title("Feld |psi(x,t_end)|^2")
            ax[0].legend()
            ax[1].plot(x[core], v_analytic(x[core], args.t_end, sigma0), "r-", label="analytic v")
            ax[1].plot(x[core], v_num[core], "b--", label="numeric v")
            ax[1].set_title("Fuehrung v(x,t_end)")
            ax[1].legend()
            ax[2].bar(centers, first_hist, width=binw, alpha=0.6, label="Bohm-Ensemble")
            ax[2].plot(x, np.abs(psi_an) ** 2, "r-", label="analytic")
            ax[2].set_title("Enddichte t=%.3f relL2=%.4f" % (args.t_end, mean_d))
            ax[2].legend()
            fig.tight_layout()
            fig.savefig(args.plot, dpi=120)
            print("plot_saved", args.plot)
        except Exception as e:
            print("plot_skipped", repr(e))


if __name__ == "__main__":
    main()
