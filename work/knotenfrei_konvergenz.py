"""Knotenfreier Kernkontrollfall - Abnahme mit echter Loesungsdistanz.

Selbststaendige Referenz: freies Gauss-Paket (p0=0, hbar=m=1) hat eine
geschlossene analytische Loesung fuer Feld, Fuehrung, Trajektorie und Dichte.

  psi(x,0)   = (pi sigma0^2)^(-1/4) exp(-x^2/(2 sigma0^2))
  psi(x,t)   = (pi sigma0^2)^(-1/4) (1+i t/sigma0^2)^(-1/2)
               exp(-x^2/(2(sigma0^2+i t)))
  |psi(x,t)|^2= (pi sigma_t^2)^(-1/2) exp(-x^2/sigma_t^2),
               sigma_t^2=(sigma0^4+t^2)/sigma0^2
  S(x,t)     = x^2 t / (2(sigma0^4+t^2))
  v(x,t)     = dS/dx = x t/(sigma0^4+t^2)
  x(t)       = x0 sqrt(1+(t/sigma0^2)^2)

Pipeline: Feld numerisch (Spektral-Split, fuer V=0 exakt), Fuehrung aus
numerischem Feld (Spektralableitung), Trajektorie per RK4, Enddichte als
Histogramm. Konvergenz ist die Loesungsdistanz der Dichte-Histogramme
zwischen Verfeinerungsstufen mit IDENTISCHEN Startpunkten und IDENTISCHEN
Bins; NICHT die Differenz zweier skalarer Fehlernormen.
"""
import json
import numpy as np


def gauss_dichte(x, sigma0, t):
    sigma_t2 = (sigma0**4 + t**2) / sigma0**2
    return np.exp(-x**2 / sigma_t2) / np.sqrt(np.pi * sigma_t2)


def anfangspunkte(rng, sigma0, n):
    # rho(x,0) ~ exp(-x^2/sigma0^2) => Var = sigma0^2/2
    return rng.normal(0.0, sigma0 / np.sqrt(2.0), size=n)


def feld_loesung(Nx, Lx, sigma0):
    x = (np.arange(Nx) - Nx / 2) * (Lx / Nx)
    psi = (np.pi * sigma0**2) ** (-0.25) * np.exp(-x**2 / (2 * sigma0**2))
    return x, psi


def spektral_d1(psi, dx):
    Nx = psi.size
    k = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)
    return np.fft.fftshift(np.fft.ifft(1j * k * np.fft.fft(np.fft.fftshift(psi))))


def propagiere_freies_feld(psi0, dx, Lx, dt, nt):
    Nx = psi0.size
    k = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)
    phase = np.exp(-0.5j * k**2 * dt)
    psi = psi0.copy()
    speicher = [psi.copy()]
    for _ in range(nt):
        psi = np.fft.ifft(phase * np.fft.fft(psi))
        speicher.append(psi.copy())
    return speicher


def fuehrungsfeld(psi, dx):
    dpsi = spektral_d1(psi, dx)
    rho = (psi.conj() * psi).real
    num = (psi.conj() * dpsi).imag
    return num / np.maximum(rho, 1e-300)


def trajektorie_ende(x0, speicher, dx, nt, dt, t_end):
    # RK4 mit linearer Zeitinterpolation des Geschwindigkeitsfelds.
    x = x0.copy()
    Nx = speicher[0].size
    xg = (np.arange(Nx) - Nx / 2) * dx

    def vel(t, pos):
        i0 = min(int(t / dt), nt - 1)
        w = (t - i0 * dt) / dt
        psi = (1 - w) * speicher[i0] + w * speicher[i0 + 1]
        v = fuehrungsfeld(psi, dx)
        return np.interp(pos, xg, v)

    t = 0.0
    while t < t_end - 1e-12:
        h = min(dt, t_end - t)
        k1 = vel(t, x)
        k2 = vel(t + h / 2, x + h / 2 * k1)
        k3 = vel(t + h / 2, x + h / 2 * k2)
        k4 = vel(t + h, x + h * k3)
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return x


def dichte_histogramm(xend, bins):
    counts, edges = np.histogram(xend, bins=bins)
    rho = counts.astype(float)
    rho /= rho.sum()
    zentren = 0.5 * (edges[:-1] + edges[1:])
    db = edges[1] - edges[0]
    return zentren, rho / db, rho


def loesungsdistanz(a, b):
    num = np.sqrt(np.sum((a - b) ** 2))
    den = np.sqrt(np.sum(b**2))
    return num / den if den > 0 else float("inf")


def dichte_gegen_analytisch(zentren, rho_hist, sigma0, t, db):
    ref = gauss_dichte(zentren, sigma0, t)
    num = np.sqrt(np.sum((rho_hist - ref) ** 2) * db)
    den = np.sqrt(np.sum(ref**2) * db)
    return num / den


def feld_gegen_analytisch(psi_num, x, sigma0, t):
    psi_an = (np.pi * sigma0**2) ** (-0.25) * (1 + 1j * t / sigma0**2) ** (-0.5) \
        * np.exp(-x**2 / (2 * (sigma0**2 + 1j * t)))
    num = np.sqrt(np.sum(np.abs(psi_num - psi_an) ** 2))
    den = np.sqrt(np.sum(np.abs(psi_an) ** 2))
    return num / den


def konfiguration(Nx, Lx, sigma0, t_end, dt, Np, NBins, xlo, xhi, seed):
    dx = Lx / Nx
    rng = np.random.default_rng(seed)
    x0 = anfangspunkte(rng, sigma0, Np)
    xg, psi0 = feld_loesung(Nx, Lx, sigma0)
    nt = int(round(t_end / dt))
    speicher = propagiere_freies_feld(psi0, dx, Lx, dt, nt)
    xend = trajektorie_ende(x0, speicher, dx, nt, dt, t_end)
    bins = np.linspace(xlo, xhi, NBins + 1)
    zentren, rho_hist, rho_mass = dichte_histogramm(xend, bins)
    db = bins[1] - bins[0]
    acc = dichte_gegen_analytisch(zentren, rho_hist, sigma0, t_end, db)
    feld = feld_gegen_analytisch(speicher[nt], xg, sigma0, t_end)
    edge = max(np.max(np.abs(speicher[0][:Nx // 10])),
               np.max(np.abs(speicher[0][-Nx // 10:])))
    return dict(zentren=zentren, rho_hist=rho_hist, acc=acc, feld=feld,
                edge=float(edge), mass=rho_mass)


def main():
    sigma0, t_end, NBins, xlo, xhi, Np = 1.0, 1.0, 60, -8.0, 8.0, 40000
    Lx = 40.0
    seeds = [1, 2, 3, 4, 5]

    dt_stufen = [0.02, 0.01, 0.005]
    dt_laufe = {}
    for dt in dt_stufen:
        dt_laufe[dt] = konfiguration(4096, Lx, sigma0, t_end, dt, Np,
                                     NBins, xlo, xhi, seed=1)
    conv_dt_001_002 = loesungsdistanz(dt_laufe[0.02]["rho_hist"],
                                      dt_laufe[0.01]["rho_hist"])
    conv_dt_0005_001 = loesungsdistanz(dt_laufe[0.01]["rho_hist"],
                                       dt_laufe[0.005]["rho_hist"])

    dx_stufen = [2048, 4096, 8192]
    dx_laufe = {}
    for Nx in dx_stufen:
        dx_laufe[Nx] = konfiguration(Nx, Lx, sigma0, t_end, 0.005, Np,
                                     NBins, xlo, xhi, seed=1)
    conv_dx_2048_4096 = loesungsdistanz(dx_laufe[2048]["rho_hist"],
                                        dx_laufe[4096]["rho_hist"])
    conv_dx_4096_8192 = loesungsdistanz(dx_laufe[4096]["rho_hist"],
                                        dx_laufe[8192]["rho_hist"])

    acc_seed = []
    for s in seeds:
        r = konfiguration(4096, Lx, sigma0, t_end, 0.005, Np, NBins, xlo,
                          xhi, seed=s)
        acc_seed.append(r["acc"])
    acc_seed = np.array(acc_seed)

    acc_fein = dx_laufe[8192]["acc"]
    feld_fein = dx_laufe[8192]["feld"]
    edge_fein = dx_laufe[8192]["edge"]

    out = {
        "observable": "bohm_pipeline_vs_closed_form_free_gauss_packet",
        "reference": "geschlossene analytische Loesung freies Gauss-Paket p0=0: "
                     "Feld, Trajektorie x(t)=x0*sqrt(1+(t/sigma0^2)^2), Born-Dichte",
        "config": {
            "sigma0": sigma0, "t_end": t_end, "NBins": NBins,
            "xlo": xlo, "xhi": xhi, "Np": Np,
            "dt_stufen": dt_stufen, "dx_Nx_stufen": dx_stufen, "seeds": seeds,
        },
        "metrics": {
            "accuracy": float(acc_fein),
            "convergence_dt": float(max(conv_dt_001_002, conv_dt_0005_001)),
            "convergence_dx": float(max(conv_dx_2048_4096, conv_dx_4096_8192)),
            "uncertainty": float(acc_seed.std(ddof=1)),
        },
        "detail": {
            "feld_rel_l2_fein": float(feld_fein),
            "edge_amplitude_fein": float(edge_fein),
            "conv_dt_002_vs_001": float(conv_dt_001_002),
            "conv_dt_001_vs_0005": float(conv_dt_0005_001),
            "conv_dx_2048_vs_4096": float(conv_dx_2048_4096),
            "conv_dx_4096_vs_8192": float(conv_dx_4096_8192),
            "acc_per_seed": [float(v) for v in acc_seed],
            "acc_dt002": float(dt_laufe[0.02]["acc"]),
            "acc_dt001": float(dt_laufe[0.01]["acc"]),
            "acc_dt0005": float(dt_laufe[0.005]["acc"]),
        },
    }
    with open("knotenfrei_konvergenz_metrics.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out["metrics"], indent=2))
    print("detail:", json.dumps(out["detail"], indent=2))


if __name__ == "__main__":
    main()
