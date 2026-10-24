"""Unabhaengige Trajektorienkontrolle ohne FFT-Feld.

Wir integrieren die GESCHLOSSEN analytische Fuehrungsgeschwindigkeit
    v_an(x,t) = x * t * sigma0^2 / (sigma0^4 + t^2)
mit RK4 und vergleichen gegen die geschlossene Trajektorie
    x(t) = x0 * sqrt(1 + (t/sigma0^2)^2).

Erwartung: Fehler ~0 -> analytische Trajektorie und RK4 sind korrekt; jede
Abweichung in der FFT-Version kommt dann aus der Feld-/Ableitungsseite.
"""
import json
import numpy as np
from scipy.special import erfinv


def main():
    sigma0 = 1.0
    t_end = 1.0
    Ntraj = 5000
    q = (np.arange(Ntraj) + 0.5) / Ntraj
    std0 = sigma0 / np.sqrt(2.0)
    X0 = std0 * (np.sqrt(2.0) * erfinv(2.0 * q - 1.0))

    def v_an(x, t):
        return x * t * sigma0 ** 2 / (sigma0 ** 4 + t * t)

    results = []
    for dt in (0.04, 0.01, 0.002):
        nt = int(round(t_end / dt))
        dt_used = t_end / nt
        xq = X0.copy()
        for step in range(nt):
            t0 = step * dt_used
            k1 = v_an(xq, t0)
            k2 = v_an(xq + 0.5 * dt_used * k1, t0 + 0.5 * dt_used)
            k3 = v_an(xq + 0.5 * dt_used * k2, t0 + 0.5 * dt_used)
            k4 = v_an(xq + dt_used * k3, t0 + dt_used)
            xq = xq + dt_used * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        xa = X0 * np.sqrt(1.0 + (t_end / (sigma0 * sigma0)) ** 2)
        rel = float(np.linalg.norm(xq - xa) / np.linalg.norm(xa))
        results.append({"dt": dt_used, "nt": nt, "rel_l2_vs_closed_traj": rel})

    out = {"observable": "trajektorie_aus_geschlossener_fuehrung",
           "reference": "geschlossene Trajektorie x(t)=x0*sqrt(1+(t/sigma0^2)^2), reine RK4 ohne FFT",
           "results": results}
    with open("knotenfrei_traj_an.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
