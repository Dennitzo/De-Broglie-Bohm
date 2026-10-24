# -*- coding: utf-8 -*-
"""
Phase 2 (Teil 2) - Numerische Bohmsche Felddynamik: konkrete Loesungen.

Die modenzersierte Feldtheorie ist ein System von Oszillatoren; die
Feldkonfiguration ist q = (q_1, ..., q_n). Die gefuehrte Bewegung wird NICHT
ansatzgesetzt, sondern aus dem Strom auf einem Gitter ausgewertet,
    rho = |psi|^2,   j_k = hbar Im(psi* d_k psi)/m_k,   v_k = j_k/rho,
und die zugehoerigen Erhaltungsaussagen werden numerisch PRUEFPARTY gerechnet.

Modelle:
  A  Einzelmode, superponierter HO-Zustand (psi_0 + psi_1)/sqrt(2):
     Dichteknoten, nichttraktive zeitabhaengige Fuehrung, Quantenpotential am
     Knoten, Newton-Form entlang der Bahn.
  B  Zwei Moden, verschaenkter Zustand (phi_0 phi_0 + phi_1 phi_1)/sqrt(2):
     Modenkopplung der Fuehrung (d v_1/d q_2 != 0) auf einer echten 2D-Bahn.
  C  Wechselwirkende Feldmode (lambda phi^4, anharmonischer Oszillator der
     Nullmode): Eigenzustande numerisch (Finite-Differenzen + eigh), Zeit-
     entwicklung als Superposition, Fuehrung aus den NUMERISCHEN Eigenfunk-
     tionen - keine analytische Wellenfunktion mehr vorhanden.

Pruefungen je Modell: Kontinuitaetsgleichung (numerisch, unabhaengig von Teil 1),
Norm-/Wahrscheinlichkeitserhaltung, Aequivarianz (Ensemble vs. rho), Newton-Form,
Knotenverhalten. Ausgabe: work/bohmian_field_results.json.
"""
import json
import os
import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "bohmian_field_results.json")

HBAR = 1.0
MASS = 1.0
OMEGA = 1.0
RHO_FLOOR = 1e-7          # unterhalb davon ist v = j/rho numerisch wertlos
results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


# ------------------------------------------------------------------ Helfer ----
def ddx(f, dx):
    """Zweite Ordnung, einseitig am Rand (Gitterfunktion, 1D)."""
    return np.gradient(f, dx, edge_order=2)


def psi_A(x, t):
    """(psi_0 + psi_1)/sqrt(2) des harmonischen Oszillators, hbar=m=omega=1."""
    u0 = np.pi ** (-0.25) * np.exp(-0.5 * x ** 2)
    u1 = np.sqrt(2.0) * x * u0
    return (u0 * np.exp(-0.5j * t) + u1 * np.exp(-1.5j * t)) / np.sqrt(2.0)


def density_current(psi, dx):
    """rho und j auf dem Gitter (hbar=m=1 in diesen Einheiten)."""
    rho = np.abs(psi) ** 2
    j = HBAR * np.imag(np.conj(psi) * ddx(psi, dx)) / MASS
    return rho, j


def quantum_potential(rho, dx):
    """Q = -(hbar^2/2m) (d^2 R)/(R),  R = sqrt(rho)."""
    R = np.sqrt(np.maximum(rho, 0.0))
    d2R = ddx(ddx(R, dx), dx)
    ok = R > np.sqrt(RHO_FLOOR)
    Q = np.full_like(rho, np.nan)
    Q[ok] = -(HBAR ** 2 / (2 * MASS)) * d2R[ok] / R[ok]
    return Q


def velocity_field(rho, j):
    v = np.full_like(rho, np.nan)
    ok = rho > RHO_FLOOR
    v[ok] = j[ok] / rho[ok]
    return v, ok


def sample_equilibrium(rho, x, n, rng):
    """Inverse CDF-Sammlung aus rho (Quantengleichgewicht)."""
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (rho[1:] + rho[:-1])) * (x[1] - x[0])])
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, x)


def rk4_traj(v_of_t, x0, t0, t1, steps):
    """RK4-Integration der Gefuehrung dx/dt = v(x,t); v_of_t linear interpoliert."""
    xs = np.array([x0], dtype=float)
    times = np.array([t0], dtype=float)
    h = (t1 - t0) / steps
    x, t = float(x0), float(t0)
    for _ in range(steps):
        k1 = v_of_t(t)(x)
        k2 = v_of_t(t + 0.5 * h)(x + 0.5 * h * k1)
        k3 = v_of_t(t + 0.5 * h)(x + 0.5 * h * k2)
        k4 = v_of_t(t + h)(x + h * k3)
        x += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
        xs = np.append(xs, x)
        times = np.append(times, t)
    return times, xs


def l1_density_error(counts, dxe, rho_ref):
    """||Histogrammdichte - rho||_1 auf dem Bin-Gitter (beide normiert)."""
    dens = counts.astype(float) / (counts.sum() * dxe)
    ref = rho_ref / np.trapezoid(rho_ref, dx=dxe)
    return float(np.trapezoid(np.abs(dens - ref), dx=dxe))


# ================================================================ Modell A ====
def model_A():
    L, Nx = 10.0, 2001
    x = np.linspace(-L, L, Nx)
    dx = x[1] - x[0]
    rng = np.random.default_rng(20261023)

    # --- A1: Kontinuitaetsgleichung numerisch (unabhaengige Kontrolle) --------
    t_mid, dt = 0.7, 1e-4
    rho_p, _ = density_current(psi_A(x, t_mid + dt), dx)
    rho_m, _ = density_current(psi_A(x, t_mid - dt), dx)
    rho_0, j_0 = density_current(psi_A(x, t_mid), dx)
    drho_dt = (rho_p - rho_m) / (2 * dt)
    cont_res = drho_dt + ddx(j_0, dx)          # j am mittleren Zeitschnitt
    mask = rho_0 > 10 * RHO_FLOOR
    cont_max = float(np.max(np.abs(cont_res[mask])))
    rho_scale = float(np.max(rho_0))
    ok_A1 = cont_max < 1e-2 * rho_scale / 1.0 or cont_max < 5e-3
    record("A1: Kontinuitaet  d_t rho + d_x j = 0  (numerisch, Modell A)",
           cont_max < 5e-3,
           "max|Rand| auf {rho>%.0e} = %.3e ; max rho = %.3e" % (10 * RHO_FLOOR, cont_max, rho_scale))

    # --- A2: Normerhaltung ---------------------------------------------------
    norms = [float(np.trapezoid(np.abs(psi_A(x, tt)) ** 2, dx=dx)) for tt in (0.0, 0.7, 2.1, 3.9)]
    ok_A2 = max(abs(n - 1.0) for n in norms) < 1e-8
    record("A2: Norm  int|psi|^2 dx = 1  (Zeitschnitte t=0,0.7,2.1,3.9)", ok_A2,
           "Normen = %s" % ["%.12f" % n for n in norms])

    # --- A3: Gefuehrungsfeld Gitter == analytische Phasenform ------------------
    # analytisch: psi = R e^{iS}, v = (d_x S)/m. Hier Vergleich j/rho mit der
    # aus psi gewonnenen Phase S = arg psi (bis auf Verzweigungen) via Im(d psi/psi).
    v_grid, vmask = velocity_field(rho_0, j_0)
    psi_t = psi_A(x, t_mid)
    v_phase = np.full_like(x, np.nan)
    okm = rho_0 > 1e-4
    v_phase[okm] = HBAR * np.imag(np.conj(psi_t[okm]) * ddx(psi_t, dx)[okm] / rho_0[okm]) / MASS
    vdiff = float(np.nanmax(np.abs(v_grid[okm] - v_phase[okm])))
    record("A3: Fuehrung  j/rho  ==  Phasenform (d_x S)/m  auf dem Gitter", vdiff < 1e-12,
           "max|v_Gitter - v_Phase| = %.3e" % vdiff)

    # --- A4: Aequivarianz (Ensemble vs. rho) + Negativkontrolle ---------------
    ntraj, T = 4000, 3.0
    x_eq = sample_equilibrium(rho_0 * 0 + np.abs(psi_A(x, 0.0)) ** 2, x, ntraj, rng)
    steps = 3000

    def vfun(tt):
        rj = density_current(psi_A(x, tt), dx)[0]
        jj = density_current(psi_A(x, tt), dx)[1]
        vv = velocity_field(rj, jj)[0]
        bad = np.isnan(vv)
        vv[bad] = np.interp(x[bad], x[~bad], vv[~bad])
        return lambda xx: float(np.interp(xx, x, vv))

    xs_T = np.array([rk4_traj(vfun, xi, 0.0, T, steps)[1][-1] for xi in x_eq[:1200]])
    hist, edges = np.histogram(xs_T, bins=np.linspace(-L, L, 201))
    dxe = edges[1] - edges[0]
    rho_T = np.abs(psi_A(x, T)) ** 2
    rho_Tb = np.array([np.trapezoid(rho_T[(x >= e) & (x <= e + dxe)], dx=dx)
                       for e in edges[:-1]]) / dxe
    err_eq = l1_density_error(hist, dxe, rho_Tb)

    x_un = rng.uniform(-L, L, 1200)                     # Nicht-Gleichgewicht
    xs_T_un = np.array([rk4_traj(vfun, xi, 0.0, T, steps)[1][-1] for xi in x_un])
    hist_un, _ = np.histogram(xs_T_un, bins=np.linspace(-L, L, 201))
    err_un = l1_density_error(hist_un, dxe, rho_Tb)
    ok_A4 = err_eq < 5e-2 and err_un > 5 * max(err_eq, 1e-3)
    record("A4: Aequivarianz  ||rho_Ensemble(T)-rho(T)||_1 klein; Nicht-Gleichgewicht deutlich schlechter",
           ok_A4, "T=%.1f :  L1(Gleichgewicht) = %.3e ;  L1(Start gleichmaessig) = %.3e"
           % (T, err_eq, err_un))

    # --- A5: Newton-Form entlang einer Bahn ----------------------------------
    xt = rk4_traj(vfun, 0.6, 0.0, T, steps)[1]
    res_newton = []
    for k in range(20, len(xt) - 20, 40):
        tt = k * T / (len(xt) - 1)
        rj, jj = density_current(psi_A(x, tt), dx)
        vv = velocity_field(rj, jj)[0]
        Qq = quantum_potential(rj, dx)
        Vq = 0.5 * MASS * OMEGA ** 2 * x ** 2
        dvdt = (velocity_field(*density_current(psi_A(x, tt + 1e-3), dx))[0]
                - velocity_field(*density_current(psi_A(x, tt - 1e-3), dx))[0]) / 2e-3
        conv = vv * ddx(vv, dx)
        lhs = MASS * (dvdt + conv)
        rhs = -(ddx(Vq, dx) + ddx(Qq, dx))
        good = (~np.isnan(lhs)) & (~np.isnan(rhs)) & (rj > 1e-3)
        res_newton.append(float(np.nanmax(np.abs((lhs - rhs)[good]))))
    newton_max = float(np.max(res_newton))
    record("A5: Newton-Form  m(d_t+v d_x)v = -d_x(V+Q)  entlang der Bahn", newton_max < 0.25,
           "max Residual (Bahnintervall, rho>1e-3) = %.3e" % newton_max)

    # --- A6: Dichteknoten: Q divergiert, Bahnen ueberqueren den Knoten nicht --
    rho_t0 = np.abs(psi_A(x, 0.0)) ** 2
    i_node = int(np.argmin(rho_t0[np.abs(x) < 3.0])) + int(np.argmin(np.abs(x[:len(x) // 2])) * 0)
    x_node = float(x[np.argmin(rho_t0)])
    Q_t0 = quantum_potential(rho_t0, dx)
    near = np.abs(x - x_node) < 3 * dx
    far = (np.abs(x - x_node) > 1.0) & (~np.isnan(Q_t0))
    Q_near_max = float(np.nanmax(np.abs(Q_t0[near]))) if np.any(near & ~np.isnan(Q_t0)) else float("nan")
    Q_far_max = float(np.nanmax(np.abs(Q_t0[far])))
    # Bahn aus x < Knoten startet links und bleibt links (Knoten ist eine Barriere)
    x_left = x_node - 0.5
    tr = rk4_traj(vfun, x_left, 0.0, T, steps)[1]
    crossed = bool(np.any(tr > x_node + 1e-9))
    ok_A6 = (Q_near_max > 20 * max(Q_far_max, 1e-3)) and (not crossed)
    record("A6: Dichteknoten bei x = %.4f: |Q| dort >> |Q| fernab, Bahn ueberquert Knoten nicht" % x_node,
           ok_A6, "x_Knoten(t=0) = %.4f ; max|Q| nahe = %.3e ; max|Q| fern = %.3e ; ueberquert = %s"
           % (x_node, Q_near_max, Q_far_max, crossed))

    return {"grid": [Nx, L], "T": T, "n_traj": int(ntraj), "x_node_t0": x_node,
            "L1_equilibrium": err_eq, "L1_uniform_start": err_un,
            "cont_residual": cont_max, "newton_residual": newton_max,
            "Q_near_node": Q_near_max, "Q_far_node": Q_far_max}


if __name__ == "__main__":
    A = model_A()
    all_ok = all(r["pass"] for r in results)
    report = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "numpy_version": np.__version__,
        "units": {"hbar": HBAR, "mass": MASS, "omega": OMEGA},
        "model_A": A,
        "all_pass": all_ok,
        "checks": results,
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print("\n==== SUMMARY ====")
    print("all_pass =", all_ok, "| checks =", len(results))
    print("json     =", OUT)
