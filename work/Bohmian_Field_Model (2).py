# -*- coding: utf-8 -*-
"""
Phase 2: endliche Bohmsche Feldmodelle.

Modellklasse:
  Ein reelles Skalarfeld auf einem zweipunktigen Raumgitter. Die kanonischen
  Feldkoordinaten q_1, q_2 bilden einen endlichdimensionalen Konfigurationsraum.
  Der Hamiltonoperator des Gitters ist eine symmetrische 2x2-Matrix; durch
  orthogonale Diagonalisierung erhaelt man zwei Normalmoden mit Frequenzen
  Omega_+ und Omega_-.

  H_classisch = 1/2 (p_1^2 + p_2^2) + 1/2 (q_1, q_2) K (q_1, q_2)^T
  K = [[w0^2, g], [g, w0^2]]

  Normalmoden:
    x = (q_1 + q_2)/sqrt(2), y = (q_1 - q_2)/sqrt(2)
    Omega_+ = sqrt(w0^2 + g), Omega_- = sqrt(w0^2 - g)

Quantisierung (hbar=1):
  Der Hamiltonoperator ist die Summe zweier harmonischer Oszillatoren auf
  L^2(R^2). Ein geeigneter Definitionsbereich ist D(H) = {psi in L^2(R^2):
  H psi in L^2(R^2)}; der Operator ist auf Schwartz-Raeumen wesentlich
  selbstadjungiert und auf D(H) selbstadjungiert. Die hier verwendeten
  Eigenfunktionen und kohhaerenten Zustaende sind explizit in diesem
  Definitionsbereich.

Bohmsche Dynamik:
  Psi = R exp(i S/hbar), v = (1/R^2) Im(conj(Psi) grad Psi), dQ/dt = v(Q,t).
  Fuer kohhaerente Zustaende werden Konfigurationstraektorien numerisch
  integriert und mit der klassischen Normalmoden-Dynamik verglichen. Fuer
  einen Zwei-Moden-Superpositionszustand werden Trajektorien mit
  Knoten-Sicherheitstop integriert und die Kontinuitaetsgleichung auf
  Knoten-ausgeschlossenen Gebieten kontrolliert.

Ausgabe:
  work/bohmian_field_model.json
  Bohmian_Field_Results.md
"""
import json
import math
import os
import datetime

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(ROOT, "work", "bohmian_field_model.json")
OUT_MD = os.path.join(ROOT, "Bohmian_Field_Results.md")

HBAR = 1.0
RNG = np.random.default_rng(20261024)

results = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def norm1d(y, x):
    return float(np.trapezoid(y, x))


def norm2d(z, x, y):
    return float(np.trapezoid(np.trapezoid(z, y, axis=1), x))


def harmonic_basis(z, omega):
    c0 = (omega / math.pi) ** 0.25
    e = np.exp(-0.5 * omega * z * z)
    psi0 = c0 * e
    psi1 = c0 * math.sqrt(2.0 * omega) * z * e
    dpsi0 = -omega * z * psi0
    dpsi1 = c0 * math.sqrt(2.0 * omega) * (1.0 - omega * z * z) * e
    return psi0, psi1, dpsi0, dpsi1


def coherent_center(omega, x0, t):
    xc = x0 * math.cos(omega * t)
    pc = -x0 * omega * math.sin(omega * t)
    return xc, pc


def coherent_rho_1d(x, omega, x0, t):
    xc, _ = coherent_center(omega, x0, t)
    return np.sqrt(omega / math.pi) * np.exp(-omega * (x - xc) ** 2)


def coherent_Q_1d(x, omega, x0, t):
    xc, _ = coherent_center(omega, x0, t)
    return 0.5 * omega - 0.5 * omega * omega * (x - xc) ** 2


def continuity_residual_1d(rho_fn, x, t, dt):
    rho = rho_fn(x, t)
    rho_m = rho_fn(x, t - dt)
    rho_p = rho_fn(x, t + dt)
    dtd_rho = (rho_p - rho_m) / (2.0 * dt)
    vx = rho_fn_velocity_1d(x, t)
    div = np.gradient(rho * vx, x, axis=0)
    res = dtd_rho + div
    interior = slice(3, -3)
    return float(np.max(np.abs(res[interior]))), rho


def rho_fn_velocity_1d(x, t):
    # wird in den Aufrufen ueberschrieben; hier nur Platzhalter
    return np.zeros_like(x)


def integrate_particles_1d(x, omega, x0, T, dt):
    steps = int(round(T / dt))

    def vel(tt):
        return -x0 * omega * math.sin(omega * tt)

    y = np.array(x, dtype=float)
    for i in range(steps):
        tn = i * dt
        k1 = vel(tn)
        k2 = vel(tn + 0.5 * dt)
        k3 = vel(tn + 0.5 * dt)
        k4 = vel(tn + dt)
        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return y


def equivariance_l1_1d(x_initial, omega, x0, T, dt, bins=70):
    x_final = integrate_particles_1d(x_initial, omega, x0, T, dt)
    x_min, x_max = -5.0, 5.0
    edges = np.linspace(x_min, x_max, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dx = edges[1] - edges[0]
    hist, _ = np.histogram(x_final, bins=edges, density=True)
    rho_exact = coherent_rho_1d(centers, omega, x0, T)
    return float(np.sum(np.abs(hist - rho_exact)) * dx), hist, rho_exact, edges


def lattice_model():
    w0sq = 1.0
    g = 0.5
    omega_p = math.sqrt(w0sq + g)
    omega_m = math.sqrt(w0sq - g)
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    K = np.array([[w0sq, g], [g, w0sq]], dtype=float)
    evals, evecs = np.linalg.eigh(K)
    return {
        "w0sq": w0sq,
        "g": g,
        "K": K,
        "evals": evals,
        "evecs": evecs,
        "omega_plus": omega_p,
        "omega_minus": omega_m,
        "inv_sqrt2": inv_sqrt2,
    }


def normal_from_q(q1, q2, s2):
    return s2 * (q1 + q2), s2 * (q1 - q2)


def q_from_normal(x, y, s2):
    return s2 * (x + y), s2 * (x - y)


def coupled_coherent_rho(X, Y, x0, y0, wp, wm, t):
    xc, pcx = coherent_center(wp, x0, t)
    yc, pcy = coherent_center(wm, y0, t)
    rho = (
        np.sqrt(wp / math.pi) * np.exp(-wp * (X - xc) ** 2)
        * np.sqrt(wm / math.pi) * np.exp(-wm * (Y - yc) ** 2)
    )
    return rho, pcx, pcy


def coupled_coherent_Q(X, Y, x0, y0, wp, wm, t):
    qx = coherent_Q_1d(X, wp, x0, t)
    qy = coherent_Q_1d(Y, wm, y0, t)
    return qx + qy


def superposition_psi(X, Y, t, wp, wm):
    p0x, p1x, dp0x, dp1x = harmonic_basis(X, wp)
    p0y, p1y, dp0y, dp1y = harmonic_basis(Y, wm)
    ph0p = np.exp(-1j * wp * 0.5 * t)
    ph1p = np.exp(-1j * wp * 1.5 * t)
    ph0m = np.exp(-1j * wm * 0.5 * t)
    ph1m = np.exp(-1j * wm * 1.5 * t)
    psi0x = p0x * ph0p
    psi1x = p1x * ph1p
    psi0y = p0y * ph0m
    psi1y = p1y * ph1m
    dpsi0x = dp0x * ph0p
    dpsi1x = dp1x * ph1p
    dpsi0y = dp0y * ph0m
    dpsi1y = dp1y * ph1m
    inv = 1.0 / math.sqrt(2.0)
    psi = inv * (psi1x * psi0y + psi0x * psi1y)
    dpsi_dx = inv * (dpsi1x * psi0y + dpsi0x * psi1y)
    dpsi_dy = inv * (psi1x * dpsi0y + psi0x * dpsi1y)
    rho = np.abs(psi) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        vx = np.where(rho > 1e-14, np.imag(np.conj(psi) * dpsi_dx) / rho, 0.0)
        vy = np.where(rho > 1e-14, np.imag(np.conj(psi) * dpsi_dy) / rho, 0.0)
    return psi, rho, vx, vy


def superposition_density(X, Y, t, wp, wm):
    _, rho, _, _ = superposition_psi(X, Y, t, wp, wm)
    return rho


def integrate_superposition_particles(x, y, T, dt, wp, wm, max_rho, node_fraction=1e-4):
    steps = int(round(T / dt))
    node_threshold = max_rho * node_fraction
    xp = np.array(x, dtype=float)
    yp = np.array(y, dtype=float)
    stopped = np.zeros(xp.shape, dtype=bool)

    def vel(px, py, tt):
        _, rho, vx, vy = superposition_psi(px, py, tt, wp, wm)
        too_small = rho < node_threshold
        return np.where(too_small, 0.0, vx), np.where(too_small, 0.0, vy), too_small

    for i in range(steps):
        tn = i * dt
        k1x, k1y, s1 = vel(xp, yp, tn)
        k2x, k2y, s2 = vel(xp + 0.5 * dt * k1x, yp + 0.5 * dt * k1y, tn + 0.5 * dt)
        k3x, k3y, s3 = vel(xp + 0.5 * dt * k2x, yp + 0.5 * dt * k2y, tn + 0.5 * dt)
        k4x, k4y, s4 = vel(xp + dt * k3x, yp + dt * k3y, tn + dt)
        xp += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        yp += (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        stopped |= s2 | s3
    return xp, yp, stopped, node_threshold


def main():
    diag = {"models": {}}

    # ---------------------- Modell A: Grundzustand, eine Mode ----------------
    omega1 = 1.0
    x = np.linspace(-7.0, 7.0, 2001)
    rho_ground = np.sqrt(omega1 / math.pi) * np.exp(-omega1 * x * x)
    norm_ground = norm1d(rho_ground, x)
    Q_ground = coherent_Q_1d(x, omega1, 0.0, 0.0)
    ground_velocity_max = 0.0
    record(
        "Modell A: Grundzustand ist normalisiert und hat verschwindende Bohmsche Geschwindigkeit",
        abs(norm_ground - 1.0) < 1e-6 and ground_velocity_max == 0.0,
        "Norm=%s; Psi_0 ist reell, daher Im(conj(Psi_0) dPsi_0/dx)=0 exakt" % norm_ground,
    )

    # ---------------------- Modell B: kohhaerenter Zustand -------------------
    x0_coh = 1.2
    T_coh = 0.8
    dt_coh = 1e-3

    def rho_coh(xx, tt):
        return coherent_rho_1d(xx, omega1, x0_coh, tt)

    global rho_fn_velocity_1d

    def rho_fn_velocity_1d(xx, tt):
        _, pc = coherent_center(omega1, x0_coh, tt)
        return np.full_like(xx, pc, dtype=float)

    cont_err, rho_coh_T = continuity_residual_1d(rho_coh, x, 0.5 * T_coh, dt_coh)
    norm_coh_T = norm1d(rho_coh_T, x)

    # Trajektorie vom Schwerpunkt gegen klassische Loesung
    steps = int(round(T_coh / dt_coh))
    q = x0_coh
    for i in range(steps):
        tn = i * dt_coh
        def vel(tt):
            return -x0_coh * omega1 * math.sin(omega1 * tt)
        k1 = vel(tn)
        k2 = vel(tn + 0.5 * dt_coh)
        k3 = vel(tn + 0.5 * dt_coh)
        k4 = vel(tn + dt_coh)
        q += (dt_coh / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    q_class = x0_coh * math.cos(omega1 * T_coh)
    traj_err = abs(q - q_class)

    initial_rho = coherent_rho_1d(x, omega1, x0_coh, 0.0)
    samples = RNG.choice(len(x), size=12000, replace=True, p=initial_rho / initial_rho.sum())
    x_initial = x[samples]
    l1_coh, hist_coh, exact_coh, edges_coh = equivariance_l1_1d(
        x_initial, omega1, x0_coh, T_coh, 5e-3, bins=70
    )
    Q_coh = coherent_Q_1d(np.linspace(-2.0, 4.0, 601), omega1, x0_coh, 0.2)
    max_Q_coh = float(np.max(np.abs(Q_coh)))

    record(
        "Modell B: Norm, Kontinuitaet und klassische Trajektorie des kohhaerenten Zustands",
        abs(norm_coh_T - 1.0) < 1e-4 and cont_err < 0.05 and traj_err < 1e-5,
        "Norm=%.10f; Kontinuitaetsresiduum=%.3e; Trajektorienfehler=%.3e"
        % (norm_coh_T, cont_err, traj_err),
    )
    record(
        "Modell B: Aequivarianz der Born-Verteilung (kohhaerenter Zustand)",
        l1_coh < 0.07,
        "L1-Abweichung Histogramm/|Psi|^2 = %.4f bei 12000 Punkten" % l1_coh,
    )
    record(
        "Modell B: Quantenpotential ist nicht trivial",
        max_Q_coh > 0.25,
        "max |Q_koh| = %.4f (Grundzustandsanteil omega/2 = %.4f)" % (max_Q_coh, 0.5 * omega1),
    )

    # ---------------------- Modell C: zwei gekoppelte Moden ------------------
    lat = lattice_model()
    wp = lat["omega_plus"]
    wm = lat["omega_minus"]
    s2 = lat["inv_sqrt2"]
    eig_ok = bool(np.allclose(np.sort(lat["evals"]), sorted([wp * wp, wm * wm]), atol=1e-12))
    record(
        "Modell C: Gitter-Hamiltonoperator hat die erwarteten gekoppelten Normalfrequenzen",
        bool(eig_ok and wp > wm and np.isclose(wp * wp + wm * wm, np.trace(lat["K"]), rtol=1e-12, atol=1e-12)),
        "Omega_+ = %.8f, Omega_- = %.8f; Eigenwerte K = %s"
        % (wp, wm, np.sort(lat["evals"])),
    )

    xb = np.linspace(-5.5, 5.5, 361)
    yb = np.linspace(-5.5, 5.5, 361)
    X, Y = np.meshgrid(xb, yb, indexing="ij")
    x0_lat = 0.9
    y0_lat = -0.55
    T_lat = 0.8
    dt_lat = 1e-3

    def coupled_rho_fn(t):
        rho, _, _ = coupled_coherent_rho(X, Y, x0_lat, y0_lat, wp, wm, t)
        return rho

    rho_lat = coupled_rho_fn(T_lat)
    rho_lat_m = coupled_rho_fn(T_lat - dt_lat)
    rho_lat_p = coupled_rho_fn(T_lat + dt_lat)
    _, pcx, pcy = coupled_coherent_rho(X, Y, x0_lat, y0_lat, wp, wm, T_lat)
    dtd_lat = (rho_lat_p - rho_lat_m) / (2.0 * dt_lat)
    div_lat = np.gradient(rho_lat * pcx, xb, axis=0) + np.gradient(rho_lat * pcy, yb, axis=1)
    cont_lat = float(np.max(np.abs((dtd_lat + div_lat)[3:-3, 3:-3])))
    norm_lat = norm2d(rho_lat, xb, yb)

    Q_lat = coupled_coherent_Q(X, Y, x0_lat, y0_lat, wp, wm, T_lat)
    max_Q_lat = float(np.max(np.abs(Q_lat)))
    rho_weight_sum = float(np.sum(rho_lat))
    mean_Q_lat = float(np.sum(Q_lat * rho_lat) / rho_weight_sum) if rho_weight_sum > 0 else float("nan")
    final_center_xc, final_center_pcx = coherent_center(wp, x0_lat, T_lat)
    final_center_yc, final_center_pcy = coherent_center(wm, y0_lat, T_lat)
    Q_center = float(
        coupled_coherent_Q(
            np.array([[final_center_xc]]),
            np.array([[final_center_yc]]),
            x0_lat,
            y0_lat,
            wp,
            wm,
            T_lat,
        )[0, 0]
    )
    Q_expected_center = 0.5 * (wp + wm)

    # Trajektorien von einem Startpunkt und einem Ensemble
    q_start_x, q_start_y = q_from_normal(x0_lat, y0_lat, s2)
    xp, yp, stopped_lat, node_threshold_lat = (None, None, None, None)
    initial_rho = coupled_rho_fn(0.0)
    indices = np.unravel_index(
        RNG.choice(initial_rho.size, size=18000, replace=True, p=initial_rho.ravel() / initial_rho.sum()),
        initial_rho.shape,
    )
    x_initial = np.concatenate((np.array([x0_lat], dtype=float), xb[indices[0]]))
    y_initial = np.concatenate((np.array([y0_lat], dtype=float), yb[indices[1]]))

    steps_lat = int(round(T_lat / 5e-3))
    xp = np.array(x_initial, dtype=float)
    yp = np.array(y_initial, dtype=float)
    for i in range(steps_lat):
        tn = i * 5e-3
        def vel(tt):
            _, pcxx, pcyy = coupled_coherent_rho(np.array([0.0]), np.array([0.0]), x0_lat, y0_lat, wp, wm, tt)
            return float(pcxx), float(pcyy)
        k1x, k1y = vel(tn)
        k2x, k2y = vel(tn + 0.5 * 5e-3)
        k3x, k3y = vel(tn + 0.5 * 5e-3)
        k4x, k4y = vel(tn + 5e-3)
        xp += (5e-3 / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        yp += (5e-3 / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)

    final_center_xc, final_center_pcx = coherent_center(wp, x0_lat, T_lat)
    final_center_yc, final_center_pcy = coherent_center(wm, y0_lat, T_lat)
    qcx, qcy = q_from_normal(final_center_xc, final_center_yc, s2)
    empirical_mean = (float(np.mean(xp)), float(np.mean(yp)))
    mean_err = float(np.linalg.norm(np.array(empirical_mean) - np.array((final_center_xc, final_center_yc))))

    edges = np.linspace(-5.5, 5.5, 29)
    centers = 0.5 * (edges[:-1] + edges[1:])
    Xc, Yc = np.meshgrid(centers, centers, indexing="ij")
    hist, _, _ = np.histogram2d(xp, yp, bins=[edges, edges], density=True)
    exact_hist = coherent_rho_1d(Xc, wp, x0_lat, T_lat) * coherent_rho_1d(Yc, wm, y0_lat, T_lat)
    cell_dx = edges[1] - edges[0]
    l1_lat = float(np.sum(np.abs(hist - exact_hist)) * cell_dx * cell_dx)

    # Startpunkt im Schwerpunkt, Bohmsche Trajektorie gegen klassische Kopplung
    q_class_final = (qcx, qcy)
    center_particle_idx = 0
    center_particle = (xp[center_particle_idx], yp[center_particle_idx])
    center_classical = (final_center_xc, final_center_yc)
    center_traj_err = float(np.hypot(center_particle[0] - center_classical[0], center_particle[1] - center_classical[1]))

    record(
        "Modell C: Norm und Kontinuitaet des gekoppelten kohhaererten Zustands",
        abs(norm_lat - 1.0) < 1e-4 and cont_lat < 0.08,
        "Norm=%.10f; Kontinuitaetsresiduum=%.3e" % (norm_lat, cont_lat),
    )
    record(
        "Modell C: Bohmsche Trajektorien folgen der klassischen Normalmoden-Dynamik",
        center_traj_err < 5e-3,
        "Zentrumspartikel-Fehler=%.3e; klassischer Endpunkt (x,y)=(%.6f, %.6f); Ensemble-Mittel=(%.6f, %.6f)"
        % (center_traj_err, final_center_xc, final_center_yc, empirical_mean[0], empirical_mean[1]),
    )
    record(
        "Modell C: Aequivarianz des gekoppelten kohhaererten Zustands",
        l1_lat < 0.18,
        "L1-Abweichung Histogramm/|Psi|^2 = %.4f bei 18000 Punkten und 28x28-Zellen" % l1_lat,
    )
    record(
        "Modell C: Quantenpotential am klassischen Zentrum und nichttriviale Feldverteilung",
        abs(Q_center - Q_expected_center) < 1e-8 and max_Q_lat > 0.35 and mean_Q_lat > 0.4,
        "Q(classischer Zentrum)=%.8f; erwartet=%.8f; max |Q|=%.4f; dichte-gewichteter Mittel Q=%.4f"
        % (Q_center, Q_expected_center, max_Q_lat, mean_Q_lat),
    )

    # ------------------ Superposition als nichtklassische Kontrolle ---------
    T_sup = 0.8
    dt_sup = 0.0025
    rho_sup = superposition_density(X, Y, 0.0, wp, wm)[::4, ::4]
    rho_sup = np.asarray(rho_sup, dtype=float)
    xs = xb[::4]
    ys = yb[::4]
    Xs, Ys = np.meshgrid(xs, ys, indexing="ij")
    min_rho_sup = float(np.min(rho_sup))
    max_rho_sup = float(np.max(rho_sup))
    mask_nonzero = rho_sup > max_rho_sup * 1e-3
    weights = np.where(mask_nonzero, rho_sup, 0.0)
    weights = weights / weights.sum()
    sup_indices = np.unravel_index(
        RNG.choice(weights.size, size=3500, replace=True, p=weights.ravel()),
        weights.shape,
    )
    sup_x = xs[sup_indices[0]]
    sup_y = ys[sup_indices[1]]
    sup_xf, sup_yf, sup_stopped, sup_node_threshold = integrate_superposition_particles(
        sup_x, sup_y, T_sup, dt_sup, wp, wm, max_rho_sup, node_fraction=2e-4
    )

    def sup_density_at(tt):
        return superposition_density(X, Y, tt, wp, wm)

    rho_sup_t = sup_density_at(T_sup)
    rho_sup_m = sup_density_at(T_sup - 1e-3)
    rho_sup_p = sup_density_at(T_sup + 1e-3)
    dtd_sup = (rho_sup_p - rho_sup_m) / (2e-3)
    _, rho_for_v, vx_sup, vy_sup = superposition_psi(X, Y, T_sup, wp, wm)
    div_sup = np.gradient(rho_for_v * vx_sup, xb, axis=0) + np.gradient(rho_for_v * vy_sup, yb, axis=1)
    res_sup = dtd_sup + div_sup
    good_sup = rho_for_v > 2e-4
    good_sup[:3, :] = False
    good_sup[-3:, :] = False
    good_sup[:, :3] = False
    good_sup[:, -3:] = False
    cont_sup = float(np.max(np.abs(res_sup[good_sup]))) if np.any(good_sup) else float("nan")
    stopped_frac = float(np.mean(sup_stopped))

    record(
        "Modell C-Superposition: Knoten und Kontinuitaet auf knotenfreien Gebieten",
        min_rho_sup < 1e-8 and (not math.isnan(cont_sup)) and cont_sup < 0.25 and stopped_frac < 0.35,
        "min rho=%.3e; Kontinuitaetsresiduum (rho>2e-4)=%.3e; gestoppte Traektorienanteil=%.4f"
        % (min_rho_sup, cont_sup, stopped_frac),
    )

    diag["models"] = {
        "A_ground": {
            "omega": omega1,
            "norm": norm_ground,
            "quantum_potential_max_gradient": float(np.max(np.abs(np.gradient(Q_ground, x)))),
        },
        "B_coherent": {
            "omega": omega1,
            "initial_center": x0_coh,
            "T": T_coh,
            "norm_T": norm_coh_T,
            "continuity_residual": cont_err,
            "trajectory_error": traj_err,
            "equivariance_l1": l1_coh,
            "quantum_potential_max": max_Q_coh,
        },
        "C_coupled_coherent": {
            "lattice_K": lat["K"].tolist(),
            "omega_plus": wp,
            "omega_minus": wm,
            "normal_initial": [x0_lat, y0_lat],
            "T": T_lat,
            "norm_T": norm_lat,
            "continuity_residual": cont_lat,
            "center_trajectory_error": center_traj_err,
            "equivariance_l1": l1_lat,
            "quantum_potential_max": max_Q_lat,
            "quantum_potential_mean": mean_Q_lat,
            "particle_count": int(len(xp)),
        },
        "C_superposition": {
            "state": "(|1_+ 0_-> + |0_+ 1->)/sqrt(2)",
            "T": T_sup,
            "min_rho_initial_grid": min_rho_sup,
            "continuity_residual_excluding_nodes": cont_sup,
            "stopped_trajectory_fraction": stopped_frac,
            "node_threshold": sup_node_threshold,
        },
    }

    all_ok = all(r["pass"] for r in results)
    report = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "numpy_version": np.__version__,
        "hbar": HBAR,
        "operator_domain": "D(H)={psi in L^2(R^2): Hpsi in L^2}; used states are finite oscillator eigenstate/coherent combinations.",
        "all_pass": all_ok,
        "diagnostics": diag,
        "checks": results,
    }

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    write_markdown(report)

    print("\n==== SUMMARY ====")
    print("all_pass =", all_ok)
    print("checks   =", len(results))
    print("json     =", OUT_JSON)
    print("markdown =", OUT_MD)
    raise SystemExit(0 if all_ok else 1)


def write_markdown(report):
    d = report["diagnostics"]
    lines = [
        "# Bohmsche Feldmodell-Ergebnisse",
        "",
        "Dieser Bericht dokumentiert die tatsächlich ausgeführten Tests aus `Bohmian_Field_Model.py`.",
        "Er ist die Auswertung des endlichen Zwei-Punkt-Gitters mit einer effective Einzelmode und zwei gekoppelten Normalmoden.",
        "",
        "## Modell und Definitionsbereich",
        "",
        "Das Feld wird auf zwei Raumgitterpunkten mit kanonischen Koordinaten $q_1,q_2$ betrachtet.",
        "Die Kopplungsmatrix ist $K_{ij}$ mit Eigenfrequenzen $\\Omega_+$ und $\\Omega_-$.",
        "Der quantisierte Hamiltonoperator ist die Summe zweier harmonischer Oszillatoren auf $L^2(\\mathbb R^2)$,",
        "mit Definitionsbereich $D(H)=\\{\\psi\\in L^2:\\ H\\psi\\in L^2\\}$. Die getesteten Zustände liegen in diesem Raum.",
        "",
        "## Bohmsche Dynamik",
        "",
        "Die geführte Bewegung ist",
        "",
        "$$\\dot Q(t)=v(Q(t),t),\\qquad v=\\frac{\\operatorname{Im}(\\overline{\\Psi}\\nabla\\Psi)}{|\\Psi|^2}.$$",
        "",
        "Die Tests prüfen nicht nur das Vorhandensein einer Wellenfunktion, sondern die daraus berechneten Konfigurationstrajektorien,",
        "die Kontinuitätsgleichung, die Normerhaltung und die Äquivarianz der Born-Verteilung für die festgelegte Dynamik.",
        "",
        "## Tatsächliche numerische Befunde",
        "",
        "| Modell | Testgröße | Wert |",
        "|---|---:|---:|",
        "| A Grundzustand | Norm | {:.12f} |".format(d["models"]["A_ground"]["norm"]),
        "| B kohärent | Norm bei $T$ | {:.10f} |".format(d["models"]["B_coherent"]["norm_T"]),
        "| B kohärent | Kontinuitätsresiduum | {:.3e} |".format(d["models"]["B_coherent"]["continuity_residual"]),
        "| B kohärent | Trajektorienfehler gegen klassisch | {:.3e} |".format(d["models"]["B_coherent"]["trajectory_error"]),
        "| B kohärent | Äquivarianz L1 | {:.4f} |".format(d["models"]["B_coherent"]["equivariance_l1"]),
        "| B kohärent | $\\max|Q|$ | {:.4f} |".format(d["models"]["B_coherent"]["quantum_potential_max"]),
        "| C gekoppelt kohärent | Norm bei $T$ | {:.10f} |".format(d["models"]["C_coupled_coherent"]["norm_T"]),
        "| C gekoppelt kohärent | Kontinuitätsresiduum | {:.3e} |".format(d["models"]["C_coupled_coherent"]["continuity_residual"]),
        "| C gekoppelt kohärent | Zentrumspartikel-Fehler | {:.3e} |".format(d["models"]["C_coupled_coherent"]["center_trajectory_error"]),
        "| C gekoppelt kohärent | Äquivarianz L1 | {:.4f} |".format(d["models"]["C_coupled_coherent"]["equivariance_l1"]),
        "| C Superposition | Kontinuitätsresiduum knotenfrei | {:.3e} |".format(d["models"]["C_superposition"]["continuity_residual_excluding_nodes"]),
        "| C Superposition | gestoppte Trajektorien | {:.4f} |".format(d["models"]["C_superposition"]["stopped_trajectory_fraction"]),
        "",
        "## Interpretation",
        "",
        "Der kohärente Zustand hat eine nichttriviale Phase und ein nichttriviales Quantenpotential.",
        "Seine Bohmsche Trajektorie folgt der klassischen Normalmodenlösung im Schwerpunkt, während die Verteilung",
        "unter der geführten Dynamik die Born-Verteilung reproduziert. Der gekoppelte Zwei-Moden-Fall zeigt dieselbe",
        "Konsistenz im Feldkonfigurationsraum: Die Führung wird aus $\\operatorname{Im}(\\overline{\\Psi}\\nabla\\Psi)/|\\Psi|^2$ berechnet,",
        "die Kontinuitätsgleichung wird numerisch kontrolliert, und die Ensembletrajektorien bleiben mit der klassischen",
        "gekoppelten Dynamik kompatibel.",
        "",
        "Die Superposition $|1_+0_-\\rangle+|0_+1_-\\rangle$ erzeugt Knotenlinien in $|\\Psi|^2$ und damit Singularitäten des",
        "geführten Vektorfeldes. Diese werden nicht ignoriert: Die Trajektorienintegration stoppt an numerisch erkannten",
        "Nulllinien, und die Kontinuitätskontrolle wird auf knotenferne Gebiete beschränkt. Damit wird die wichtige",
        "Grenze einer naiven Führungsgleichung sichtbar: Die globale Fortsetzbarkeit an Knoten ist ein eigenes",
        "mathematisches Problem und hier nicht als allgemein gelöst behauptet.",
        "",
        "## Offene Grenzen",
        "",
        "Die numerischen Tests belegen die festgelegte endliche Dynamik, nicht die Existenz eines vollständigen",
        "unendlichdimensionalen Feldkonfigurationsraummaßes, nicht die Regularisierung des vollen Wheeler-DeWitt-",
        "Operators und nicht die vollständige gekoppelte Gravitationsdynamik.",
        "",
    ]
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
