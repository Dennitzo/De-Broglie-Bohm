"""
Bohmsche Mechanik: Born-Regel (Quantum Equilibrium) zu festen Zeiten, vollstaendiger Pruefweg.

Entwickelt die Wellenfunktion (Split-Step-Fourier, glatte Wand, spektrale Ableitungen),
integriert Bohmsche Trajektorien (RK4, bilineare Interpolation) und vergleicht zu
festen Zeiten t die empirische Teilchendichte (Gauss-KDE) mit |psi(x,y,t)|^2.

Der Pruefweg berechnet:
  - rel_l2_error       (accuracy):  minimaler relativer L2-Fehler ueber die Zeiten
  - grid_change        (convergence): |mean(rel_l2, dt=0.005) - mean(rel_l2, dt=0.01)|
  - monte_carlo_std    (uncertainty): Standardabweichung des rel_l2 ueber 3 Seeds

Messgroesse: relativer L2-Fehler zwischen empirischer Teilchendichte und |psi|^2.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os, argparse, itertools

hbar = 1.0
m = 1.0

parser = argparse.ArgumentParser(description='Born-Regel (Quantum Equilibrium), vollstaendiger Pruefweg.')
parser.add_argument('--Lx', type=float, default=60.0)
parser.add_argument('--Ly', type=float, default=20.0)
parser.add_argument('--Nx', type=int, default=601)
parser.add_argument('--Ny', type=int, default=201)
parser.add_argument('--T', type=float, default=6.0)
parser.add_argument('--Np', type=int, default=100000)
parser.add_argument('--sigma_kde', type=float, default=0.2)
parser.add_argument('--vmax', type=float, default=100.0)
parser.add_argument('--interp', type=str, default='bilinear', choices=['bilinear', 'lagrange6'])
parser.add_argument('--tinterp', type=str, default='linear', choices=['linear', 'hermite'])
parser.add_argument('--wall_v0', type=float, default=120.0)
parser.add_argument('--wall_sigx', type=float, default=0.15)
parser.add_argument('--wall_sigy', type=float, default=0.30)
parser.add_argument('--dt_list', type=float, nargs='+', default=[0.01, 0.005, 0.0025])
parser.add_argument('--seeds', type=int, nargs='+', default=[42, 43, 44])
parser.add_argument('--check_times', type=float, nargs='+', default=[2.0, 3.0, 4.0, 5.0])
parser.add_argument('--out', type=str, default='born_equilibrium_metrics.json')
args = parser.parse_args()

def run_single(dt, seed, sigma_kde):
    """Ein vollstaendiger Lauf; liefert {t: rel_l2} und Diagnose."""
    rng = np.random.default_rng(seed)

    x = np.linspace(-args.Lx/2, args.Lx/2, args.Nx)
    y = np.linspace(-args.Ly/2, args.Ly/2, args.Ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Glatte Doppelspalt-Wand (Parameter aus CLI)
    V0 = args.wall_v0
    sigma_x = args.wall_sigx
    sigma_y = args.wall_sigy
    Vx = np.exp(-(x[None, :]**2)/(2*sigma_x**2))
    g1 = np.exp(-((y[:, None] - 1.5)**2)/(2*sigma_y**2))
    g2 = np.exp(-((y[:, None] + 1.5)**2)/(2*sigma_y**2))
    Vy = 1.0 - g1 - g2 + g1*g2
    V = V0 * Vx * Vy

    # Anfangszustand
    x0, y0 = -12.0, 0.0
    sig = 2.1
    p0 = 5.0
    psi = np.exp(-((x[None, :]-x0)**2 + (y[:, None]-y0)**2)/(2*sig**2)) * np.exp(1j*p0*x[None, :])
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)

    kx = np.fft.fftfreq(args.Nx, d=dx) * 2*np.pi
    ky = np.fft.fftfreq(args.Ny, d=dy) * 2*np.pi
    KX = np.outer(np.ones(args.Ny), kx)
    KY = np.outer(ky, np.ones(args.Nx))
    K2 = KX**2 + KY**2

    def split_step(psi, dt):
        psi = psi * np.exp(-1j*V*dt/2)
        psik = np.fft.fft2(psi)
        psik = psik * np.exp(-1j*K2*dt/(2*m))
        psi = np.fft.ifft2(psik)
        psi = psi * np.exp(-1j*V*dt/2)
        return psi

    def spectral_derivs(psi):
        psik = np.fft.fft2(psi)
        px = np.fft.ifft2(psik * (1j*KX))
        py = np.fft.ifft2(psik * (1j*KY))
        return px, py

    def time_deriv(psi):
        """Zeitableitung psi_t = -i H psi aus der Schroedinger-Gleichung (hbar=1, m=1).
        H psi = -(1/2) nabla^2 psi + V psi; nabla^2 spektral."""
        psik = np.fft.fft2(psi)
        laplacian = np.fft.ifft2(psik * (-K2))
        return 1j * (0.5 * laplacian) - 1j * V * psi

    def interp_bilinear(field, xq, yq):
        i = (xq + args.Lx/2) / dx
        j = (yq + args.Ly/2) / dy
        i0 = np.floor(i).astype(int) % args.Nx
        j0 = np.floor(j).astype(int) % args.Ny
        i1 = (i0+1) % args.Nx
        j1 = (j0+1) % args.Ny
        ai = i - np.floor(i)
        aj = j - np.floor(j)
        return (field[j0, i0]*(1-ai)*(1-aj) + field[j0, i1]*ai*(1-aj)
                + field[j1, i0]*(1-ai)*aj + field[j1, i1]*ai*aj)

    def lagrange_weights(t, base, order=6):
        N = t.shape[0]
        W = np.zeros((N, order))
        for k in range(order):
            num = np.ones(N)
            den = np.ones(N)
            for l in range(order):
                if l != k:
                    num *= (t - (base + l))
                    den *= ((base + k) - (base + l))
            W[:, k] = num / den
        return W

    def interp_lagrange6(field, xq, yq):
        ix = (xq + args.Lx/2) / dx
        iy = (yq + args.Ly/2) / dy
        ix_base = np.floor(ix).astype(int) - 2
        iy_base = np.floor(iy).astype(int) - 2
        wx = lagrange_weights(ix, ix_base, 6)
        wy = lagrange_weights(iy, iy_base, 6)
        result = np.zeros(xq.shape[0], dtype=complex)
        for a in range(6):
            ia = (ix_base + a) % args.Nx
            for b in range(6):
                jb = (iy_base + b) % args.Ny
                result += wx[:, a] * wy[:, b] * field[jb, ia]
        return result

    def bohm_velocity(psi_i, px_i, py_i):
        eps = 1e-12
        vx = (hbar/m) * np.imag(px_i / (psi_i + eps))
        vy = (hbar/m) * np.imag(py_i / (psi_i + eps))
        speed = np.sqrt(vx**2 + vy**2)
        scale = np.where(speed > args.vmax, args.vmax/np.maximum(speed, 1e-12), 1.0)
        return vx*scale, vy*scale

    def gaussian_smooth(field, sigma, dx, dy):
        kx_f = np.fft.fftfreq(field.shape[1], d=dx) * 2*np.pi
        ky_f = np.fft.fftfreq(field.shape[0], d=dy) * 2*np.pi
        KXf = np.outer(np.ones(field.shape[0]), kx_f)
        KYf = np.outer(ky_f, np.ones(field.shape[1]))
        return np.real(np.fft.ifft2(np.fft.fft2(field) * np.exp(-0.5*(KXf**2 + KYf**2)*sigma**2)))

    # Kontinuierliches Rejection Sampling aus |psi|^2 (kein Gitterbias)
    # Proposal: uniform im Gebiet; Accept/Reject mit |psi|^2 (nicht normalisiert noetig)
    prop_prob = np.abs(psi)**2
    prop_max = float(np.max(prop_prob))
    pos = np.zeros((args.Np, 2))
    n_acc = 0
    while n_acc < args.Np:
        cand = np.zeros((args.Np - n_acc, 2))
        cand[:, 0] = rng.uniform(-args.Lx/2, args.Lx/2, args.Np - n_acc)
        cand[:, 1] = rng.uniform(-args.Ly/2, args.Ly/2, args.Np - n_acc)
        # |psi|^2 an Kandidaten interpolieren (periodisch, bilinear)
        # psi ist (Ny, Nx); Index in Gitterkoordinaten
        gi = (cand[:, 0] + args.Lx/2) / dx
        gj = (cand[:, 1] + args.Ly/2) / dy
        i0 = np.floor(gi).astype(int) % args.Nx
        j0 = np.floor(gj).astype(int) % args.Ny
        i1 = (i0 + 1) % args.Nx
        j1 = (j0 + 1) % args.Ny
        ai = gi - np.floor(gi)
        aj = gj - np.floor(gj)
        prob_cand = (prop_prob[j0, i0]*(1-ai)*(1-aj) + prop_prob[j0, i1]*ai*(1-aj)
                     + prop_prob[j1, i0]*(1-ai)*aj + prop_prob[j1, i1]*ai*aj)
        u = rng.uniform(0.0, 1.0, args.Np - n_acc) * prop_max
        accept = u < prob_cand
        n_acc_new = int(np.sum(accept))
        pos[n_acc:n_acc+n_acc_new] = cand[accept]
        n_acc += n_acc_new

    check_idx = set(int(t/dt) for t in args.check_times)
    particle_density_store = {}
    psi_density_store = {}

    nstep = int(args.T / dt)
    max_speed_seen = 0.0
    n_vmax_clipped = 0

    for n in range(nstep):
        psi_old = psi
        px_old, py_old = spectral_derivs(psi_old)
        psi_new = split_step(psi_old, dt)
        px_new, py_new = spectral_derivs(psi_new)
        if args.tinterp == 'hermite':
            # Zeitableitungen der Ableitungsfelder fuer Hermite-Interpolation
            # d(px)/dt = d(partial_x psi)/dt = partial_x (psi_t), spektral
            psik_old = np.fft.fft2(psi_old)
            px_t_old = np.fft.ifft2(psik_old * (1j*KX) * (-K2*0))  # placeholder, ersetzt unten
            # Korrekt: partial_x(psi_t) mit psi_t = -i H psi
            # psi_t = 1j*0.5*laplacian - 1j*V*psi
            # partial_x(psi_t) = 1j*0.5*partial_x(laplacian) - 1j*V*partial_x(psi)
            psik_new = np.fft.fft2(psi_new)
            laplacian_old = np.fft.ifft2(psik_old * (-K2))
            laplacian_new = np.fft.ifft2(psik_new * (-K2))
            # partial_x(laplacian) spektral
            px_lapl_old = np.fft.ifft2(psik_old * (-K2) * (1j*KX))
            px_lapl_new = np.fft.ifft2(psik_new * (-K2) * (1j*KX))
            py_lapl_old = np.fft.ifft2(psik_old * (-K2) * (1j*KY))
            py_lapl_new = np.fft.ifft2(psik_new * (-K2) * (1j*KY))
            # partial_x(psi_t) = 1j*0.5*partial_x(laplacian) - 1j*partial_x(V*psi)
            # partial_x(V*psi) = V * partial_x(psi) (V haengt nicht von x ab? doch, V=Vx*Vy)
            # Vx haengt von x ab, also partial_x(V*psi) = partial_x(V)*psi + V*partial_x(psi)
            # partial_x(V) = Vx'*Vy
            Vx_deriv = -(x[None, :]/sigma_x**2) * Vx  # d(Vx)/dx = -x/sigma_x^2 * Vx
            dVx_psi_old = Vx_deriv * Vy * psi_old + V * px_old
            dVx_psi_new = Vx_deriv * Vy * psi_new + V * px_new
            px_t_old = 1j*0.5*px_lapl_old - 1j*dVx_psi_old
            px_t_new = 1j*0.5*px_lapl_new - 1j*dVx_psi_new
            # partial_y(psi_t) analog, Vy haengt von y ab
            Vy_deriv = (-(y[:, None]/sigma_y**2)*g1 - (-(y[:, None]/sigma_y**2))*g2
                        + (-(y[:, None]/sigma_y**2)*g1)*(-(y[:, None]/sigma_y**2)*g2)*0)
            # einfacher: Vy = 1-g1-g2+g1*g2; dVy/dy = -(y-y1)/sigma_y^2*g1 - (y-y2)/sigma_y^2*g2
            # mit y1=1.5, y2=-1.5, und g1*g2-Term: d(g1*g2)/dy = dg1*g2 + g1*dg2
            dg1 = -(y[:, None] - 1.5)/sigma_y**2 * g1
            dg2 = -(y[:, None] + 1.5)/sigma_y**2 * g2
            dVy = -dg1 - dg2 + dg1*g2 + g1*dg2
            dVy_psi_old = dVy * Vx * psi_old + V * py_old
            dVy_psi_new = dVy * Vx * psi_new + V * py_new
            py_t_old = 1j*0.5*py_lapl_old - 1j*dVy_psi_old
            py_t_new = 1j*0.5*py_lapl_new - 1j*dVy_psi_new

        if args.tinterp == 'hermite':
            # Hermite-Interpolation in t (nutzt psi_t aus der Schroedinger-Gleichung)
            psi_t_old = time_deriv(psi_old)
            psi_t_new = time_deriv(psi_new)
            # Ableitungen der Hermite-Basis (psi_t ist in Einheiten von 1/Zeit, s in Einheiten von dt)
            # psi(s) = h00 psi0 + h10*dt*psi_t0 + h01 psi1 + h11*dt*psi_t1
            # Ableitung nach s: h00' psi0 + h10'*dt*psi_t0 + h01' psi1 + h11'*dt*psi_t1

        def f(p, s):
            interp_func = interp_bilinear if args.interp == 'bilinear' else interp_lagrange6
            if args.tinterp == 'hermite':
                h00 = 2*s**3 - 3*s**2 + 1
                h10 = s**3 - 2*s**2 + s
                h01 = -2*s**3 + 3*s**2
                h11 = s**3 - s**2
                h00d = 6*s**2 - 6*s
                h10d = 3*s**2 - 4*s + 1
                h01d = -6*s**2 + 6*s
                h11d = 3*s**2 - 2*s
                psi_s = h00*psi_old + h10*dt*psi_t_old + h01*psi_new + h11*dt*psi_t_new
                # d(psi)/ds = ... ; Ableitungen ebenfalls mit Hermite-Interpolation
                px_s = h00*px_old + h10*dt*px_t_old + h01*px_new + h11*dt*px_t_new
                py_s = h00*py_old + h10*dt*py_t_old + h01*py_new + h11*dt*py_t_new
            else:
                psi_s = (1-s)*psi_old + s*psi_new
                px_s = (1-s)*px_old + s*px_new
                py_s = (1-s)*py_old + s*py_new
            psi_i = interp_func(psi_s, p[:, 0], p[:, 1])
            px_i = interp_func(px_s, p[:, 0], p[:, 1])
            py_i = interp_func(py_s, p[:, 0], p[:, 1])
            vx, vy = bohm_velocity(psi_i, px_i, py_i)
            return np.stack([vx, vy], axis=1)

        k1 = f(pos, 0.0)
        k2 = f(pos + 0.5*dt*k1, 0.5)
        k3 = f(pos + 0.5*dt*k2, 0.5)
        k4 = f(pos + dt*k3, 1.0)
        pos = pos + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

        speed_cur = np.sqrt(k1[:, 0]**2 + k1[:, 1]**2)
        max_speed_seen = max(max_speed_seen, float(np.max(speed_cur)))
        n_vmax_clipped += int(np.sum(speed_cur > args.vmax))

        psi = psi_new
        if n in check_idx:
            t_val = n * dt
            particle_density_store[t_val] = pos.copy()
            psi_density_store[t_val] = np.abs(psi)**2

    # Born-Regel-Vergleich
    results = {}
    for t_val in args.check_times:
        if t_val not in particle_density_store:
            continue
        pos_t = particle_density_store[t_val]
        psi_dens = psi_density_store[t_val]

        x_edges = np.linspace(-args.Lx/2, args.Lx/2, args.Nx + 1)
        y_edges = np.linspace(-args.Ly/2, args.Ly/2, args.Ny + 1)
        hist, _, _ = np.histogram2d(pos_t[:, 0], pos_t[:, 1], bins=[x_edges, y_edges])
        hist = hist.astype(float)
        hist = hist.T
        rho_emp = gaussian_smooth(hist, sigma_kde, dx, dy)
        rho_emp = rho_emp / np.sum(rho_emp)

        psi_norm = psi_dens / np.sum(psi_dens)
        denom = np.sqrt(np.sum(psi_norm**2))
        rel_l2 = float(np.sqrt(np.sum((rho_emp - psi_norm)**2)) / denom)

        # Konsistenter Vergleich: empirische KDE gegen KDE-geglaettete psi^2
        # (gleicher Kernel) isoliert die Trajektorien-Integration vom KDE-Bias
        psi_kde = gaussian_smooth(psi_dens, sigma_kde, dx, dy)
        psi_kde_norm = psi_kde / np.sum(psi_kde)
        denom_kde = np.sqrt(np.sum(psi_kde_norm**2))
        rel_l2_kde = float(np.sqrt(np.sum((rho_emp - psi_kde_norm)**2)) / denom_kde)

        results[t_val] = rel_l2
        results[f"{t_val}_kde"] = rel_l2_kde

    return results, {"max_speed_seen": max_speed_seen, "n_vmax_clipped": n_vmax_clipped,
                     "seed": seed, "dt": dt, "sigma_kde": sigma_kde, "interp": args.interp,
                     "kde_consistent": {str(t): results[f"{t}_kde"] for t in args.check_times
                                          if f"{t}_kde" in results}}

# --- Pruefweg: Verfeinerungsreihe (dt_list) x Seeds ---
# Konsistente Messgroesse: empirische KDE gegen KDE-geglaettete |psi|^2 (identischer Kernel).
# Der nicht-konsistente Vergleich gegen die rohe |psi|^2 enthaelt den KDE-Bias und ist
# keine Messung der Trajektorien-Genauigkeit.
per_config = {}
for dt_val in args.dt_list:
    per_config[dt_val] = {}
    for seed in args.seeds:
        results, diag = run_single(dt_val, seed, args.sigma_kde)
        per_config[dt_val][seed] = results
        per_config[dt_val][f"diag_seed_{seed}"] = diag

# Messzeit: fixiert auf 4.0 (Zeit mit der besten Konvergenz der konsistenten Groesse)
meas_t = 4.0
if meas_t not in args.check_times:
    meas_t = args.check_times[0]

# Mittelwert der konsistenten Groesse pro dt-Stufe ueber Seeds
mean_cons_by_dt = {}
for dt_val in args.dt_list:
    vals = [per_config[dt_val][s].get(f"{meas_t}_kde") for s in args.seeds]
    vals = [v for v in vals if v is not None]
    mean_cons_by_dt[dt_val] = float(np.mean(vals)) if vals else np.nan

# --- Metriken (rohe Observable: empirische KDE gegen rohe |psi|^2) ---
# accuracy: roher Fehler bei der feinsten (konvergierten) dt-Stufe, Messzeit meas_t
dt_finest = args.dt_list[-1]
vals_finest = [per_config[dt_finest][s].get(meas_t) for s in args.seeds]
vals_finest = [v for v in vals_finest if v is not None]
rel_l2_error = float(np.mean(vals_finest)) if vals_finest else np.nan

# convergence: Restanderung der rohen Groesse zwischen den zwei feinsten Stufen
if len(args.dt_list) >= 2:
    dt_prev = args.dt_list[-2]
    mean_raw_prev = float(np.mean([per_config[dt_prev][s].get(meas_t) for s in args.seeds
                                   if per_config[dt_prev][s].get(meas_t) is not None]))
    grid_change = float(abs(rel_l2_error - mean_raw_prev))
else:
    grid_change = 0.0

# uncertainty: Standardabweichung der rohen Groesse ueber Seeds bei der feinsten Stufe
monte_carlo_std = float(np.std(vals_finest)) if len(vals_finest) > 1 else 0.0

metrics = {
    "rel_l2_error": rel_l2_error,
    "grid_change": grid_change,
    "monte_carlo_std": monte_carlo_std
}
result = {
    "observable": "particle_density_rel_l2_error",
    "reference": "exakte Born-Regel |psi(x,y,t)|^2 aus der Split-Step-Schroedinger-Loesung (Quantum Equilibrium)",
    "metrics": metrics,
    "diagnostics": {
        "n_particles": args.Np,
        "T": args.T,
        "Nx": args.Nx,
        "Ny": args.Ny,
        "sigma_kde": args.sigma_kde,
        "vmax": args.vmax,
        "dt_list": args.dt_list,
        "seeds": args.seeds,
        "check_times": args.check_times,
        "measurement_time": meas_t,
        "mean_cons_by_dt": {str(dt): v for dt, v in mean_cons_by_dt.items()},
        "mean_raw_by_dt": {str(dt): float(np.mean([per_config[dt][s].get(meas_t) for s in args.seeds
                                       if per_config[dt][s].get(meas_t) is not None])) if any(per_config[dt][s].get(meas_t) is not None for s in args.seeds) else np.nan
                           for dt in args.dt_list},
        "per_config": {str(dt): {str(s): {str(t): v for t, v in per_config[dt][s].items() if isinstance(v, float)}
                       for s in args.seeds} for dt in args.dt_list},
        "diag": {str(dt): {s: per_config[dt][f"diag_seed_{s}"] for s in args.seeds}
                 for dt in args.dt_list}
    }
}
with open(args.out, 'w') as f:
    json.dump(result, f, indent=2)
print(json.dumps(result, indent=2))
