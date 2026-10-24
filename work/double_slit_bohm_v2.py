"""
Bohmsche Mechanik: Doppelspalt-Experiment, korrigierte und validierte Version.

Loest die Schroedinger-Gleichung (Split-Step-Fourier, spektrale Ableitungen),
integriert Bohmsche Trajektorien (RK4, bilineare Interpolation), erfasst die
Erstankunft an der Detektorebene x=det_x mit Schnittpunktinterpolation und
vergleicht die Verteilung mit dem zeitintegrierten positiven Wahrscheinlichkeitsfluss
ueber dieselben Bins. Speichert die Metriken als JSON.

Korrekturen gegenueber double_slit_bohm.py:
  - Spektrale Ableitungen statt zentraler Differenzen (behebt sinc-Fehler).
  - Erstankunft an der tatsaechlichen Detektorebene mit Schnittpunktinterpolation.
  - Fluss mit tatsaechlichem Zeitgewicht dt pro Schritt (nicht dt pro 10-Schritt-Block).
  - Konsistente Bins und Normierung: Histogrammmasse und Flussmasse je pro Bin.
  - Getrennte Pruefung von Trefferwahrscheinlichkeit und Verteilungsform.
  - Grosses Gebiet (Paket nicht abgeschnitten, Rand nicht erreicht).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os, argparse

hbar = 1.0
m = 1.0

parser = argparse.ArgumentParser(description='Doppelspalt, Bohmsche Mechanik, validiert.')
parser.add_argument('--Lx', type=float, default=60.0)
parser.add_argument('--Ly', type=float, default=20.0)
parser.add_argument('--Nx', type=int, default=601)
parser.add_argument('--Ny', type=int, default=201)
parser.add_argument('--dt', type=float, default=0.01)
parser.add_argument('--T', type=float, default=8.0)
parser.add_argument('--Np', type=int, default=2000)
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--det_x', type=float, default=10.0)
parser.add_argument('--n_bins', type=int, default=60)
parser.add_argument('--vmax', type=float, default=30.0)
parser.add_argument('--interp', type=str, default='bilinear', choices=['bilinear', 'lagrange6'])
parser.add_argument('--deriv', type=str, default='spectral', choices=['spectral', 'central'])
parser.add_argument('--wall', type=str, default='smooth', choices=['sharp', 'smooth'])
parser.add_argument('--out', type=str, default='detector_arrival_metrics.json')
args = parser.parse_args()

rng = np.random.default_rng(args.seed)

# --- Gebiet und Gitter ---
x = np.linspace(-args.Lx/2, args.Lx/2, args.Nx)
y = np.linspace(-args.Ly/2, args.Ly/2, args.Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
# Gitterindex nahe der Detektorebene
det_idx = int(np.argmin(np.abs(x - args.det_x)))
det_x_actual = x[det_idx]

# --- Doppelspalt-Wand bei x=0 ---
wall_half = 0.05
slit_y = 1.5
slit_half = 0.55
V0 = 120.0
if args.wall == 'sharp':
    V = np.zeros((args.Ny, args.Nx))
    mask = np.abs(x) < wall_half
    in_slit = (np.abs(y[:, None] - slit_y) < slit_half) | (np.abs(y[:, None] + slit_y) < slit_half)
    V = np.where(mask[None, :] & (~in_slit), V0, 0.0)
else:
    # Glatte Wand (Gauß in x) mit zwei glatten Spalten (Gauß in y)
    # Weiche Kanten: spektrale Ableitungen bleiben glatt (keine Gibbs-Artefakte)
    sigma_x = 0.05
    sigma_y = 0.30
    Vx = np.exp(-(x[None, :]**2)/(2*sigma_x**2))  # (1, Nx)
    g1 = np.exp(-((y[:, None] - slit_y)**2)/(2*sigma_y**2))  # (Ny, 1)
    g2 = np.exp(-((y[:, None] + slit_y)**2)/(2*sigma_y**2))  # (Ny, 1)
    Vy = 1.0 - g1 - g2 + g1*g2  # (Ny, 1)
    V = V0 * Vx * Vy  # (Ny, Nx) durch Broadcasting

# --- Anfangszustand: Gauß-Wellenpaket (Quantum Equilibrium) ---
x0, y0 = -12.0, 0.0
sig = 2.1
p0 = 5.0
psi = np.exp(-((x[None, :] - x0)**2 + (y[:, None] - y0)**2)/(2*sig**2)) * np.exp(1j*p0*x[None, :])
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)

# --- Frequenzen (Split-Step) ---
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
    """Spektrale Ableitungen (exakt fuer glatte periodische Funktionen)."""
    psik = np.fft.fft2(psi)
    px = np.fft.ifft2(psik * (1j*KX))
    py = np.fft.ifft2(psik * (1j*KY))
    return px, py

def central_derivs(psi):
    """Zentrale Differenzen (lokal, keine Gibbs-Artefakte an der Wand)."""
    px = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / (2*dx)
    py = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / (2*dy)
    return px, py

def interp_bilinear(field, xq, yq):
    """Periodische bilineare Interpolation (komplexe Felder)."""
    i = (xq + args.Lx/2) / dx
    j = (yq + args.Ly/2) / dy
    i0 = np.floor(i).astype(int) % args.Nx
    j0 = np.floor(j).astype(int) % args.Ny
    i1 = (i0 + 1) % args.Nx
    j1 = (j0 + 1) % args.Ny
    ai = i - np.floor(i)
    aj = j - np.floor(j)
    return (field[j0, i0]*(1-ai)*(1-aj) + field[j0, i1]*ai*(1-aj)
            + field[j1, i0]*(1-ai)*aj + field[j1, i1]*ai*aj)

def lagrange_weights(t, base, order=6):
    """Lagrange-Basisgewichte (order Punkte) an Position t mit Basisindex base."""
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
    """Periodische 6-Punkt-Lagrange-Interpolation (komplexe Felder)."""
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
    """Bohmsche Geschwindigkeit mit Regularisierung und Begrenzung."""
    eps = 1e-12
    vx = (hbar/m) * np.imag(px_i / (psi_i + eps))
    vy = (hbar/m) * np.imag(py_i / (psi_i + eps))
    speed = np.sqrt(vx**2 + vy**2)
    scale = np.where(speed > args.vmax, args.vmax/np.maximum(speed, 1e-12), 1.0)
    return vx*scale, vy*scale

# --- Teilchen: Anfangsverteilung |psi|^2 ---
prob = np.abs(psi)**2
prob_flat = prob.ravel() / prob.ravel().sum()
idx = rng.choice(args.Nx*args.Ny, size=args.Np, p=prob_flat)
iy, ix = np.unravel_index(idx, (args.Ny, args.Nx))
rx = rng.uniform(-0.5, 0.5, args.Np) * dx
ry = rng.uniform(-0.5, 0.5, args.Np) * dy
pos = np.stack([x[ix] + rx, y[iy] + ry], axis=1)

# --- Hauptschleife ---
nstep = int(args.T / args.dt)
hit_done = np.zeros(args.Np, dtype=bool)
hits_y = np.zeros(args.Np)  # y an Erstankunft (fuer hit_done Teilchen)
flux_profile = np.zeros(args.Ny)  # Summe j_x(y,t)*dt an der Detektorebene
traj_store = []  # Positionen alle 5 Schritte
max_speed_seen = 0.0
n_vmax_clipped = 0

for n in range(nstep):
    deriv_func = spectral_derivs if args.deriv == 'spectral' else central_derivs
    psi_old = psi
    px_old, py_old = deriv_func(psi_old)
    psi_new = split_step(psi_old, args.dt)
    px_new, py_new = deriv_func(psi_new)

    # RK4 mit linearer Interpolation von psi und Ableitungen zwischen t_n und t_{n+1}
    pos_old = pos.copy()
    interp_func = interp_bilinear if args.interp == 'bilinear' else interp_lagrange6
    def f(p, s):
        psi_s = (1-s)*psi_old + s*psi_new
        px_s = (1-s)*px_old + s*px_new
        py_s = (1-s)*py_old + s*py_new
        psi_i = interp_func(psi_s, p[:, 0], p[:, 1])
        px_i = interp_func(px_s, p[:, 0], p[:, 1])
        py_i = interp_func(py_s, p[:, 0], p[:, 1])
        vx, vy = bohm_velocity(psi_i, px_i, py_i)
        return np.stack([vx, vy], axis=1)

    k1 = f(pos, 0.0)
    k2 = f(pos + 0.5*args.dt*k1, 0.5)
    k3 = f(pos + 0.5*args.dt*k2, 0.5)
    k4 = f(pos + args.dt*k3, 1.0)
    pos_new = pos + (args.dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    # Erstankunft: Kreuzung von links nach rechts an der Detektorebene
    # (pos_old <= det_x_actual und pos_new > det_x_actual, erster Durchgang)
    cross = (~hit_done) & (pos_old[:, 0] <= det_x_actual) & (pos_new[:, 0] > det_x_actual)
    if np.any(cross):
        s_cross = (det_x_actual - pos_old[:, 0]) / (pos_new[:, 0] - pos_old[:, 0] + 1e-15)
        s_cross = np.clip(s_cross, 0.0, 1.0)
        y_cross = pos_old[:, 1] + s_cross * (pos_new[:, 1] - pos_old[:, 1])
        hits_y[cross] = y_cross[cross]
        hit_done[cross] = True

    # Wahrscheinlichkeitsfluss j_x an der Detektorebene (Gitterindex det_idx)
    # j_x = (hbar/m) * Im(conj(psi) * d_x psi)
    jx_det = (hbar/m) * np.imag(np.conj(psi[:, det_idx]) * px_old[:, det_idx])
    flux_profile += jx_det * args.dt

    # Diagnose: Geschwindigkeitsbegrenzung
    speed_cur = np.sqrt(k1[:, 0]**2 + k1[:, 1]**2)
    max_speed_seen = max(max_speed_seen, float(np.max(speed_cur)))
    n_vmax_clipped += int(np.sum(speed_cur > args.vmax))

    pos = pos_new
    psi = psi_new
    if n % 5 == 0:
        traj_store.append(pos.copy())

# --- Erstankunftsverteilung (nur hit_done Teilchen, kein Rueckfluss) ---
n_hits = int(np.sum(hit_done))
hits_y_used = hits_y[hit_done]

# --- Referenz: zeitintegrierter positiver Wahrscheinlichkeitsfluss ueber dieselben Bins ---
# Masse pro Bin = sum_t ( sum_{y in bin} j_x(y,t) * dy ) * dt
bin_edges = np.linspace(-args.Ly/2, args.Ly/2, args.n_bins + 1)
# Indizes der Gitterpunkte y in jedem Bin
bin_idx = np.digitize(y, bin_edges) - 1  # 0..n_bins-1, Randpunkte werden dem letzten/firsten Bin zugeordnet
bin_idx = np.clip(bin_idx, 0, args.n_bins - 1)
exp_mass = np.zeros(args.n_bins)
for b in range(args.n_bins):
    mask_b = bin_idx == b
    if np.any(mask_b):
        # summe j_x(y,t)*dt ueber Zeit und Gitterpunkte im Bin
        # flux_profile ist bereits sum_t j_x(y,t)*dt
        exp_mass[b] = np.sum(flux_profile[mask_b]) * dy
# Positiver Fluss (Teilchen, die von links nach rechts laufen)
exp_mass = np.maximum(exp_mass, 0.0)

# --- Histogramm der Erstankunfts-y-Werte ueber dieselben Bins ---
hist_counts, _ = np.histogram(hits_y_used, bins=bin_edges)
# --- Normierung beider Verteilungen auf Wahrscheinlichkeit 1 ---
obs_norm = hist_counts.astype(float) / np.sum(hist_counts) if np.sum(hist_counts) > 0 else np.zeros(args.n_bins)
exp_norm = exp_mass / np.sum(exp_mass) if np.sum(exp_mass) > 0 else np.zeros(args.n_bins)

# --- Relativer L2-Fehler zwischen beobachteter und erwarteter Verteilung ---
denom = np.sqrt(np.sum(exp_norm**2)) if np.sum(exp_norm**2) > 0 else 1e-15
arrival_rel_l2 = float(np.sqrt(np.sum((obs_norm - exp_norm)**2)) / denom)

# --- Absolute Trefferwahrscheinlichkeit (getrennte Pruefung) ---
# empirische Ankunftswahrscheinlichkeit vs. integrierte positive Flussmasse
hit_prob_emp = n_hits / args.Np
flux_total_mass = float(np.sum(exp_mass))  # unnormierte Masse (Integral j_x^+ dt dy)

# --- Randbehandlung: Pruefe, dass die Wellenfunktion am Rand klein ist ---
psi_edge_max = float(max(np.max(np.abs(psi[0, :])), np.max(np.abs(psi[-1, :])),
                          np.max(np.abs(psi[:, 0])), np.max(np.abs(psi[:, -1]))))
max_x_erreicht = float(np.max(pos[:, 0]))
# Plausibilitaetspruefung: Masse rechts der Detektorebene und links der Wand am Ende
# (dient als unabhaengige Pruefung, ob flux_total_mass die transmittierte Masse misst)
mask_right = x > det_x_actual
mask_left = x < 0.0
mass_right_end = float(np.sum(np.abs(psi[:, mask_right])**2) * dx * dy)
mass_left_end = float(np.sum(np.abs(psi[:, mask_left])**2) * dx * dy)
# Anteil der Teilchen, die noch aktiv sind (nicht am Detektor angekommen)
frac_active = float(np.mean(~hit_done))

# --- Diagnose: Verteilung und Fluss (fuer Plots und Text) ---
bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

# --- Plot: Detektor (Histogramm vs. Flussreferenz) ---
os.makedirs('plots', exist_ok=True)
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(bin_centers, obs_norm, width=(bin_edges[1]-bin_edges[0])*0.8, alpha=0.6,
       label='Bohm-Erstankunft (Histogramm)')
ax.plot(bin_centers, exp_norm, 'r-', lw=2, label='Fluss-Referenz (positive Masse)')
ax.set_xlabel('$y$ (Erstankunftsposition)')
ax.set_ylabel('normalisierte Wahrscheinlichkeit')
ax.set_title(f'Doppelspalt: Erstankunft vs. Fluss (rel. L2={arrival_rel_l2:.4f})')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig('plots/double_slit_v2_detector.png', dpi=150); plt.close()

# --- Plot: Trajektorien (alle, nicht nur Hits) ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.imshow(np.abs(psi.T)**2, extent=[-args.Lx/2, args.Lx/2, -args.Ly/2, args.Ly/2],
          origin='lower', cmap='Blues', alpha=0.8, aspect='auto')
for k in range(min(args.Np, 500)):
    xs = [p[k, 0] for p in traj_store]
    ys = [p[k, 1] for p in traj_store]
    ax.plot(xs, ys, 'r-', lw=0.4, alpha=0.5)
ax.axvline(0, color='k', lw=3)
ax.axvline(args.det_x, color='g', lw=2, ls='--')
ax.set_xlim(-args.Lx/2, args.Lx/2); ax.set_ylim(-args.Ly/2, args.Ly/2)
ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title('Bohmsche Trajektorien (korrigiert)')
plt.tight_layout(); plt.savefig('plots/double_slit_v2_trajectories.png', dpi=150); plt.close()

# --- Ergebnis-JSON im Pruefvertragsformat ---
result = {
    "observable": "detector_first_arrival_rel_l2_error",
    "reference": "zeitintegrierter positiver Wahrscheinlichkeitsfluss ueber die Detektorebene x=10 (Standard-Fluss-Vorhersage der Born-Regel)",
    "metrics": {
        "arrival_rel_l2_error": arrival_rel_l2,
        "arrival_grid_change": None,
        "arrival_monte_carlo_std": None
    },
    "diagnostics": {
        "n_particles": args.Np,
        "n_hits": n_hits,
        "hit_prob_empirical": hit_prob_emp,
        "flux_total_mass": flux_total_mass,
        "hit_prob_vs_flux": hit_prob_emp - flux_total_mass,
        "det_x_actual": float(det_x_actual),
        "det_idx": det_idx,
        "dt": args.dt,
        "T": args.T,
        "Nx": args.Nx,
        "Ny": args.Ny,
        "dx": float(dx),
        "dy": float(dy),
        "vmax": args.vmax,
        "max_speed_seen": max_speed_seen,
        "n_vmax_clipped": n_vmax_clipped,
        "psi_edge_max": psi_edge_max,
        "max_x_erreicht": max_x_erreicht,
        "mass_right_end": mass_right_end,
        "mass_left_end": mass_left_end,
        "frac_active": frac_active,
        "norm_end": float(np.sum(np.abs(psi)**2) * dx * dy),
        "n_bins": args.n_bins,
        "interp": args.interp,
        "deriv": args.deriv,
        "wall": args.wall,
        "seed": args.seed
    }
}
with open(args.out, 'w') as f:
    json.dump(result, f, indent=2)
print(json.dumps(result, indent=2))
