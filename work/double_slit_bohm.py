"""
Bohmsche Mechanik: Doppelspalt-Experiment (2D, Split-Step-Fourier).
Loest die Schroedinger-Gleichung fuer ein Wellenpaket an einem Doppelspalt,
berechnet Bohmsche Trajektorien (v = (hbar/m) Im(nabla psi/psi)) und
erzeugt drei PNG-Plots:
  1) double_slit_density.png       : |psi|^2 am Detektor (Profil + Histogramm)
  2) double_slit_trajectories.png  : Bohmsche Trajektorien auf |psi|^2
  3) double_slit_wave.png          : |psi|^2-Entwicklung (Zeitreihe)
Ausgabe: JSON mit Born-Regel-Vergleich (Chi^2) und Kennzahlen.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os

hbar = 1.0
m = 1.0

# --- Gebiet und Gitter ---
Lx, Ly = 30.0, 10.0
Nx, Ny = 512, 256
x = np.linspace(-Lx/2, Lx/2, Nx)
y = np.linspace(-Ly/2, Ly/2, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]

# --- Doppelspalt-Wand bei x=0, Spalten bei y=+-1.5 ---
wall_half = 0.05
slit_y = 1.5
slit_half = 0.55
V = np.zeros((Ny, Nx))
mask = (np.abs(x) < wall_half)
in_slit = (np.abs(y[:, None] - slit_y) < slit_half) | (np.abs(y[:, None] + slit_y) < slit_half)
V = np.where(mask[None, :] & (~in_slit), 120.0, 0.0)

# --- Anfangsbedingung: Gauß-Wellenpaket (Quantum-Equilibrium-Verteilung) ---
x0, y0 = -12.0, 0.0
sig = 2.1
p0 = 5.0
psi = np.exp(-((x[None, :]-x0)**2 + (y[:, None]-y0)**2)/(2*sig**2)) * np.exp(1j*p0*x[None, :])
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)

# --- Frequenzen (Split-Step) ---
kx = np.fft.fftfreq(Nx, d=dx) * 2*np.pi   # Achse 1 (x)
ky = np.fft.fftfreq(Ny, d=dy) * 2*np.pi   # Achse 0 (y)
KX = np.outer(np.ones(Ny), kx)
KY = np.outer(ky, np.ones(Nx))
K2 = KX**2 + KY**2

def split_step(psi, dt):
    psi = psi * np.exp(-1j*V*dt/2)
    psik = np.fft.fft2(psi)
    psik = psik * np.exp(-1j*K2*dt/(2*m))
    psi = np.fft.ifft2(psik)
    psi = psi * np.exp(-1j*V*dt/2)
    return psi

# --- Gitter-Ableitungen und Interpolation ---
def grid_derivs(psi):
    # zentrale Differenzen, periodische Randbedingung
    px = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / (2*dx)
    py = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / (2*dy)
    return px, py

def interp_field(field, xq, yq):
    # bilineare Interpolation mit periodischen Randbedingungen
    i = (xq + Lx/2) / dx
    j = (yq + Ly/2) / dy
    i0 = np.floor(i).astype(int) % Nx
    j0 = np.floor(j).astype(int) % Ny
    i1 = (i0 + 1) % Nx
    j1 = (j0 + 1) % Ny
    ai = i - np.floor(i)
    aj = j - np.floor(j)
    return (field[j0, i0]*(1-ai)*(1-aj) + field[j0, i1]*ai*(1-aj)
            + field[j1, i0]*(1-ai)*aj + field[j1, i1]*ai*aj)

def bohm_velocity(psi, px, py, xq, yq):
    psi_i = interp_field(psi, xq, yq)
    pxi = interp_field(px, xq, yq)
    pyi = interp_field(py, xq, yq)
    eps = 1e-12
    vx = (hbar/m) * np.imag(pxi / (psi_i + eps))
    vy = (hbar/m) * np.imag(pyi / (psi_i + eps))
    # Geschwindigkeitsbegrenzung gegen Divergenz an psi-Nullstellen
    vmax = 30.0
    speed = np.sqrt(vx**2 + vy**2)
    scale = np.where(speed > vmax, vmax / np.maximum(speed, 1e-12), 1.0)
    return vx * scale, vy * scale

# --- Parameter der Simulation ---
dt = 0.01
T_total = 16.0
nstep = int(T_total / dt)
N_particles = 10000

# Anfangsverteilung der Teilchen: |psi|^2 (Quantum Equilibrium)
prob = np.abs(psi)**2
prob_flat = prob.ravel()
prob_flat = prob_flat / prob_flat.sum()
idx = np.random.choice(Nx*Ny, size=N_particles, p=prob_flat)
iy, ix = np.unravel_index(idx, (Ny, Nx))
# kleine Zufallsverschiebung innerhalb der Zelle
rx = np.random.uniform(-0.5, 0.5, N_particles) * dx
ry = np.random.uniform(-0.5, 0.5, N_particles) * dy
pos_x = x[ix] + rx
pos_y = y[iy] + ry
positions = np.stack([pos_x, pos_y], axis=1)

# Detektor
det_x = 10.0
det_idx = np.argmin(np.abs(x - det_x))
traj_store = []   # Trajektorien (alle 5 Schritte)
wave_store = []   # |psi|^2 zu ausgewaehlten Zeiten
flux_store = []   # j_x am Detektor (alle 10 Schritte)
record_times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
record_idx = set(int(t/dt) for t in record_times)

def rk4_step(psi_old, px_old, py_old, psi_new, px_new, py_new, positions, dt):
    # lineare Interpolation von psi und Ableitungen zwischen t_n und t_{n+1}
    def f(pos, s):
        psi_s = (1-s)*psi_old + s*psi_new
        px_s = (1-s)*px_old + s*px_new
        py_s = (1-s)*py_old + s*py_new
        vx, vy = bohm_velocity(psi_s, px_s, py_s, pos[:,0], pos[:,1])
        return np.stack([vx, vy], axis=1)
    k1 = f(positions, 0.0)
    k2 = f(positions + 0.5*dt*k1, 0.5)
    k3 = f(positions + 0.5*dt*k2, 0.5)
    k4 = f(positions + dt*k3, 1.0)
    return positions + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

# --- Hauptschleife ---
for n in range(nstep):
    psi_old = psi
    px_old, py_old = grid_derivs(psi_old)
    psi_new = split_step(psi_old, dt)
    px_new, py_new = grid_derivs(psi_new)
    positions = rk4_step(psi_old, px_old, py_old, psi_new, px_new, py_new, positions, dt)
    psi = psi_new
    if n % 5 == 0:
        traj_store.append(positions.copy())
    if n in record_idx:
        wave_store.append((n*dt, np.abs(psi)**2))
    if n % 10 == 0:
        px_cur = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1))/(2*dx)
        jx_det = (hbar/m) * np.imag(px_cur[:, det_idx] / (psi[:, det_idx] + 1e-12)) * np.abs(psi[:, det_idx])**2
        flux_store.append((n*dt, jx_det))

# --- Auftreffpunkte am Detektor (erster Durchgang pro Teilchen) ---
hits_y = []
hit_done = np.zeros(N_particles, dtype=bool)
for n, positions in enumerate(traj_store):
    pass_det = np.where((positions[:,0] >= det_x - 0.1) & (~hit_done))[0]
    for k in pass_det:
        hits_y.append(positions[k,1])
        hit_done[k] = True
    if np.all(hit_done):
        break

# --- Born-Regel-Vergleich: integrierter Fluss am Detektor ---
flux_profile = np.zeros(Ny)
for t, jx in flux_store:
    flux_profile += jx * dt
# Nur positiver Fluss (Teilchen, die den Detektor von links treffen)
flux_profile = np.maximum(flux_profile, 0.0)
flux_profile = flux_profile / np.sum(flux_profile)

hist_counts, bin_edges = np.histogram(hits_y, bins=30)
bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
bin_width = bin_edges[1] - bin_edges[0]
# Born-Dichte am Detektor
prof_dens = np.interp(bin_centers, y, flux_profile)
prof_dens = prof_dens / np.sum(prof_dens)
# Pearson-Korrelation als robuste Formpruefung
hist_norm = hist_counts / np.sum(hist_counts)
rho = np.corrcoef(hist_norm, prof_dens)[0, 1]
# Chi^2 mit gepoolten Bins (erwartete Zaehlungen >= 5)
n_hits_total = len(hits_y)
expected = prof_dens * n_hits_total
chi2 = 0.0
obs_pooled = 0.0
exp_pooled = 0.0
for o, e in zip(hist_counts.astype(float), expected):
    obs_pooled += o
    exp_pooled += e
    if exp_pooled >= 5.0:
        if exp_pooled > 0:
            chi2 += (obs_pooled - exp_pooled)**2 / exp_pooled
        obs_pooled = 0.0
        exp_pooled = 0.0

# --- Plots ---
os.makedirs("plots", exist_ok=True)

# Plot 1: Detektor
fig, ax = plt.subplots(figsize=(8,5))
ax.bar(bin_centers, hist_counts/np.sum(hist_counts), width=bin_edges[1]-bin_edges[0], alpha=0.6,
       label='Bohm-Trajektorien (Histogramm)')
ax.plot(y, flux_profile/np.sum(flux_profile), 'r-', lw=2, label='$|\\psi|^2$ (Born)')
ax.set_xlabel('$y$ (Position am Detektor)'); ax.set_ylabel('relative Dichte')
ax.set_title('Doppelspalt: Bohmsche Auftreffpunkte vs. Born-Regel')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.savefig('plots/double_slit_detector.png', dpi=150); plt.close()

# Plot 2: Trajektorien auf |psi|^2
fig, ax = plt.subplots(figsize=(9,5))
ax.imshow(np.abs(psi.T)**2, extent=[-Lx/2, Lx/2, -Ly/2, Ly/2], origin='lower',
          cmap='Blues', alpha=0.8, aspect='auto')
for k in range(N_particles):
    xs = [p[k,0] for p in traj_store]
    ys = [p[k,1] for p in traj_store]
    ax.plot(xs, ys, 'r-', lw=0.4, alpha=0.5)
ax.axvline(0, color='k', lw=3)
ax.set_xlim(-Lx/2, Lx/2); ax.set_ylim(-Ly/2, Ly/2)
ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title('Bohmsche Trajektorien am Doppelspalt')
plt.tight_layout(); plt.savefig('plots/double_slit_trajectories.png', dpi=150); plt.close()

# Plot 3: Wellenentwicklung
fig, axes = plt.subplots(2, 3, figsize=(14,7))
for ax, (t, d) in zip(axes.ravel(), wave_store):
    ax.imshow(d.T, extent=[-Lx/2, Lx/2, -Ly/2, Ly/2], origin='lower',
              cmap='viridis', aspect='auto')
    ax.set_title(f'$t={t:.1f}$')
    ax.set_xlim(-Lx/2, Lx/2); ax.set_ylim(-Ly/2, Ly/2)
plt.tight_layout(); plt.savefig('plots/double_slit_wave.png', dpi=150); plt.close()

# Diagnose: Positionen zu ausgewählten Zeiten (alle 5 Schritte gespeichert)
diag = {}
for t, n in [(1.0, 100), (2.0, 200), (3.0, 300), (4.0, 400), (6.0, 600), (8.0, 800)]:
    idx = n // 5
    if idx < len(traj_store):
        pos = traj_store[idx]
        diag[t] = {"mean_x": float(np.mean(pos[:,0])),
                   "mean_y": float(np.mean(pos[:,1])),
                   "p95_x": float(np.percentile(pos[:,0], 95))}
max_x = float(np.max(positions[:, 0]))
result = {
    "diag_positions": diag, 
    "n_particles": N_particles,
    "n_hits": len(hits_y),
    "pearson_born": float(rho),
    "chi2_born": float(chi2),
    "norm_end": float(np.sum(np.abs(psi)**2)*dx*dy),
    "detector_density_max": float(np.max(flux_profile)),
    "simulation_time": T_total,
    "max_x_erreicht": max_x
}
print(json.dumps(result, indent=2))
