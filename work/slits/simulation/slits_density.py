"""Controlled paraxial diagnostic: i*u_x=-u_yy/(2*k), dy/dx=j/rho.

Fixed-x transverse densities are not physical-time 2D densities or detector
first arrivals. The independent reference is the closed Gaussian solution.
This corrected version does not certify the original experiment requirements.
"""
import numpy as np
import json, os
from paraxial_reference import (MODEL_VERSION, analytic_field, bin_masses, initial_particles,
    spectral_guidance, rk4_step, validate_parameters, source_hash)

def gauss(y, mu, sigma):
    return np.exp(-0.5 * ((y - mu) / sigma) ** 2)

def build_initial_field(Ny, Ymax, k, d, a):
    y = np.linspace(-Ymax, Ymax, Ny, endpoint=False)
    dy = 2.0 * Ymax / Ny
    return y, dy, analytic_field(y, 0.0, k, d, a)


def propagate_free(u0, Nx, L, dy, m):
    Ny = len(u0)
    dx = L / (Nx - 1)
    k_y = 2.0 * np.pi * np.fft.fftfreq(Ny, d=dy)
    phase_step = np.exp(-0.5j * k_y ** 2 * dx / m)
    u = np.zeros((Nx, Ny), dtype=complex)
    u[0] = u0
    for n in range(Nx - 1):
        u[n + 1] = np.fft.ifft(phase_step * np.fft.fft(u[n], axis=0))
    return u, dx

def exact_reference(y, L, k, d, a):
    return analytic_field(y, L, k, d, a)


def velocity_fields(u, dx, dy, eps, k):
    # Parameter speed along x is 1. No carrier-phase/density division twice.
    if eps != 0:
        raise ValueError('Density regularisation changes the model; eps must be 0.')
    ky = 2*np.pi*np.fft.fftfreq(u.shape[-1], d=dy)
    du = np.fft.ifft(1j*ky*np.fft.fft(u, axis=-1), axis=-1)
    rho = np.abs(u)**2
    vy = np.divide(np.imag(np.conj(u)*du), k*rho,
                   out=np.full(rho.shape, np.nan), where=rho>np.finfo(float).tiny)
    return np.ones(rho.shape), vy


def rel_l2(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    denom = np.sqrt(np.sum(b ** 2))
    if denom == 0:
        return float('nan')
    return float(np.sqrt(np.sum((a - b) ** 2)) / denom)

def fringe_spacing(y, I):
    """Observed central peak spacing, never a proof of Fraunhofer validity."""
    I = np.asarray(I)
    indices = np.where((I[1:-1] > I[:-2]) & (I[1:-1] >= I[2:]) &
                       (I[1:-1] > .05*np.max(I)))[0] + 1
    if len(indices) < 2:
        return float('nan')
    # Restrict to central lobes; envelope/DC spectral power is not fringes.
    central = indices[np.argsort(np.abs(y[indices]))[:5]]
    return float(np.median(np.diff(np.sort(y[central]))))


def run(params):
    validate_parameters(params)
    Ny = params['Ny']; Ymax = params['Ymax']; k = params['k']; d = params['d']; a = params['a']
    L = params['L']; Nx = params['Nx']; Np = params['Np']; eps = params['eps']
    x_start = params.get('x_start', 0.0); seed = params['seed']
    if x_start != 0.0:
        raise ValueError('This pilot fixes x_start=0.0; use a separate fixed-plane study for nonzero start.')
    check_times = params.get('check_times', [0.25, 0.5, 0.75])  # Bruchteile von L
    lam = 2.0 * np.pi / k
    m = params.get('m', k)

    y, dy, u0 = build_initial_field(Ny, Ymax, k, d, a)
    u, dx = propagate_free(u0, Nx, L, dy, m)

    # --- 1. Propagation vs. exakte Referenz ---
    uL_ref = exact_reference(y, L, k, d, a)
    I_num = np.abs(u[-1, :]) ** 2
    I_ref = np.abs(uL_ref) ** 2
    I_num = I_num / (np.sum(I_num) * dy)
    I_ref = I_ref / (np.sum(I_ref) * dy)
    prop_rel_l2 = rel_l2(I_num, I_ref)

    # --- Geschwindigkeitsfeld ---
    vx, vy = velocity_fields(u, dx, dy, eps, k)

    # --- 2. Teilchendichte zu festen Zeiten vs. |psi|^2 ---
    y_i = initial_particles(Np, seed, x_start, Ymax, k, d, a)
    f = spectral_guidance(y, u0, dy, k)
    # Fixed physical bins independent of solver grid; compare integrated masses.
    bin_edges = np.linspace(-Ymax, Ymax, int(params.get('bins', 128)) + 1)

    n_clipped = 0
    time_frac_to_idx = {}
    probability_checks = {}
    density_cells = {}
    target_positions = [float(tf)*L for tf in check_times]
    if any(x <= x_start or x > L for x in target_positions):
        raise ValueError('Check planes must lie after fixed x_start and before L.')
    positions = np.unique(np.r_[np.linspace(x_start, L, Nx), target_positions])
    for left, right in zip(positions[:-1], positions[1:]):
        y_i = rk4_step(f, float(left), y_i, float(right-left))
        if np.any(y_i < y[0]) or np.any(y_i > y[-1]):
            raise ValueError('Trajectory leaves domain; enlarge domain. No clipping.')
        if right in target_positions:
            H_t = np.histogram(y_i, bins=bin_edges)[0] / Np
            I_t = bin_masses(bin_edges, right, k, d, a)
            if abs(I_t.sum()-1.0)>1e-8:
                raise ValueError('Reference probability outside grid exceeds 1e-8; enlarge domain.')
            probability_checks[str(right)] = {
                'observed':float(H_t.sum()), 'reference':float(I_t.sum()),
                'iid_expected_relative_rms':float(np.sqrt((1-np.sum(I_t**2))/(Np*np.sum(I_t**2))))}
            density_cells[str(right)] = {'bin_edges': bin_edges.tolist(), 'observed_mass': H_t.tolist(), 'reference_mass': I_t.tolist()}
            time_frac_to_idx[right] = (H_t, I_t)

    density_errors = {}
    for tf in check_times:
        H_t, I_t = time_frac_to_idx[float(tf)*L]
        density_errors['x_frac_%.2f' % tf] = rel_l2(H_t, I_t)

    max_density_error = max(v for v in density_errors.values() if v is not None) if any(v is not None for v in density_errors.values()) else float('nan')

    # --- 3. Streifenabstand ---
    s_num = fringe_spacing(y, I_num)
    s_ref = fringe_spacing(y, I_ref)
    s_ana = lam * L / d

    # --- Diagnose ---
    issues = []
    if prop_rel_l2 > 0.03:
        issues.append('propagation_vs_reference_exceeds_0.03')
    if max_density_error > 0.05:
        issues.append('density_error_exceeds_0.05: %.4f' % max_density_error)
    if n_clipped > 0:
        issues.append('trajectory_slope_clipped=%d' % n_clipped)
    max_v = float(np.nanmax(np.sqrt(vx ** 2 + vy ** 2)))
    fraunhofer_valid = bool(L > 10.0 * max(a ** 2, d ** 2) / lam)

    metrics = {
        'prop_rel_l2': prop_rel_l2,
        'density_errors': density_errors,
        'density_cells': density_cells,
        'probability_checks': probability_checks,
        'max_density_error': max_density_error,
        'fringe_spacing_num': s_num,
        'fringe_spacing_ref': s_ref,
        'fringe_spacing_analytic': s_ana,
    }
    diagnostics = {
        'n_clipped': n_clipped,
        'model_version': MODEL_VERSION,
        'acceptance': 'not_established',
        'physical_time_density': False,
        'first_arrival_validation': False,
        'convergence_established': False,
        'independent_seeds': False,
        'max_velocity': max_v,
        'fraunhofer_valid': fraunhofer_valid,
        'norm_I_num': float(np.sum(I_num) * dy),
    }
    return metrics, diagnostics, issues, (y, I_num, I_ref), {
        'k': k, 'lambda': lam, 'd': d, 'a': a, 'L': L, 'Ny': Ny, 'Ymax': Ymax,
        'dy': dy, 'Nx': Nx, 'dx': dx, 'sigma_slit': a / 2.0, 'Np': Np, 'eps': eps,
        'x_start': x_start, 'seed': seed,
        'm': m, 'check_times': check_times, 'bins': len(bin_edges)-1,
    }

def main():
    tag = os.environ.get('SLITS_TAG', 'density')
    params = {
        'Ny': 1024, 'Ymax': 30.0, 'k': 2.0 * np.pi, 'd': 4.0, 'a': 1.0,
        'L': 20.0, 'Nx': 500, 'Np': 10000, 'eps': 0.0, 'x_start': 0.0,
        'seed': 42,
        'check_times': [0.25, 0.5, 0.75],
    }
    from paraxial_reference import output_paths
    json_path, npz_path = output_paths(tag, __file__)
    metrics, diagnostics, issues, (y, I_num, I_ref), p = run(params)
    here = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(here, 'corrected_v2_result_%s.json' % tag)
    npz_path = os.path.join(here, 'corrected_v2_arrays_%s.npz' % tag)
    payload = {
        'schema': MODEL_VERSION,
        'sourceSha256': source_hash(),
        'observable': 'paraxial fixed-x transverse density; not physical time or first arrival',
        'reference': 'Closed Gaussian superposition on the infinite line (independent of FFT)',
        'params': p,
        'metrics': metrics,
        'diagnostics': diagnostics,
        'issues': issues,
    }
    if os.path.exists(json_path) or os.path.exists(npz_path):
        raise FileExistsError('Use a new SLITS_TAG; do not overwrite measured results.')
    with open(json_path, 'w') as fp:
        json.dump(payload, fp, indent=2)
    np.savez_compressed(npz_path, y=y, I_num=I_num, I_ref=I_ref)
    import hashlib
    payload['arrays'] = os.path.basename(npz_path)
    payload['arraysSha256'] = hashlib.sha256(open(npz_path, 'rb').read()).hexdigest()
    with open(json_path, 'w', encoding='utf-8') as fp:
        json.dump(payload, fp, indent=2)
    m = metrics
    print('prop_rel_l2        = %.3e' % m['prop_rel_l2'])
    for k2, v in m['density_errors'].items():
        print('%-20s = %s' % (k2, '%.4f' % v if v is not None else 'None'))
    print('max_density_error  = %.4f' % m['max_density_error'])
    print('fringe num         = %.4f' % m['fringe_spacing_num'])
    print('fringe ref         = %.4f' % m['fringe_spacing_ref'])
    print('fringe analytic    = %.4f' % m['fringe_spacing_analytic'])
    print('n_clipped          = %d' % diagnostics['n_clipped'])
    print('max_velocity       = %.3f' % diagnostics['max_velocity'])
    print('fraunhofer_valid   = %s' % diagnostics['fraunhofer_valid'])
    print('issues             = %s' % issues)
    print('saved              = %s' % json_path)

if __name__ == '__main__':
    main()
