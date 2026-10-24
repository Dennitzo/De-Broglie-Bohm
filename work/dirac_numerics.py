"""Periodic central-difference Dirac Hamiltonian with exact CN link continuity.

Natural units hbar=c=1. FFT diagonalizes the same circulant spatial Hamiltonian;
it avoids constructing or factorizing a dense 2N by 2N matrix. No run on import.
"""
import numpy as np
from bohm_numerics import positive, step_count


def run(N=512, dt=0.005, m=1.0, L=40.0, T=5.0, x0=-10.0, p0=2.0, sig=1.2):
    for value, name in ((m,"m"), (L,"L"), (sig,"sig")):
        positive(value, name)
    if N < 8:
        raise ValueError("N must be >= 8")
    count = step_count(T, dt)
    dx = L/N
    x = -L/2 + (np.arange(N)+0.5)*dx
    k = 2*np.pi*np.fft.fftfreq(N, d=dx)
    momentum_lattice = np.sin(k*dx)/dx
    energy = np.sqrt(momentum_lattice**2+m**2)
    gaussian = np.exp(-(x-x0)**2/(2*sig**2)+1j*p0*x)
    seed_k = np.fft.fft(gaussian)
    # P_+(k) applied to (Gaussian,0): positive energy for THIS lattice H.
    psi_k = np.column_stack(((1+m/energy)*seed_k/2,
                            momentum_lattice*seed_k/(2*energy)))
    psi = np.fft.ifft(psi_k, axis=0)
    psi /= np.sqrt(np.sum(abs(psi)**2)*dx)
    psi_k = np.fft.fft(psi, axis=0)
    half = dt/2
    denominator = 1+(half*energy)**2
    scalar = (1-(half*energy)**2)/denominator
    matrix_factor = -2j*half/denominator
    max_link = max_continuum = max_charge = 0.0
    density_initial = np.sum(abs(psi)**2, axis=1)
    for _ in range(count):
        old = psi
        hpsi = np.column_stack((m*psi_k[:,0]+momentum_lattice*psi_k[:,1],
                               momentum_lattice*psi_k[:,0]-m*psi_k[:,1]))
        psi_k = scalar[:,None]*psi_k + matrix_factor[:,None]*hpsi
        psi = np.fft.ifft(psi_k, axis=0)
        mid = (old+psi)/2
        old_density = np.sum(abs(old)**2, axis=1)
        new_density = np.sum(abs(psi)**2, axis=1)
        alpha_mid = mid[:, ::-1]
        # Link i+1/2, and its divergence at node i, exactly match central H.
        link = np.real(np.sum(np.conj(mid)*np.roll(alpha_mid,-1,axis=0),axis=1))
        link_residual = (new_density-old_density)/dt+(link-np.roll(link,1))/dx
        onsite = np.real(np.sum(np.conj(mid)*alpha_mid,axis=1))
        continuum_residual = (new_density-old_density)/dt+(np.roll(onsite,-1)-np.roll(onsite,1))/(2*dx)
        max_link = max(max_link, float(np.max(abs(link_residual))))
        max_continuum = max(max_continuum, float(np.max(abs(continuum_residual))))
        max_charge = max(max_charge, abs(float(np.sum(new_density)*dx)-1))
    density_final = np.sum(abs(psi)**2,axis=1)
    return {"N": N, "dt": dt, "dx": dx, "T": count*dt,
            "max_discrete_link_residual": max_link,
            "max_continuum_onsite_residual": max_continuum,
            "max_charge_error": max_charge,
            "charge_initial": float(np.sum(density_initial)*dx),
            "charge_final": float(np.sum(density_final)*dx)}, {
                "x": x, "density_initial": density_initial, "density_final": density_final,
                "psi_final": psi}
