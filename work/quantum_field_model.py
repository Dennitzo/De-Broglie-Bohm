"""Finite-cutoff real scalar Bohmian field, with reproducible physical 3D fields.

python work/quantum_field_model.py --tag 20261017

This is a free 3+1 dimensional scalar field on a spatial torus, with an explicit
foliation and Fourier cutoff. Gaussian modal coordinates define the actual field;
the wavefunctional is a different object on the configuration space of all modes.
Only new quantum_field_* outputs are written, never replacing historical runs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import scipy
from scipy.integrate import quad

L = 2 * np.pi
VOLUME = L**3
MASS = 1.0
ELL = .55
SQUEEZE = .65
SEED = 20261017
FINAL_TIME = 2.4
GRID_N = 48
PLOT_TIMES = (0.0, 1.2)
ISO_LEVELS = (-.075, .075)
THRESHOLDS = {
    "orthonormality": 2e-12,
    "fft_energy_relative": 2e-12,
    "analytic_pde_relative": 2e-11,
    "finite_difference_pde_relative": 2e-5,
    "finite_difference_continuity_relative": 2e-5,
    "finite_difference_force_absolute": 2e-5,
    "finite_difference_guidance_absolute": 2e-8,
    "energy_relative": 2e-12,
    "flow_log_probability_error": 2e-10,
    "rk4_fine_absolute": 2e-6,
    "phase_quadrature_error": 2e-10,
    "phase_branch_jump": 1e-5,
    "vacuum_static_error": 2e-12,
    "smoothed_K6_K8_l2_difference": 2e-3,
    "smoothed_K8_K16_variance_tail": 2e-10,
    "isosurface_sphere_radius_error": .01,
}


def make_modes(cutoff):
    """One canonical half of nonzero lattice, two real modes per vector."""
    half = []
    for nx in range(-cutoff, cutoff + 1):
        for ny in range(-cutoff, cutoff + 1):
            for nz in range(-cutoff, cutoff + 1):
                n = (nx, ny, nz)
                nonzero = next((j for j in n if j), 0)
                if nonzero > 0:
                    half.append(n)
    half.sort(key=lambda n: (sum(j*j for j in n), n))
    vectors, kinds = [(0, 0, 0)], ["zero"]
    for n in half:
        vectors.extend([n, n])
        kinds.extend(["cos", "sin"])
    return np.array(vectors, dtype=int), np.array(kinds)


class GaussianField:
    def __init__(self, cutoff=8, vacuum=False):
        self.cutoff = cutoff
        self.n, self.kind = make_modes(cutoff)
        self.n2 = np.sum(self.n**2, axis=1)
        self.omega = np.sqrt(self.n2 + MASS**2)
        self.r = np.where(self.n2 <= 2, SQUEEZE, 0.0)
        self.mu0 = np.zeros(len(self.n))
        self.pi0 = np.zeros(len(self.n))
        labels = {(tuple(n), kind): k for k, (n, kind) in enumerate(zip(self.n, self.kind))}
        for label, value in [(((0, 0, 0), "zero"), 1.0), (((1, 0, 0), "cos"), 2.0),
                             (((0, 1, 0), "sin"), 1.5), (((0, 0, 1), "cos"), -1.0)]:
            self.mu0[labels[label]] = value
        for label, value in [(((1, 0, 0), "sin"), 1.25), (((0, 1, 0), "cos"), -.75)]:
            self.pi0[labels[label]] = value
        if vacuum:
            self.r[:] = 0
            self.mu0[:] = 0
            self.pi0[:] = 0
        self.sigma0 = np.exp(self.r) / np.sqrt(2 * self.omega)
        self.smoothing = np.exp(-ELL**2 * self.n2 / 2)

    def parameters(self, t):
        w = self.omega
        c, s = np.cos(w*t), np.sin(w*t)
        mu = self.mu0*c + self.pi0*s/w
        pi = -self.mu0*w*s + self.pi0*c
        aa = self.sigma0**2
        bb = 1 / (4*w*w*aa)
        variance = aa*c*c + bb*s*s
        variance_dot = 2*w*(bb-aa)*s*c
        # Exactly static vacuum coordinates, without trigonometric roundoff.
        variance = np.where(self.r == 0, aa, variance)
        variance_dot = np.where(self.r == 0, 0.0, variance_dot)
        sigma = np.sqrt(variance)
        sigma_dot = variance_dot / (2*sigma)
        b = sigma_dot / sigma
        sigma_ddot = 1/(4*sigma**3) - w*w*sigma
        bdot = sigma_ddot/sigma - b*b
        return mu, pi, sigma, sigma_dot, b, bdot

    def theta(self, t):
        angle = self.omega*t
        kappa = np.exp(-2*self.r)
        c, s = np.cos(angle), np.sin(angle)
        # The atan2 denominator is strictly positive. This pointwise formula
        # supplies the continuous argument at arbitrary time order and at pi/w.
        chi = angle + np.arctan2((kappa-1)*s*c, c*c+kappa*s*s)
        mu, pi, *_ = self.parameters(t)
        return .5*(mu*pi - self.mu0*self.pi0) - .5*chi

    def exact(self, t, xi):
        mu, _, sigma, *_ = self.parameters(t)
        return mu + sigma*xi

    def velocity(self, t, q):
        mu, pi, _, _, b, _ = self.parameters(t)
        return pi + b*(q-mu)

    def energy(self, t):
        mu, pi, sigma, sigma_dot, *_ = self.parameters(t)
        coherent = .5*(pi*pi + self.omega**2*mu*mu)
        fluctuation = .5*sigma_dot**2 + .5*self.omega**2*sigma**2 + 1/(8*sigma**2)
        return coherent + fluctuation

    def logpsi(self, t, q, ids):
        mu, pi, sigma, _, b, _ = (p[ids] for p in self.parameters(t))
        y = q-mu
        logamp = -.25*np.log(2*np.pi*sigma**2).sum() - np.sum(y*y/(4*sigma**2))
        phase = np.sum(pi*y + .5*b*y*y + self.theta(t)[ids])
        return logamp + 1j*phase


def rk4(model, initial, dt):
    steps = int(round(FINAL_TIME/dt))
    assert abs(steps*dt-FINAL_TIME) < 1e-13
    out = np.empty((steps+1,) + initial.shape)
    out[0] = initial
    for i in range(steps):
        q, t = out[i], i*dt
        a = model.velocity(t, q)
        b = model.velocity(t+dt/2, q+dt*a/2)
        c = model.velocity(t+dt/2, q+dt*b/2)
        d = model.velocity(t+dt, q+dt*c)
        out[i+1] = q + dt*(a+2*b+2*c+d)/6
    return np.arange(steps+1)*dt, out


def fft_field(model, coefficients, grid_n=GRID_N, smoothed=True, cutoff=None):
    spectrum = np.zeros((grid_n,)*3, dtype=complex)
    coefficients = coefficients * (model.smoothing if smoothed else 1)
    spectrum[0,0,0] = coefficients[0]/np.sqrt(VOLUME)
    for i in range(1, len(model.n), 2):
        n = model.n[i]
        if cutoff is not None and np.max(abs(n)) > cutoff:
            continue
        index, negative = tuple(n % grid_n), tuple((-n) % grid_n)
        value = (coefficients[i]-1j*coefficients[i+1])/np.sqrt(2*VOLUME)
        spectrum[index] = value
        spectrum[negative] = value.conjugate()
    field_complex = np.fft.ifftn(spectrum)*grid_n**3
    return field_complex.real, spectrum, float(np.max(abs(field_complex.imag)))


def basis_values(model, points, ids):
    result = np.empty((len(points), len(ids)))
    for j, index in enumerate(ids):
        if model.kind[index] == "zero":
            result[:, j] = 1/np.sqrt(VOLUME)
        else:
            angle = points @ model.n[index]
            fun = np.cos if model.kind[index] == "cos" else np.sin
            result[:, j] = np.sqrt(2/VOLUME)*fun(angle)
    return result


def analytic_checks(model, xi, arrays):
    errors = {"schrodinger_ratio_max": 0.0, "continuity_log_ratio_max": 0.0,
              "guidance_force_max": 0.0, "flow_log_probability_max": 0.0}
    energies, excitation, logdet = [], [], []
    ids = np.where(model.n2 <= 2)[0]
    for t in np.linspace(0, FINAL_TIME, 25):
        mu, pi, sigma, sd, b, bd = model.parameters(t)
        q = model.exact(t, xi)
        y = q-mu
        velocity = pi+b*y
        quantum = 1/(4*sigma*sigma)-y*y/(8*sigma**4)
        potential = .5*model.omega**2*q*q
        theta_dot = .5*(pi*pi-model.omega**2*mu*mu)-1/(4*sigma*sigma)
        phase_dt = -model.omega**2*mu*y-pi*pi+.5*bd*y*y-b*pi*y+theta_dot
        logRdt = -.5*b+y*pi/(2*sigma*sigma)+b*y*y/(2*sigma*sigma)
        gradR = -y/(2*sigma*sigma)
        real_residual = phase_dt+.5*velocity**2+potential+quantum
        imag_residual = logRdt+gradR*velocity+.5*b
        errors["schrodinger_ratio_max"] = max(errors["schrodinger_ratio_max"], float(np.max(abs(real_residual+1j*imag_residual))))
        logrho_dt = -b+y*pi/(sigma*sigma)+b*y*y/(sigma*sigma)
        continuity = logrho_dt+b-y*velocity/(sigma*sigma)
        errors["continuity_log_ratio_max"] = max(errors["continuity_log_ratio_max"], float(np.max(abs(continuity))))
        acceleration = -model.omega**2*mu+(1/(4*sigma**3)-model.omega**2*sigma)*xi
        force = -model.omega**2*q+y/(4*sigma**4)
        errors["guidance_force_max"] = max(errors["guidance_force_max"], float(np.max(abs(acceleration-force))))
        logJ = np.log(sigma/model.sigma0).sum()
        logp = -.5*np.log(2*np.pi*sigma*sigma).sum()-.5*np.sum(((q-mu)/sigma)**2)
        q0 = model.exact(0, xi)
        logp0 = -.5*np.log(2*np.pi*model.sigma0**2).sum()-.5*np.sum(((q0-model.mu0)/model.sigma0)**2)
        errors["flow_log_probability_max"] = max(errors["flow_log_probability_max"], abs(float(logp+logJ-logp0)))
        total_energy = model.energy(t)
        energies.append(float(total_energy.sum()))
        excitation.append(float((total_energy-.5*model.omega).sum()))
        logdet.append(float(logJ))
    arrays.update(energy_times=np.linspace(0, FINAL_TIME, 25), energy_total=np.array(energies),
                  energy_excitation=np.array(excitation), flow_logdet=np.array(logdet))
    errors["energy_total_relative_drift"] = float(np.ptp(energies)/energies[0])
    errors["energy_excitation_relative_drift"] = float(np.ptp(excitation)/excitation[0])
    errors["energy_total"] = energies[0]
    errors["energy_vacuum_subtracted"] = excitation[0]
    errors["squeezed_modes"] = len(ids)
    return errors


def finite_difference_checks(model, xi):
    ids = np.arange(7)
    t = .37
    mu, pi, sigma, _, b, _ = model.parameters(t)
    q = model.exact(t, xi)[ids]
    psi = np.exp(model.logpsi(t, q, ids))
    rows = []
    for h in (.00025, .000125, .0000625):
        dpsi_t = (np.exp(model.logpsi(t+h,q,ids))-np.exp(model.logpsi(t-h,q,ids)))/(2*h)
        laplacian, divergence = 0j, 0.0
        guidance_error, qforce_error = 0.0, 0.0
        rho = abs(psi)**2
        def quantum_at(qvalue):
            y = qvalue-mu[ids]
            return np.sum(1/(4*sigma[ids]**2)-y*y/(8*sigma[ids]**4))
        for j in range(len(ids)):
            step = np.zeros(len(ids)); step[j] = h
            plus, minus = model.logpsi(t,q+step,ids), model.logpsi(t,q-step,ids)
            laplacian += (np.exp(plus)-2*psi+np.exp(minus))/h**2
            vp = pi[ids[j]]+b[ids[j]]*(q[j]+h-mu[ids[j]])
            vm = pi[ids[j]]+b[ids[j]]*(q[j]-h-mu[ids[j]])
            divergence += (np.exp(2*plus.real)*vp-np.exp(2*minus.real)*vm)/(2*h)
            guidance_error = max(guidance_error, abs((plus.imag-minus.imag)/(2*h)-(vp+vm)/2))
            gradient_q = (quantum_at(q+step)-quantum_at(q-step))/(2*h)
            exact_gradient = -(q[j]-mu[ids[j]])/(4*sigma[ids[j]]**4)
            qforce_error = max(qforce_error, abs(gradient_q-exact_gradient))
        hamiltonian = -.5*laplacian+.5*np.sum(model.omega[ids]**2*q*q)*psi
        rho_dt = (np.exp(2*model.logpsi(t+h,q,ids).real)-np.exp(2*model.logpsi(t-h,q,ids).real))/(2*h)
        qplus, qzero, qminus = (model.exact(t+shift,xi)[ids] for shift in (h,0,-h))
        acceleration_fd = (qplus-2*qzero+qminus)/h**2
        force = -model.omega[ids]**2*q+(q-mu[ids])/(4*sigma[ids]**4)
        rows.append({"h":h,"dimensions":len(ids),
                     "schrodinger_relative_residual":float(abs(1j*dpsi_t-hamiltonian)/abs(psi)),
                     "continuity_relative_residual":float(abs(rho_dt+divergence)/rho),
                     "phase_guidance_absolute_error":float(guidance_error),
                     "quantum_potential_gradient_absolute_error":float(qforce_error),
                     "trajectory_force_absolute_error":float(np.max(abs(acceleration_fd-force)))})
    phase_errors=[]
    for index in ids:
        for tau in (.37, 1.2, 2.4, 4.5, np.pi/model.omega[index], 2*np.pi/model.omega[index]):
            def integrand(s):
                mu, pi, sigma, *_ = model.parameters(s)
                return .5*(pi[index]**2-model.omega[index]**2*mu[index]**2)-1/(4*sigma[index]**2)
            integrated=quad(integrand,0,tau,epsabs=2e-12,epsrel=2e-12,limit=150)[0]
            phase_errors.append(abs(integrated-model.theta(tau)[index]))
    branch_errors=[]
    for index in ids:
        tau=np.pi/model.omega[index]
        branch_errors.append(abs(model.theta(tau+1e-7)[index]-model.theta(tau-1e-7)[index]))
    return {"finite_differences":rows,"phase_vs_adaptive_integral_max":float(max(phase_errors)),
            "continuous_phase_branch_max_change":float(max(branch_errors))}


def variance_cutoffs():
    rows=[]
    for cutoff in (2,4,6,8,12,16):
        n=np.indices((2*cutoff+1,)*3).reshape(3,-1).T-cutoff
        n2=np.sum(n*n,axis=1)
        omega=np.sqrt(n2+MASS**2)
        raw=float(np.sum(1/(2*omega))/VOLUME)
        smooth=float(np.sum(np.exp(-ELL**2*n2)/(2*omega))/VOLUME)
        rows.append({"K":cutoff,"real_modes":int(len(n)),"vacuum_point_variance_raw":raw,
                     "vacuum_point_variance_smoothed":smooth})
    return rows


def marching_tetrahedra(field, level):
    """Actual xyz isosurface triangles; no image-generation or height-map proxy."""
    g=field.shape[0]
    extended=np.pad(field,((0,1),)*3,mode="wrap")
    offsets=np.array([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)])
    corners=np.stack([extended[a:a+g,b:b+g,c:c+g].ravel() for a,b,c in offsets],axis=1)
    crossing=(corners.min(axis=1)<level)&(corners.max(axis=1)>level)
    base=np.indices((g,g,g)).reshape(3,-1).T[crossing]
    corners=corners[crossing]
    vertices=(base[:,None,:]+offsets[None,:,:])*(L/g)
    tetrahedra=[(0,5,1,6),(0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6)]
    triangles=[]
    for tetra in tetrahedra:
        values=corners[:,tetra]
        coords=vertices[:,tetra]
        code=np.sum((values>level)*np.array([1,2,4,8]),axis=1)
        for case in range(1,15):
            selected=np.where(code==case)[0]
            if not len(selected): continue
            inside=[i for i in range(4) if case&(1<<i)]
            outside=[i for i in range(4) if not case&(1<<i)]
            def edge(a,b):
                fraction=(level-values[selected,a])/(values[selected,b]-values[selected,a])
                return coords[selected,a]+fraction[:,None]*(coords[selected,b]-coords[selected,a])
            if len(inside)==1 or len(inside)==3:
                one,many=(inside,outside) if len(inside)==1 else (outside,inside)
                triangles.append(np.stack([edge(one[0],j) for j in many],axis=1))
            else:
                a,b=inside;c,d=outside
                ac,ad,bc,bd=edge(a,c),edge(a,d),edge(b,c),edge(b,d)
                triangles.extend([np.stack([ac,ad,bc],axis=1),np.stack([ad,bd,bc],axis=1)])
    return np.concatenate(triangles) if triangles else np.empty((0,3,3))


def write_figures(root, stem, fields, arrays, metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    fig=plt.figure(figsize=(13,8),layout="constrained")
    maximum=max(float(np.max(abs(v))) for pair in fields for v in pair)
    color_limit=np.ceil(maximum*10)/10
    meshes=[]
    for row,(mean,actual) in enumerate(fields):
        for col,(label,field) in enumerate((("Mittelwert",mean),("Feldrealisierung",actual))):
            ax=fig.add_subplot(2,3,row*3+col+1,projection="3d")
            this=[]
            for level,color in zip(ISO_LEVELS,("#2563a6","#dc5534")):
                triangles=marching_tetrahedra(field,level)
                collection=Poly3DCollection(triangles,facecolor=color,edgecolor="none",alpha=.52,rasterized=True)
                ax.add_collection3d(collection)
                this.append({"level":level,"triangles":triangles})
                arrays[f"isosurface_t{PLOT_TIMES[row]:g}_{label}_level{level:g}"]=triangles
            meshes.append(this)
            ax.set(xlim=(0,L),ylim=(0,L),zlim=(0,L),xlabel="x",ylabel="y",zlabel="z",
                   title=f"{label}, t={PLOT_TIMES[row]:g}")
            ax.set_box_aspect((1,1,1));ax.view_init(elev=23,azim=-53)
            ax.set_xticks([0,np.pi,L],["0","π","2π"])
            ax.set_yticks([0,np.pi,L],["0","π","2π"])
            ax.set_zticks([0,np.pi,L],["0","π","2π"])
        ax=fig.add_subplot(2,3,row*3+3)
        im=ax.imshow(actual[:,:,GRID_N//2].T,origin="lower",extent=(0,L,0,L),
                     cmap="RdBu_r",vmin=-color_limit,vmax=color_limit)
        ax.set(xlabel="x",ylabel="y",title=f"Realisierung: Schnitt z=π, t={PLOT_TIMES[row]:g}")
        ax.set_xticks([0,np.pi,L],["0","π","2π"])
        ax.set_yticks([0,np.pi,L],["0","π","2π"])
    fig.colorbar(im,ax=fig.axes[2::3],shrink=.8,label="Geglätteter Feldwert φℓ")
    fig.suptitle("Bohmsches Skalarfeld im dreidimensionalen Raum\n"
                 "Isowert −0,075: blau; +0,075: rot. K=8, ℓ=0,55; periodischer Würfel.",fontsize=13)
    png=root/"figures"/(stem+".png");svg=png.with_suffix(".svg")
    fig.savefig(png,dpi=155);fig.savefig(svg,metadata={"Date":None});plt.close(fig)
    # The offline viewer uses a coarser, still physical 24^3 mesh of the same raw
    # field. Full 48^3 meshes remain available in NPZ and in the static figures.
    payload=[]
    for i,pair in enumerate(meshes):
        coarse_field=fields[i//2][i%2][::2,::2,::2]
        payload.append({"title":f"{'Mittelwert' if i%2==0 else 'Feldrealisierung'}, t={PLOT_TIMES[i//2]}",
                        "surfaces":[{"level":item["level"],"triangles":np.round(marching_tetrahedra(coarse_field,item["level"]),5).tolist()} for item in pair]})
    html="""<!doctype html><html lang="de"><meta charset="utf-8"><title>Bohmsches Skalarfeld: echte 3D-Isosurfaces</title>
<style>body{margin:0;background:#f5f7fa;color:#162238;font:16px system-ui}header{padding:20px}main{display:grid;grid-template-columns:1fr 1fr;gap:12px;padding:0 20px 20px}section{background:white;border-radius:12px;padding:10px}h2{font-size:16px;margin:8px}canvas{width:100%;height:400px;touch-action:none}p{margin:6px 0}</style>
<header><h1>Bohmsches Skalarfeld im dreidimensionalen Raum</h1><p>Ziehen dreht die tatsächlichen xyz-Isosurfaces; blau φℓ=−0,075, rot φℓ=+0,075.</p><p>Periodischer Würfel [0,2π]³; m=1, K=8, ℓ=0,55. Die dargestellte Feldrealisierung ist keine Wellenfunktion. Interaktive Geometrie: 24³; Rohfelder und statische Abbildungen: 48³.</p></header><main id="panels"></main><script>
const DATA=PAYLOAD;const L=2*Math.PI;
DATA.forEach(data=>{const box=document.createElement('section');box.innerHTML='<h2>'+data.title+'</h2><canvas width="720" height="540"></canvas>';document.querySelector('#panels').appendChild(box);const canvas=box.querySelector('canvas'),ctx=canvas.getContext('2d');let yaw=-.7,pitch=.5,last=null;
function transform(p){let x=p[0]-L/2,y=p[1]-L/2,z=p[2]-L/2;const a=x*Math.cos(yaw)-y*Math.sin(yaw),b=x*Math.sin(yaw)+y*Math.cos(yaw);return [360+56*a,270-56*(z*Math.cos(pitch)-b*Math.sin(pitch)),b*Math.cos(pitch)+z*Math.sin(pitch)];}
function draw(){ctx.clearRect(0,0,720,540);const triangles=[];for(const s of data.surfaces)for(const tri of s.triangles){const p=tri.map(transform);triangles.push({p,depth:p.reduce((a,b)=>a+b[2],0),color:s.level<0?'rgba(37,99,166,.50)':'rgba(220,85,52,.50)'});}triangles.sort((a,b)=>a.depth-b.depth);for(const t of triangles){ctx.beginPath();ctx.moveTo(t.p[0][0],t.p[0][1]);ctx.lineTo(t.p[1][0],t.p[1][1]);ctx.lineTo(t.p[2][0],t.p[2][1]);ctx.closePath();ctx.fillStyle=t.color;ctx.fill();}ctx.strokeStyle='#334155';ctx.lineWidth=1;for(let k=0;k<8;k++)for(let bit=0;bit<3;bit++)if(!(k&(1<<bit))){let a=[k&1?L:0,k&2?L:0,k&4?L:0],j=k|(1<<bit),b=[j&1?L:0,j&2?L:0,j&4?L:0];a=transform(a);b=transform(b);ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();}ctx.font='18px system-ui';ctx.fillStyle='#162238';[['x',[L,0,0]],['y',[0,L,0]],['z',[0,0,L]]].forEach(([s,p])=>{const q=transform(p);ctx.fillText(s,q[0]+7,q[1]);});}
canvas.onpointerdown=e=>{last=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId)};canvas.onpointerup=()=>last=null;canvas.onpointermove=e=>{if(last){yaw+=(e.clientX-last[0])*.01;pitch+=(e.clientY-last[1])*.01;last=[e.clientX,e.clientY];draw()}};draw();});
</script></html>""".replace("PAYLOAD",json.dumps(payload,separators=(",",":")))
    (root/"figures"/(stem+".html")).write_text(html,encoding="utf-8")
    metrics["plots"]={"grid_n":GRID_N,"times":list(PLOT_TIMES),"isosurface_levels":list(ISO_LEVELS),
                       "cutplane_z":float(np.pi),"color_limit":float(color_limit),
                       "mesh_triangle_counts":[[len(s["triangles"]) for s in pair] for pair in meshes],
                       "physical_coordinates":"x,y,z in [0,2pi]; duplicate periodic face used for marching tetrahedra",
                       "html":"offline canvas projection of 24^3 subsampled same physical field, rotatable by pointer",
                       "html_mesh_grid_n":GRID_N//2,"svg_3d_surfaces":"rasterized surfaces with vector axes and text"}
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag",default="20261017")
    parser.add_argument("--pilot",action="store_true")
    args=parser.parse_args()
    if not args.tag.replace("_","").replace("-","").isalnum(): parser.error("invalid tag")
    root=Path(__file__).resolve().parents[1];stem="quantum_field_"+args.tag
    paths=[root/"work"/(stem+s) for s in (".json",".npz",".md")]+[root/"figures"/(stem+s) for s in (".png",".svg",".html")]
    if not args.pilot and any(p.exists() for p in paths):parser.error("Historical output exists; choose a new tag.")
    start=time.perf_counter();model=GaussianField(8)
    generator=np.random.Generator(np.random.PCG64(SEED));xi=generator.standard_normal(len(model.n))
    arrays={"mode_vectors":model.n,"mode_kinds":model.kind,"xi":xi,"mu0":model.mu0,"pi0":model.pi0,
            "squeeze_r":model.r,"omega":model.omega,"smoothing":model.smoothing}
    metrics={"model":{"space_dimension":3,"spacetime_dimension":4,"L":float(L),"m":MASS,"hbar":1,"c":1,
                       "cutoffs":[4,6,8],"final_K":8,"real_modes":len(model.n),"ell":ELL,
                       "squeeze_low":SQUEEZE,"squeezed_condition":"|n|^2<=2, including zero mode",
                       "seed":SEED,"rng":"numpy.PCG64; largest K=8 sample, nested subsetting for smaller K",
                       "mean_modes":[{"n":n.tolist(),"kind":str(k),"mu0":float(m),"pi0":float(p)}
                                     for n,k,m,p in zip(model.n,model.kind,model.mu0,model.pi0) if m or p],
                       "foliation":"fixed torus time; no cutoff Lorentz-invariance claim"},
             "acceptance_thresholds":THRESHOLDS,"analytic":analytic_checks(model,xi,arrays),
             "finite_difference":finite_difference_checks(model,xi),"uv_cutoff":variance_cutoffs()}
    print(json.dumps({k:metrics[k] for k in ("analytic","finite_difference")},indent=2),flush=True)
    reference_xi=np.stack([xi,np.sin(np.arange(len(xi))*.71),np.cos(np.arange(len(xi))*.31)])
    initial=model.exact(0,reference_xi)
    rows=[]
    for dt in (.08,.04,.02):
        ts,trajectory=rk4(model,initial,dt)
        exact=np.stack([model.exact(t,reference_xi) for t in ts])
        error=float(np.max(abs(trajectory-exact)))
        rows.append({"dt":dt,"max_absolute_error":error})
        arrays[f"rk4_t_dt{dt:g}"]=ts
        # Dynamic low modes carry all errors; all 4913 modes were integrated.
        arrays[f"rk4_dynamic_modes_dt{dt:g}"]=trajectory[:,:,model.n2<=2]
        print("RK4",rows[-1],flush=True)
    metrics["rk4_convergence"]=rows
    metrics["rk4_orders"]=[float(np.log2(rows[i]["max_absolute_error"]/rows[i+1]["max_absolute_error"])) for i in (0,1)]
    vacuum=GaussianField(8,vacuum=True)
    vacuum_error=max(float(np.max(abs(vacuum.exact(t,xi)-vacuum.exact(0,xi)))) for t in np.linspace(0,FINAL_TIME,17))
    vacuum_velocity=max(float(np.max(abs(vacuum.velocity(t,vacuum.exact(t,xi))))) for t in np.linspace(0,FINAL_TIME,17))
    metrics["vacuum"]={"static_trajectory_error":vacuum_error,"guidance_velocity_max":vacuum_velocity,
                       "interpretation":"Stationary vacuum wavefunctional yields static actual mode coordinates, not temporal noise."}
    # Exact degree-two Gaussian moments from a deterministic Gauss-Hermite rule.
    gh_nodes, gh_weights=np.polynomial.hermite.hermgauss(5)
    gh_nodes*=np.sqrt(2);gh_weights/=np.sqrt(np.pi)
    mean_error=variance_error=twotime_error=0.0
    mu_reference,_,sigma_reference,*_=model.parameters(.23)
    q_reference=mu_reference[:,None]+sigma_reference[:,None]*gh_nodes
    for tau in np.linspace(0,FINAL_TIME,13):
        mu,_,sigma,*_=model.parameters(tau)
        cubature=mu[:,None]+sigma[:,None]*gh_nodes
        measured_mean=cubature@gh_weights
        measured_variance=((cubature-measured_mean[:,None])**2)@gh_weights
        measured_twotime=((cubature-mu[:,None])*(q_reference-mu_reference[:,None]))@gh_weights
        mean_error=max(mean_error,float(np.max(abs(measured_mean-mu))))
        variance_error=max(variance_error,float(np.max(abs(measured_variance-sigma**2))))
        twotime_error=max(twotime_error,float(np.max(abs(measured_twotime-sigma*sigma_reference))))
    metrics["deterministic_gaussian_moments"]={"hermite_nodes":5,"times":13,
        "mean_absolute_error":mean_error,"variance_absolute_error":variance_error,
        "bohmian_two_time_covariance_error":twotime_error,
        "two_time_interpretation":"Covariance of actual Bohmian coordinates; distinct from operator Wightman covariance."}
    arrays["gauss_hermite_nodes"]=gh_nodes;arrays["gauss_hermite_weights"]=gh_weights
    # Product quadrature verifies real basis and FFT normalization independently.
    small=GaussianField(2);coords=np.arange(7)*L/7
    xyz=np.array(np.meshgrid(coords,coords,coords,indexing="ij")).reshape(3,-1).T
    basis=basis_values(small,xyz,np.arange(len(small.n)))
    gram=basis.T@basis*VOLUME/len(xyz)
    metrics["spatial"]={"orthonormality_max":float(np.max(abs(gram-np.eye(len(small.n))))),"orthonormality_K":2,"quadrature_points":len(xyz)}
    q=model.exact(.73,xi)
    field,spectrum,imaginary=fft_field(model,q,smoothed=False)
    frequencies=np.fft.fftfreq(GRID_N)*GRID_N
    wavevectors=np.meshgrid(frequencies,frequencies,frequencies,indexing="ij")
    gradient=[(np.fft.ifftn(1j*k*spectrum)*GRID_N**3).real for k in wavevectors]
    energy_fft=.5*VOLUME*np.mean(sum(g*g for g in gradient)+MASS**2*field**2)
    energy_modes=.5*np.sum(model.omega**2*q*q)
    points=np.array([[.1,.3,1.7],[1.1,3.2,.8],[5.9,.2,4.8]])
    direct_values=basis_values(model,points,np.arange(len(model.n)))@q
    metrics["spatial"].update(fft_energy=float(energy_fft),modal_potential_energy=float(energy_modes),
                               fft_energy_relative_error=float(abs(energy_fft-energy_modes)/energy_modes),
                               fft_imaginary_max=imaginary,direct_basis_values=direct_values.tolist())
    arrays["spatial_gram_K2"]=gram
    # Parseval coupling gives exact integrated squared differences of sampled
    # smoothed cutoff fields, independent of plotting-grid interpolation.
    convergence=[]
    for tau in PLOT_TIMES:
        q=model.exact(tau,xi)
        for a,b in ((4,6),(6,8)):
            shell=(np.max(abs(model.n),axis=1)>a)&(np.max(abs(model.n),axis=1)<=b)
            l2sq=float(np.sum((q[shell]*model.smoothing[shell])**2))
            fa,_,_=fft_field(model,q,cutoff=a);fb,_,_=fft_field(model,q,cutoff=b)
            convergence.append({"time":tau,"Ka":a,"Kb":b,"L2_norm_difference":float(np.sqrt(l2sq)),
                                "L2_squared":l2sq,"fft_parseval_error":float(abs(VOLUME*np.mean((fb-fa)**2)-l2sq)),
                                "grid_sup_difference":float(np.max(abs(fb-fa)))})
    metrics["smoothed_sample_cutoff_convergence"]=convergence
    tail=metrics["uv_cutoff"][-1]["vacuum_point_variance_smoothed"]-metrics["uv_cutoff"][3]["vacuum_point_variance_smoothed"]
    metrics["smoothed_vacuum_K8_to_K16_point_variance_tail"]=float(tail)
    fields=[]
    for tau in PLOT_TIMES:
        q=model.exact(tau,xi);mu=model.parameters(tau)[0]
        mean,_,_=fft_field(model,mu);actual,_,_=fft_field(model,q)
        arrays[f"field_mean_t{tau:g}"]=mean;arrays[f"field_actual_t{tau:g}"]=actual
        arrays[f"mode_coordinates_t{tau:g}"]=q
        fields.append((mean,actual))
    arrays["spatial_axis"]=np.arange(GRID_N)*L/GRID_N
    sphere_axis=np.arange(GRID_N)*L/GRID_N
    sx,sy,sz=np.meshgrid(sphere_axis,sphere_axis,sphere_axis,indexing="ij")
    sphere_radius=1.6
    sphere_mesh=marching_tetrahedra((sx-np.pi)**2+(sy-np.pi)**2+(sz-np.pi)**2-sphere_radius**2,0)
    radius_error=float(np.max(abs(np.linalg.norm(sphere_mesh-np.pi,axis=2)-sphere_radius)))
    metrics["isosurface_calibration"]={"sphere_radius":sphere_radius,"grid_n":GRID_N,
                                        "triangles":len(sphere_mesh),"maximum_radius_error":radius_error}
    if args.pilot:
        print(json.dumps({k:metrics[k] for k in ("rk4_convergence","rk4_orders","spatial","uv_cutoff","smoothed_sample_cutoff_convergence","smoothed_vacuum_K8_to_K16_point_variance_tail")},indent=2),flush=True)
        print("Field ranges",[[[float(f.min()),float(f.max())] for f in pair] for pair in fields],flush=True)
        return
    th=THRESHOLDS;a=metrics["analytic"];fd=metrics["finite_difference"];fine=fd["finite_differences"][-1]
    checks={
        "real_fourier_orthonormality":metrics["spatial"]["orthonormality_max"]<th["orthonormality"],
        "fft_gradient_energy":metrics["spatial"]["fft_energy_relative_error"]<th["fft_energy_relative"],
        "analytic_schrodinger":a["schrodinger_ratio_max"]<th["analytic_pde_relative"],
        "analytic_continuity":a["continuity_log_ratio_max"]<th["analytic_pde_relative"],
        "analytic_quantum_force":a["guidance_force_max"]<th["finite_difference_force_absolute"],
        "finite_difference_schrodinger":fine["schrodinger_relative_residual"]<th["finite_difference_pde_relative"],
        "finite_difference_continuity":fine["continuity_relative_residual"]<th["finite_difference_continuity_relative"],
        "phase_gradient_guidance":fine["phase_guidance_absolute_error"]<th["finite_difference_guidance_absolute"],
        "quantum_potential_gradient":fine["quantum_potential_gradient_absolute_error"]<th["finite_difference_guidance_absolute"],
        "trajectory_acceleration_force":fine["trajectory_force_absolute_error"]<th["finite_difference_force_absolute"],
        "total_energy_conservation":a["energy_total_relative_drift"]<th["energy_relative"],
        "excitation_energy_conservation":a["energy_excitation_relative_drift"]<th["energy_relative"],
        "equivariant_flow_density":a["flow_log_probability_max"]<th["flow_log_probability_error"],
        "deterministic_gaussian_moments":max(mean_error,variance_error,twotime_error)<th["orthonormality"],
        "rk4_trajectory_accuracy":rows[-1]["max_absolute_error"]<th["rk4_fine_absolute"],
        "phase_adaptive_integral":fd["phase_vs_adaptive_integral_max"]<th["phase_quadrature_error"],
        "continuous_phase_branch":fd["continuous_phase_branch_max_change"]<th["phase_branch_jump"],
        "static_vacuum_guidance":max(vacuum_error,vacuum_velocity)<th["vacuum_static_error"],
        "smoothed_sample_convergence":max(r["L2_norm_difference"] for r in convergence if r["Ka"]==6)<th["smoothed_K6_K8_l2_difference"],
        "smoothed_variance_tail":tail<th["smoothed_K8_K16_variance_tail"],
        "raw_vacuum_variance_growth":all(metrics["uv_cutoff"][i+1]["vacuum_point_variance_raw"]>metrics["uv_cutoff"][i]["vacuum_point_variance_raw"] for i in range(5)),
        "physical_isosurface_sphere_calibration":radius_error<th["isosurface_sphere_radius_error"],
    }
    metrics["acceptance_checks"]=checks;metrics["all_checks_pass"]=all(checks.values())
    version=write_figures(root,stem,fields,arrays,metrics)
    paths[0].parent.mkdir(exist_ok=True)
    np.savez_compressed(paths[1],**arrays)
    metrics["provenance"]={"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"scipy":scipy.__version__,
                           "matplotlib":version,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                           "arrays_sha256":hashlib.sha256(paths[1].read_bytes()).hexdigest(),
                           "command":f"python work/quantum_field_model.py --tag {args.tag}","elapsed_seconds":time.perf_counter()-start}
    paths[0].write_text(json.dumps(metrics,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    report=f"""# Bohmsches Skalarfeld: Numerik und dreidimensionale Felder

Verbindlicher Lauf: **{args.tag}**. Alle {len(checks)} festen Prüfungen bestanden: **{all(checks.values())}**.
Reproduktion mit neuem Tag: `python work/quantum_field_model.py --tag REPRODUCTION`.

Das Modell verwendet ein reelles freies Skalarfeld auf dem dreidimensionalen Torus
mit L=2π, m=ħ=c=1, einer festen Zeitfoliation und endlichem Fouriercutoff K=8
(4913 reelle Moden). Die modalen Gaußwellenfunktionen bilden eine Wellenfunktional
auf dem Konfigurationsraum. Die Abbildungen zeigen dagegen die tatsächliche
räumliche Feldrealisierung φℓ(x,y,z) und deren Erwartungswert.

Squeezing r=.65 gilt genau für |n|²≤2 einschließlich der Nullmode, sonst r=0.
Anfangsmittel: μ0=1, μc(1,0,0)=2, μs(0,1,0)=1.5, μc(0,0,1)=−1.
Anfangsimpulse: πs(1,0,0)=1.25, πc(0,1,0)=−.75. Alle übrigen Mittelwerte sind null.
Ein PCG64-Gaußvektor mit Seed20261017 auf der K8-Basis wird für alle Cutoffs
durch Teilmengenbildung gekoppelt. Die konkreten Koeffizienten stehen in der NPZ.
Der räumliche Filter exp(−ell²|n|²/2), ell=.55, macht die dargestellten Felder glatt.
Der ungeglättete Vakuumwert hat mit dem Cutoff wachsende Punktvarianz; eine
punktweise gewöhnliche Kontinuumsfeldfunktion wird nicht behauptet.

Schrödingergleichung, Kontinuität, Quantenpotential und Führung werden analytisch
und unabhängig durch finite Differenzen geprüft. Die komplexe Gesamtphase benutzt
ein stetig entfaltetes Argument und wird gegen adaptive Zeitintegration geprüft.
Der exakte Gaußfluss q(t)=μ(t)+σ(t)ξ wird gegen RK4, seinen Dichte-Jakobian und die
Energieerhaltung getestet. Im vollständigen Vakuum ist der Bohmsche Modenfluss
stationär; ein Vakuumensemble ist kein künstlich zeitlich erneuertes Rauschen.

Die Isosurfaces in PNG, SVG und der interaktiven Offline-HTML haben echte xyz-
Koordinaten. Sie wurden aus demselben 48³-Feld durch Marching Tetrahedra gewonnen;
die interaktive Darstellung verwendet für flüssige Rotation ein 24³-Teilgitter.
Die SVG-Datei rasterisiert nur die 3D-Flächen, während Achsen und Text Vektoren bleiben.
Blau ist der Isowert −.075, rot +.075; beide Zeiten und beide Felder verwenden dieselben
Isowerte. Daneben stehen wirkliche z=π-Schnitte mit gemeinsamer Farbskala.

- Analytisches Schrödingerresiduum: {a['schrodinger_ratio_max']:.6e}.
- Finite-Differenzen-Schrödingerresiduum: {fine['schrodinger_relative_residual']:.6e}.
- Feinster RK4-Fehler: {rows[-1]['max_absolute_error']:.6e}.
- Logarithmischer Dichte-Jakobian-Fehler: {a['flow_log_probability_max']:.6e}.
- Vakuumsubtrahierte Erwartungsenergie: {a['energy_vacuum_subtracted']:.12g}.
- Numerischer geglätteter Varianztail K8 bis K16: {tail:.6e}.

QuelleSHA256: `{metrics['provenance']['source_sha256']}`.
RoharraySHA256: `{metrics['provenance']['arrays_sha256']}`.
Vollständige Parameter, Prüfschwellen und Werte: `{paths[0].name}`.
Es werden keine synthetischen Werte als experimentelle Messdaten ausgegeben.
"""
    paths[2].write_text(report,encoding="utf-8")
    print(json.dumps({"all_checks_pass":all(checks.values()),"checks":checks,"files":[str(p.relative_to(root)) for p in paths]},indent=2),flush=True)
    if not all(checks.values()):raise SystemExit(1)


if __name__=="__main__":main()
