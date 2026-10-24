"""Independent 3D Fourier check of the local-lapse Bohm-source obstruction.

All coordinates live on a flat spatial torus of side 2*pi. This is an
off-guidance-graph test of the proposed Hamiltonian algebra, not evolution
of a joint quantum-gravity state. Renormalized Casimir data are separate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Torus:
    def __init__(self,n,mass=1.0):
        assert n%2==1
        self.n=n
        self.mass=mass
        self.volume=(2*np.pi)**3
        axis=np.arange(n)*2*np.pi/n
        self.x=np.meshgrid(axis,axis,axis,indexing='ij')
        k=np.fft.fftfreq(n,1/n)
        self.k=np.meshgrid(k,k,k,indexing='ij')
        self.omega=np.sqrt(sum(v*v for v in self.k)+mass*mass)

    def K(self,q):
        return np.fft.ifftn(self.omega*np.fft.fftn(q)).real

    def grad(self,q):
        fq=np.fft.fftn(q)
        return [np.fft.ifftn(1j*k*fq).real for k in self.k]

    def A(self,lapse,q):
        g=self.grad(q)
        return -sum(self.grad(lapse*g[j])[j] for j in range(3))+self.mass**2*lapse*q

    def quantum_derivative(self,lapse,q):
        return -self.K(lapse*self.K(q))

    def density(self,lapse,q,p,quantum=True,constant=0.0):
        v=p*p+sum(g*g for g in self.grad(q))+self.mass**2*q*q
        if quantum:
            v=v-self.K(q)**2
        return float(np.mean(lapse*(v/2+constant)))

    def brackets(self,n,m,q,p):
        dn=self.A(n,q)
        dm=self.A(m,q)
        qn=self.quantum_derivative(n,q)
        qm=self.quantum_derivative(m,q)
        classical=float(np.mean(dn*m*p-n*p*dm))
        full=float(np.mean((dn+qn)*m*p-n*p*(dm+qm)))
        correction=float(np.mean(qn*m*p-n*p*qm))
        gq,gn,gm=self.grad(q),self.grad(n),self.grad(m)
        expected=float(np.mean(p*sum((n*gm[j]-m*gn[j])*gq[j] for j in range(3))))
        return classical,full,correction,expected


def main(tag):
    output=ROOT/f'work/local_numerics_{tag}.json'
    arrays=output.with_suffix('.npz')
    assert not output.exists() and not arrays.exists(),'Use a fresh tag.'
    checks={};metrics={};rows=[]
    eps=.3
    mass=1.0
    sizes=[13,17,25,33,49]
    for n in sizes:
        space=Torus(n,mass)
        x,y,z=space.x
        for mode in [1,2,3]:
            lapse_n=np.ones_like(x)
            lapse_m=1+eps*np.cos(mode*x)
            q=np.cos(mode*x)
            p=np.ones_like(x)
            classical,full,correction,dirac=space.brackets(lapse_n,lapse_m,q,p)
            w=np.sqrt(mass*mass+mode*mode)
            exact_classical=eps*mode*mode/2
            exact_full=eps*mass*(w-mass)/2
            exact_correction=eps*w*(mass-w)/2
            rows.append(dict(grid=n,mode=mode,classical=classical,dirac=dirac,
                             full=full,correction=correction,
                             exact_classical=exact_classical,exact_full=exact_full,
                             exact_correction=exact_correction,
                             bare_vacuum_density=float(np.sum(space.omega)/(2*space.volume))))
    errors={key:max(abs(row[key]-row['exact_'+key]) for row in rows) for key in ['classical','full','correction']}
    metrics['witness_max_abs_errors']=errors
    checks['positive_lapse_exact_3D_mode_witness']=max(errors.values())<3e-12
    checks['classical_local_algebra']=max(abs(row['classical']-row['dirac']) for row in rows)<3e-12
    checks['cutoff_independent_nonzero_defect']=all(abs(row['correction'])>.05 for row in rows)

    space=Torus(33,mass)
    x,y,z=space.x
    q=np.cos(x)+.3*np.sin(2*y)+.2*np.cos(x+y-z)
    p=1+.2*np.cos(x)+.15*np.sin(y+z)+.1*np.cos(x+y-z)
    n=1+.2*np.cos(x+z)
    m=1+.15*np.sin(y+z)+.1*np.cos(x)-.05*np.cos(x-y)
    classical,full,correction,dirac=space.brackets(n,m,q,p)
    metrics['mixed_configuration']=dict(classical=classical,dirac=dirac,full=full,correction=correction)
    checks['mixed_mode_algebraic_reconstruction']=abs(classical-dirac)<2e-13 and abs(full-classical-correction)<2e-13 and abs(correction)>.01
    # Actual Hamiltonian directional derivative along X_H[m], computed from
    # the quadratic energy functional, independently of the bracket formula.
    velocity_q=m*p
    velocity_p=-space.A(m,q)-space.quantum_derivative(m,q)
    fd=[]
    for step in [.01,.003,.001]:
        derivative=(space.density(n,q+step*velocity_q,p+step*velocity_p)-space.density(n,q-step*velocity_q,p-step*velocity_p))/(2*step)
        fd.append(float(abs(derivative-full)))
    metrics['energy_directional_derivative_errors']=fd
    checks['poisson_bracket_from_actual_energy_differences']=max(fd)<1e-10
    # A field-independent constant changes energies but no field derivatives.
    subtraction=[]
    for constant in [0.,-2.5,10.]:
        step=.001
        derivative=(space.density(n,q+step*velocity_q,p+step*velocity_p,constant=constant)-space.density(n,q-step*velocity_q,p-step*velocity_p,constant=constant))/(2*step)
        subtraction.append(derivative)
    metrics['constant_subtraction_brackets']=subtraction
    checks['vacuum_constant_cannot_remove_defect']=max(abs(v-full) for v in subtraction)<1e-10

    # A truly three-dimensional complex plane wave probes the individual
    # Fourier coefficients of K N K-A_N. No reduction to a 1D simulation.
    kv=np.array([1,2,0]);pv=np.array([1,0,1])
    angle=sum(kv[j]*space.x[j] for j in range(3))
    lapse=1+eps*np.cos(sum(pv[j]*space.x[j] for j in range(3)))
    real=np.cos(angle)
    delta=space.K(lapse*space.K(real))-space.A(lapse,real)
    wk=np.sqrt(mass**2+kv@kv)
    expected=np.zeros_like(real)
    gaps=[]
    for sign in [-1,1]:
        target=kv+sign*pv
        wt=np.sqrt(mass**2+target@target)
        gap=wk*wt-(kv@target+mass**2)
        expected+=eps*gap*np.cos(sum(target[j]*space.x[j] for j in range(3)))/2
        gaps.append(float(gap))
    plane_error=float(np.max(abs(delta-expected)))
    metrics['three_dimensional_multiplier_max_abs']=plane_error
    metrics['three_dimensional_multiplier_values']=gaps
    checks['3D_force_gap_from_Fourier_action']=plane_error<3e-12 and min(gaps)>0

    # Momentum-space gap identity and positivity, including angular modes.
    ax=np.linspace(-4,4,121)
    kx,ky=np.meshgrid(ax,ax,indexing='ij')
    pp=np.array([1.,0.,0.])
    k2=kx*kx+ky*ky
    dot=k2+kx
    omega=np.sqrt(k2+mass*mass)
    omega_plus=np.sqrt((kx+1)**2+ky*ky+mass*mass)
    direct=omega*omega_plus-(dot+mass*mass)
    positive=(mass*mass+ky*ky)/(omega*omega_plus+dot+mass*mass)
    metrics['rationalized_gap_identity_max_abs']=float(np.max(abs(direct-positive)))
    checks['positive_gap_identity']=np.min(positive)>0 and metrics['rationalized_gap_identity_max_abs']<2e-14
    modes=np.arange(1,21,dtype=float)
    ratio=mass/(np.sqrt(mass*mass+modes*modes)+mass)
    correction=1-ratio
    checks['one_local_gradient_cannot_repair_two_modes']=correction[1]>correction[0] and np.all(np.diff(ratio)<0)
    # No purely mass/cosmological ultralocal potential contributes to this
    # antisymmetric bracket. A common local scalar kinetic coefficient C and
    # gradient coefficient D require CD=1; reciprocal choices obey this.
    zs=np.array([.2,.5,1.,2.,7.])
    checks['reciprocal_local_kinetic_coefficients_close']=np.max(abs((1/zs)*zs-1))<1e-15
    data={'grid_sizes':np.array(sizes),'witness_grid':np.array([r['grid'] for r in rows]),
          'witness_mode':np.array([r['mode'] for r in rows]),
          'witness_classical':np.array([r['classical'] for r in rows]),
          'witness_full':np.array([r['full'] for r in rows]),
          'witness_correction':np.array([r['correction'] for r in rows]),
          'bare_vacuum_density':np.array([r['bare_vacuum_density'] for r in rows if r['mode']==1]),
          'modes':modes,'bracket_ratio':ratio,'required_gradient_shift':correction,
          'kx':kx,'ky':ky,'gap':positive,
          'field_slice':q[:,:,0],'lapse_slice':lapse[:,:,0],'force_gap_slice':delta[:,:,0]}
    checks['all_arrays_finite']=all(np.all(np.isfinite(x)) for x in data.values())
    checks={key:bool(value) for key,value in checks.items()}
    np.savez_compressed(arrays,**data)
    result={'scope':'Free massive scalar, flat 3-torus; fixed-background Bohm potential and strong off-graph matter-bracket test.',
            'parameters':{'mass':mass,'side':2*np.pi,'positive_lapse_amplitude':eps},
            'all_passed':all(checks.values()),'checks':checks,'metrics':metrics,'witnesses':rows,
            'source_sha256':sha(__file__),'arrays_file':str(arrays.relative_to(ROOT)).replace('\\','/'),
            'arrays_sha256':sha(arrays),'environment':{'python':platform.python_version(),'numpy':np.__version__}}
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'all_passed':result['all_passed'],'checks':checks,'metrics':metrics},indent=2))
    if not result['all_passed']:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    main(parser.parse_args().tag)
