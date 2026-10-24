"""Flux-conserving test scattering through a specified finite exterior segment.

The segment is 0 <= proper distance <= 8m, with constant zero-potential
leads on both sides. This defines a regulated scattering experiment on
the computed potential, not a black-hole greybody factor or an RSET.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.interpolate import CubicSpline,CubicHermiteSpline
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def transmission(ell,rho,V,omega,rtol=3e-10,rho_deriv=None,V_deriv=None):
    p=CubicSpline(ell,rho) if rho_deriv is None else CubicHermiteSpline(ell,rho,rho_deriv)
    pot=CubicSpline(ell,V) if V_deriv is None else CubicHermiteSpline(ell,V,V_deriv)
    def rhs(l,Y):
        return np.exp(-p(l))*np.array([Y[1],(pot(l)-omega**2)*Y[0]])
    sol=solve_ivp(rhs,(ell[0],ell[-1]),np.array([1,-1j*omega]),
                  method='DOP853',rtol=rtol,atol=rtol*.01,dense_output=True)
    assert sol.success
    psi,chi=sol.y[:,-1]
    incoming=(psi-chi/(1j*omega))/2;outgoing=(psi+chi/(1j*omega))/2
    T=1/abs(incoming)**2;R=abs(outgoing/incoming)**2
    sample=sol.sol(np.linspace(ell[0],ell[-1],81))
    flux=np.imag(np.conjugate(sample[0])*sample[1])
    defect=float(np.max(abs(flux/(-omega)-1)))
    return float(R),float(T),defect

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='20261018')
    parser.add_argument('--input-tag',default='20261018_final')
    args=parser.parse_args()
    if not all(t.replace('_','').replace('-','').isalnum() for t in (args.tag,args.input_tag)):
        parser.error('Invalid tag')
    rawpath=ROOT/f'work/swave_numerics_{args.input_tag}.npz'
    dest=ROOT/f'work/swave_scattering_{args.tag}.json'
    if dest.exists():parser.error('Historical report exists; choose a new tag')
    raw=np.load(rawpath,allow_pickle=False);data={};models=[];checks={}
    r,p,u,v=sp.symbols('r rho u v',real=True)
    U=r/3*(1-u*u-v*v-(2*r+9/r)*u*v+6*p*u*u/r**2)
    V=(r+6/r)*u*v-6*p*u*u/r**2-r*U
    pot=sp.exp(2*p)*(U+u*v)/r
    dpot=sp.lambdify((r,p,u,v),sum(sp.diff(pot,z)*f for z,f in {r:u,p:v,u:U,v:V}.items()),'numpy')
    frequencies=np.geomspace(.015,.8,48)
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.1),layout='constrained')
    for m,col in zip((1,2,5),('#126782','#d97706','#7950a4')):
        pre=f'm{m}_';l=raw[pre+'outer_x_from_throat'];p=raw[pre+'outer_rho'];v=raw[pre+'outer_potential']
        pd=raw[pre+'outer_v'];vd=dpot(raw[pre+'outer_r'],p,raw[pre+'outer_u'],pd)
        values=np.array([transmission(l,p,v,f/m,rho_deriv=pd,V_deriv=vd) for f in frequencies])
        reduced=np.array([transmission(l[::2],p[::2],v[::2],f/m,rho_deriv=pd[::2],V_deriv=vd[::2]) for f in frequencies])
        refined=np.array([transmission(l,p,v,f/m,3e-12,rho_deriv=pd,V_deriv=vd) for f in frequencies])
        grid_error=float(np.max(abs(reduced[:,:2]-values[:,:2])))
        solver_error=float(np.max(abs(refined[:,:2]-values[:,:2])))
        zero=np.array([transmission(l,p,np.zeros_like(v),f/m) for f in frequencies[::12]])
        unit=float(np.max(abs(values[:,0]+values[:,1]-1)))
        checks.update({pre+'flux_unitarity':unit<2e-8,pre+'Wronskian':bool(np.max(values[:,2])<2e-7),
                       pre+'grid_refinement':grid_error<2e-6,pre+'solver_refinement':solver_error<2e-8,
                       pre+'zero_potential_transparent':bool(np.max(zero[:,0])<2e-18 and np.max(abs(zero[:,1]-1))<2e-8),
                       pre+'nonzero_backscatter':bool(np.max(values[:,0])>.1)})
        lapse_interpolant=CubicHermiteSpline(l,p,pd)
        optical=quad(lambda z:float(np.exp(-lapse_interpolant(z))),float(l[0]),float(l[-1]),epsabs=2e-10)[0]
        models.append({'mass':m,'proper_segment':[float(l[0]),float(l[-1])],
            'optical_segment_length':optical,'max_flux_unitarity_error':unit,
            'max_Wronskian_error':float(np.max(values[:,2])),
            'half_grid_probability_difference':grid_error,'refined_solver_difference':solver_error,
            'low_frequency_reflection':float(values[0,0]),'high_frequency_reflection':float(values[-1,0])})
        data.update({pre+'m_omega':frequencies,pre+'reflection':values[:,0],pre+'transmission':values[:,1],pre+'flux_defect':values[:,2]})
        axes[0].semilogx(frequencies,values[:,0],color=col,label=f'm/√α = {m}')
        axes[1].loglog(frequencies,np.maximum(values[:,2],1e-16),color=col)
    axes[0].set(xlabel='Dimensionslose Frequenz mω',ylabel='Reflexionswahrscheinlichkeit R',ylim=(-.02,1.02),title='Regulierter Test: Eigendistanz 0 bis 8m')
    axes[1].set(xlabel='Dimensionslose Frequenz mω',ylabel='Maximaler relativer Flussrest',title='Erhaltung des Wronskians')
    axes[0].legend(frameon=False)
    for ax in axes:
        ax.grid(alpha=.2)
        ax.xaxis.set_major_locator(FixedLocator([.02,.05,.1,.2,.5]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle('Radiale Rückstreuung mit ausdrücklich festgelegten Nullpotential-Anschlüssen',fontsize=11)
    fbase=ROOT/f'figures/swave_{args.tag}_scattering'
    for ext in ('png','svg'):fig.savefig(fbase.with_suffix('.'+ext),dpi=185)
    arrays=dest.with_suffix('.npz');np.savez_compressed(arrays,**data)
    report={'scope':'Finite exterior potential segment with zero-potential leads; endpoint jumps are part of this regulated potential and its reflection. Not a global singular-geometry scattering matrix, Hawking flux or spectral RSET',
       'source_npz':rawpath.relative_to(ROOT).as_posix(),'source_npz_sha256':sha(rawpath),
       'script_sha256':sha(__file__),'arrays_sha256':sha(arrays),'models':models,'checks':checks,
       'all_checks_pass':all(checks.values()),'frequencies_per_model':len(frequencies)}
    dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    if not report['all_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
