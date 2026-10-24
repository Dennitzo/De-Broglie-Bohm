"""Analytic regular black-bounce geometry and supporting null-source plots."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def invariants(r,a=1.,M=1.):
    R=np.sqrt(r*r+a*a);v=np.sqrt(2*M/R)
    c1=M*(2*r*r-a*a)/R**5;c2=-M*r*r/R**5
    c3=(M*(2*a*a-r*r)-a*a*R)/R**5;c4=(2*M*r*r+a*a*R)/R**5
    K=4*(c1*c1+2*c2*c2+2*c3*c3+c4*c4)
    nec=-2*a*a*(1+v)**2/R**4
    return R,v,K,nec
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tag',required=True);args=ap.parse_args()
    assert args.tag.replace('_','').replace('-','').isalnum()
    out=ROOT/f'work/nonsing_geometry_{args.tag}.json';assert not out.exists()
    r=np.linspace(-8,8,3201);R,v,K,nec=invariants(r)
    bound=4*(6+3*9)
    checks={'finite_arrays':all(np.isfinite(vv).all() for vv in [R,v,K,nec]),
      'positive_radius':min(R)>=1,'finite_throat_K':abs(K[len(r)//2]-16)<1e-13,
      'strict_ingoing_NEC_violation':max(nec)<0,'uniform_K_bound':max(K)<=bound}
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    plt.rcParams.update({'savefig.dpi':170})
    axs[0,0].plot(r,R,label='R = √(r² + a²)');axs[0,0].plot(r,abs(r),'k:',label='|r|')
    axs[0,0].axhline(1,color='gray',ls='--');axs[0,0].set(ylabel='Flächenradius R / M',title='Positives Minimum statt Zentrum');axs[0,0].legend(fontsize=9)
    axs[0,1].plot(r,K,color='#267b9c');axs[0,1].scatter([0],[16],color='black',s=20,label='K(0) M⁴ = 16')
    axs[0,1].set(ylabel='Kretschmann-Skalar K M⁴',title='Krümmung bleibt am Hals endlich');axs[0,1].legend(fontsize=9)
    axs[1,0].plot(r,1-v,label='dr/dT = 1 − v');axs[1,0].plot(r,-1-v,label='dr/dT = −1 − v')
    axs[1,0].axhline(0,color='gray',lw=.7);axs[1,0].axvspan(-np.sqrt(3),np.sqrt(3),alpha=.1,color='gray')
    axs[1,0].set(ylabel='Radiale Nullgeschwindigkeit',title='Beide Nullrichtungen laufen innen zu kleinerem r');axs[1,0].legend(fontsize=9)
    axs[1,1].plot(r,nec,color='#a85932');axs[1,1].axhline(0,color='gray',lw=.7)
    axs[1,1].set(ylabel='M² Gₐᵦ k₋ᵃ k₋ᵇ',title='Erforderliche zusätzliche Quelle verletzt die NEC')
    for ax in axs.flat:
        ax.set_xlabel('Signierte Koordinate r / M');ax.axvline(0,color='gray',ls=':',lw=.6)
        for x in [-np.sqrt(3),np.sqrt(3)]:ax.axvline(x,color='#719db4',ls='--',lw=.7)
    fig.suptitle('Vorgegebene reguläre Simpson–Visser-Geometrie: a / M = 1')
    figures={}
    for ext in ['png','svg']:
        p=ROOT/f'figures/nonsing_{args.tag}_geometry.{ext}';assert not p.exists();fig.savefig(p,bbox_inches='tight');figures[p.relative_to(ROOT).as_posix()]=sha(p)
    plt.close(fig)
    np.savez_compressed(out.with_suffix('.npz'),r=r,R=R,v=v,K=K,NEC=nec)
    report={'schema':'regular-black-bounce-geometry-v1','source_sha256':sha(__file__),
      'symbolic_receipt_sha256':sha(ROOT/'work/nonsing_symbolic_20261018.json'),
      'parameters':{'M':1,'a':1},'K_at_throat':16,'uniform_K_bound':bound,
      'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':bool(all(checks.values())),
      'arrays_file':out.with_suffix('.npz').relative_to(ROOT).as_posix(),'arrays_sha256':sha(out.with_suffix('.npz')),'figures':figures}
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');print(json.dumps(report,indent=2))
    assert report['all_checks_pass']
if __name__=='__main__':main()
