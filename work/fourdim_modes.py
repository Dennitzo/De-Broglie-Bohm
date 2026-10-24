"""Exact test-field multipole potentials on the self-consistent OR geometry.

These endpoint diagnostics do not compute a renormalized stress tensor or
choose a self-adjoint boundary condition. No raw angular cutoff is called
renormalization. Tagged artifacts are immutable.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-tag',default='20261018_final')
    p.add_argument('--tag',required=True)
    a=p.parse_args()
    assert all(re.fullmatch('[A-Za-z0-9_-]+',s) for s in [a.input_tag,a.tag])
    out=ROOT/f'work/fourdim_modes_{a.tag}.json'
    assert not out.exists() and not out.with_suffix('.npz').exists()
    source=ROOT/f'work/fourdim_numerics_{a.input_tag}.json'
    d=json.loads(source.read_text(encoding='utf8'))
    assert d['all_checks_pass']
    data=np.load(source.with_suffix('.npz'),allow_pickle=False)
    assert d['provenance']['arrays_sha256']==sha(source.with_suffix('.npz'))
    L=d['model']['L_GNhbar_over_16pi']
    raw={};checks={};models=[]
    for row in d['models']:
        pre=row['prefix'];xi=row['xi'];delta=11/60-xi
        ar=lambda key:data[pre+'_tail_'+key]
        w,r,g,x,R=ar('w'),ar('r'),ar('g'),ar('optical'),ar('Ricci')
        f=np.exp(g)*w**(-1/3)
        aa=delta+(2*delta+1/15)*w+(5*delta+7/15)*w*w
        wr=(1-w)/r+L*(1-w)**2*aa/(r**3*w*w)
        J=(1-w)/(r*w)+L*(1-w)**2*(-delta-6*delta*w+(2/5-9*delta)*w*w)/(3*r**3*w**3)
        Vgeo=f*(w*J+wr)/(2*r)
        base=f/r**2
        Vxi=f*xi*R
        raw.update({pre+'_w':w,pre+'_x':x,pre+'_Vgeo':Vgeo,
                    pre+'_angular_unit':base,pre+'_Vxi':Vxi})
        rD=row['endpoint']['rD'];C=np.exp(row['endpoint']['gD'])
        rad_co=C*L*delta/(3*rD**4)
        x_co=3*rD**3/(8*L*delta*np.sqrt(C))
        geo_error=float(abs(Vgeo[-1]*w[-1]**(7/3)/rad_co-1))
        x_error=float(abs(x[-1]/(x_co*w[-1]**(8/3))-1))
        checks[pre+'_radial_coefficient']=geo_error<2e-5
        checks[pre+'_optical_coefficient']=x_error<2e-5
        last=[]
        for ell in (0,1,5,20):
            V=Vgeo+ell*(ell+1)*base+Vxi
            scaled=x*x*V
            raw[pre+f'_V_l{ell}']=V
            raw[pre+f'_scaled_l{ell}']=scaled
            checks[pre+f'_l{ell}_scaled_limit']=abs(scaled[-1]+xi/8)<2e-7
            if ell:
                ratio=ell*(ell+1)*base/Vgeo
                coeff=3*ell*(ell+1)*rD*rD/(L*delta)
                checks[pre+f'_l{ell}_angular_subleading']=abs(ratio[-1]/(coeff*w[-1]**2)-1)<2e-5
                raw[pre+f'_angular_ratio_l{ell}']=ratio
            last.append({'ell':ell,'scaled_V_last':float(scaled[-1]),'limit':-xi/8})
        # Integrability in the minimal case follows from the analytic exponent
        # V~x^(-7/8), not from any finite numerical integration.
        models.append({'prefix':pre,'xi':xi,'endpoint_radial_coefficient_error':geo_error,
          'endpoint_optical_coefficient_error':x_error,'limits':last,
          'minimal_L1_exponent':-7/8 if xi==0 else None,
          'frobenius_powers_x':[(1-np.sqrt(1-xi/2))/2,(1+np.sqrt(1-xi/2))/2],
          'classification':'limit-circle at inner endpoint, for each fixed finite ell; analytic proof in manuscript'})
    checks['all_arrays_finite']=all(np.isfinite(v).all() for v in raw.values())
    checks={key:bool(value) for key,value in checks.items()}
    assert all(checks.values()),[k for k,v in checks.items() if not v]
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,3,figsize=(13.8,4.1),layout='constrained')
    colors=['#0072B2','#D55E00','#009E73']
    for pre,label,col in [('m1_minimal','M=1, ξ=0',colors[0]),
                          ('m1_conformal','M=1, ξ=1/6',colors[1])]:
        w=raw[pre+'_w'];xx=raw[pre+'_x']
        axs[0].loglog(xx,raw[pre+'_Vgeo'],color=col,label=label)
        axs[1].loglog(w,raw[pre+'_angular_ratio_l5'],color=col,label=label)
    axs[0].set(xlabel='Optischer Abstand x − xD [√α]',ylabel='Radialpotential α Vgeo',
               title='Radiale Rückstreuung: Vgeo ∝ (x − xD)⁻⁷ᐟ⁸')
    axs[0].legend(fontsize=8)
    axs[1].set(xlabel='w = 1/h',ylabel='Winkelpotential / Radialpotential',
               title='Feste Winkelmode l=5: Verhältnis ∝ w²')
    for i,m in enumerate((1,2,5)):
        pre=f'm{m}_conformal';w=raw[pre+'_w']
        residual=np.abs(raw[pre+'_scaled_l20']+1/48)
        axs[2].loglog(w,residual,color=colors[i],label=f'M={m}, l=20')
    axs[2].set(xlabel='w = 1/h',ylabel='|(x − xD)² V₂₀ + 1/48|',
               title='Konforme Kopplung: Grenzwert −1/48')
    axs[2].legend(fontsize=8)
    fig.suptitle('Exakte Modenpotentiale auf der OR-Geometrie · α=1 · keine renormierte Modensumme')
    figures=[]
    for ext in ('png','svg'):
        dst=ROOT/f'figures/fourdim_{a.tag}_modes.{ext}'
        assert not dst.exists()
        fig.savefig(dst,dpi=190,bbox_inches='tight');figures.append(dst)
    plt.close(fig)
    np.savez_compressed(out.with_suffix('.npz'),**raw)
    report={'scope':'Exact test-field potentials on an approximate self-consistent geometry; no spectral RSET, no Hawking flux, no global boundary choice',
      'input_receipt':source.relative_to(ROOT).as_posix(),'input_receipt_sha256':sha(source),
      'input_arrays_sha256':sha(source.with_suffix('.npz')),'script_sha256':sha(__file__),
      'arrays_sha256':sha(out.with_suffix('.npz')),'models':models,'checks':checks,
      'all_checks_pass':True,'raw_array_count':len(raw),
      'fixed_l_limit_is_not_uniform_in_l':True,
      'figures':{str(f.relative_to(ROOT)):sha(f) for f in figures}}
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'output':str(out),'checks':len(checks),'models':models},indent=2))

if __name__=='__main__':main()
