"""Exact conformal trace obstruction on the self-consistent OR geometries.

This is not a full renormalized mode-sum computation.  It tests an exact,
state-independent necessary condition for any full conformal source, using
the complete four-dimensional curvature and a fixed effective Box(R) term.
Only xi=1/6 OR backgrounds are tested for self-consistency of the same field.
No assertion about a minimally coupled full source is inferred.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
BSTAR=4/345
BOUNDS={'independent_geometry_scaled':2e-8,'curvature_relative':2e-8,
        'FD_BoxR_scaled':2e-5,'leading_asymptotic_relative':2e-5,
        'next_order_relative':2e-5,'required_b_shift_relative':2e-5,
        'stable_vs_direct_relative':2e-9}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def symbolic_functions():
    r,w,L,d=sp.symbols('r w L delta',positive=True)
    aa=d+(2*d+sp.Rational(1,15))*w+(5*d+sp.Rational(7,15))*w*w
    bb=-d-6*d*w+(sp.Rational(2,5)-9*d)*w*w
    wr=(1-w)/r+L*(1-w)**2*aa/(r**3*w*w)
    J=(1-w)/(r*w)+L*(1-w)**2*bb/(3*r**3*w**3)
    dr=lambda f:sp.diff(f,r)+sp.diff(f,w)*wr
    A=sp.factor(w*(dr(J)/2+J*J/4)+wr*J/4)
    B=w*J/(2*r);C=wr/(2*r);D=(1-w)/r**2
    Ricci=sp.factor(-2*A-4*B-4*C+2*D)
    K=sp.factor(4*(A*A+2*B*B+2*C*C+D*D))
    Ricci2=sp.factor((A+2*B)**2+(A+2*C)**2+2*(D-B-C)**2)
    box=sp.factor(w*dr(dr(Ricci))+(2*w/r+(wr+w*J)/2)*dr(Ricci))
    Q=sp.factor(K-Ricci2)
    qs=sp.factor(Q+sp.Rational(4,345)*box)
    pQ=sp.Poly(sp.expand(Q*w**10),w)
    pBox=sp.Poly(sp.expand(box*w**10),w)
    pStar=sp.Poly(sp.expand(qs*w**9),w)
    checks={
       'Q_w_minus10':sp.factor(pQ.nth(0)-32*L**4*d**4/(81*r**12))==0,
       'Q_w_minus9':sp.factor(pQ.nth(1)-4*L**4*d**3*(44*d+1)/(81*r**12))==0,
       'BoxR_w_minus10':sp.factor(pBox.nth(0)+920*L**4*d**4/(27*r**12))==0,
       'BoxR_w_minus9':sp.factor(pBox.nth(1)+460*L**4*d**3*(12*d+1)/(81*r**12))==0,
       'unique_leading_cancellation_b':sp.factor(pQ.nth(0)+sp.Rational(4,345)*pBox.nth(0))==0,
       'exceptional_b_next_order':sp.factor(pStar.nth(0)-4*L**4*d**3*(84*d-1)/(243*r**12))==0,
       'conformal_next_order_strictly_nonzero':sp.factor(pStar.nth(0).subs(d,sp.Rational(1,60))-L**4/(32805000*r**12))==0,
    }
    # alpha=4L/3.  Solve the exact trace for the coefficient it would require.
    alpha=4*L/3
    shift=sp.factor(-(30*Ricci/alpha+qs)/box)
    shift_co=sp.factor(sp.limit(shift/w,w,0,dir='+'))
    checks['required_b_shift_coefficient']=sp.factor(shift_co-(84*d-1)/(2070*d))==0
    checks['conformal_required_b_shift_equals_bstar']=sp.factor(shift_co.subs(d,sp.Rational(1,60))-sp.Rational(4,345))==0
    # BoxR derived from the true divergence form, not a classical reduction.
    F,h,hp,j,Rp,Rpp=sp.symbols('F h hp j Rp Rpp',positive=True)
    box_general=Rpp/h+(2/(h*r)+j/(2*h)-hp/(2*h*h))*Rp
    checks['BoxR_divergence_conversion']=sp.factor(box_general.subs({h:1/w,hp:-wr/w**2,j:J,Rp:dr(Ricci),Rpp:dr(dr(Ricci))})-box)==0
    fun=sp.lambdify((r,w,L,d),(Ricci,K,Ricci2,Q,box,qs,shift,wr,J),'numpy',cse=True)
    return checks,fun,{'Q_minus10':str(sp.factor(pQ.nth(0))),
        'Q_minus9':str(sp.factor(pQ.nth(1))),
        'BoxR_minus10':str(sp.factor(pBox.nth(0))),
        'BoxR_minus9':str(sp.factor(pBox.nth(1))),
        'bstar_residual_minus9':str(sp.factor(pStar.nth(0))),
        'required_b_linear_shift':str(shift_co)}

def independent_background(row,raw,L):
    pre=row['prefix'];w=raw[pre+'_w'];r=raw[pre+'_r'];g=raw[pre+'_g'];d=11/60-row['xi']
    def rhs(z,y):
        invr2,gg=y
        a=d+(2*d+1/15)*z+(5*d+7/15)*z*z
        den=z*z*(1-z)+L*invr2*(1-z)**2*a
        return [-2*invr2*z*z/den,
                (4*z*(1-z)+L*invr2*(1-z)**2*((1/15-4*d)+(13/15-4*d)*z))/(3*den)]
    sol=solve_ivp(rhs,(float(w[0]),0),(float(r[0]**-2),float(g[0])),
         method='DOP853',dense_output=True,rtol=8e-13,atol=(2e-18,8e-15))
    if not sol.success or sol.t[-1]!=0:raise RuntimeError(sol.message)
    rr,gg=sol.sol(w);rr=rr**-.5
    err=float(max(np.max(abs(rr-r)/(1+abs(r))),np.max(abs(gg-g)/(1+abs(g)))))
    return sol,err

def diagnostics(rr,ww,L,d,fun):
    R,K,R2,Q,box,qs,shift,wr,J=[np.asarray(v) for v in fun(rr,ww,L,d)]
    alpha=4*L/3
    # The tuned sum is evaluated after symbolic cancellation of w^-10.
    full0=alpha*Q/30;full1=alpha*(Q+box)/30;fullstar=alpha*qs/30
    return {'r':rr,'w':ww,'Ricci':R,'Kretschmann':K,'Ricci_squared':R2,
            'Q_K_minus_Ricci2':Q,'BoxRicci':box,'Q_plus_bstar_BoxRicci':qs,
            'required_b_shift':shift,'required_b':BSTAR+shift,
            'full_trace_8piG_b0':full0,'full_trace_8piG_b1':full1,
            'full_trace_8piG_bstar':fullstar,'Einstein_trace':-R,
            'Einstein_trace_defect_b0':R+full0,'Einstein_trace_defect_b1':R+full1,
            'Einstein_trace_defect_bstar':R+fullstar,'w_r':wr,'logf_r':J}

def one_model(row,raw,L,fun):
    d=11/60-row['xi'];pre=row['prefix'];sol,err=independent_background(row,raw,L)
    w=raw[pre+'_w'];r=raw[pre+'_r'];ar=diagnostics(r,w,L,d,fun)
    curv=max(float(np.max(abs(ar['Ricci']-raw[pre+'_Ricci'])/(1+abs(raw[pre+'_Ricci'])))),
             float(np.max(abs(ar['Kretschmann']-raw[pre+'_Kretschmann'])/(1+abs(raw[pre+'_Kretschmann'])))))
    # Numerical fourth-derivative check: five-point derivatives of R(w) and
    # w_r(w), taken on an independently integrated background.
    wf=np.geomspace(.65,.004,73);rf=sol.sol(wf)[0]**-.5;df=diagnostics(rf,wf,L,d,fun)
    fd=[]
    for step in (.008,.004,.002,.001):
        vals={}
        for off in (-2,-1,0,1,2):
            z=wf*(1+off*step);rz=sol.sol(z)[0]**-.5
            vals[off]=diagnostics(rz,z,L,d,fun)
        Rw=(vals[-2]['Ricci']-8*vals[-1]['Ricci']+8*vals[1]['Ricci']-vals[2]['Ricci'])/(12*step*wf)
        Rww=(-vals[-2]['Ricci']+16*vals[-1]['Ricci']-30*vals[0]['Ricci']+16*vals[1]['Ricci']-vals[2]['Ricci'])/(12*(step*wf)**2)
        wrw=(vals[-2]['w_r']-8*vals[-1]['w_r']+8*vals[1]['w_r']-vals[2]['w_r'])/(12*step*wf)
        wr=df['w_r'];J=df['logf_r']
        boxfd=wf*wr*wr*Rww+(wf*wr*wrw+(2*wf/rf+(wr+wf*J)/2)*wr)*Rw
        rel=float(np.max(abs(boxfd-df['BoxRicci'])/(1+abs(df['BoxRicci']))))
        fd.append({'relative_w_step':step,'BoxR_scaled_error':rel})
        ar[f'FD_BoxR_step_{step:g}']=boxfd
    ar['FD_w']=wf;ar['FD_BoxR_exact']=df['BoxRicci']
    tw=raw[pre+'_tail_w'];tr=raw[pre+'_tail_r'];tail=diagnostics(tr,tw,L,d,fun)
    rD=row['endpoint']['rD'];a=L*d/rD**3
    coeffQ=32*a**4/81;coeffBox=-920*a**4/27
    coeffStar=4*L**4*d**3*(84*d-1)/(243*rD**12)
    ratios={'Q_leading_ratio':tail['Q_K_minus_Ricci2']*tw**10/coeffQ,
            'BoxR_leading_ratio':tail['BoxRicci']*tw**10/coeffBox,
            'bstar_next_order_ratio':tail['Q_plus_bstar_BoxRicci']*tw**9/coeffStar,
            'required_b_shift_ratio':tail['required_b_shift']/(BSTAR*tw)}
    ar.update({f'tail_{key}':val for key,val in tail.items()})
    ar.update(ratios)
    ar['tail_radius_gap']=raw[pre+'_tail_radius_gap']
    mask=(w>=1e-3)&(w<.3)
    direct=ar['Q_K_minus_Ricci2'][mask]+BSTAR*ar['BoxRicci'][mask]
    stable=ar['Q_plus_bstar_BoxRicci'][mask]
    directerr=float(np.max(abs(direct-stable)/(1+abs(stable))))
    # Necessary trace mismatch divided by Einstein scale; no full-RSET claim.
    mismatches={}
    for key in ('b0','b1','bstar'):
        ratio=np.abs(tail['Einstein_trace_defect_'+key])/np.abs(tail['Einstein_trace'])
        ar['tail_relative_Einstein_defect_'+key]=ratio
        mismatches[key]=float(ratio[-1])
    checks={
      'independent_background':err<BOUNDS['independent_geometry_scaled'],
      'curvature_matches_full_4D_archive':curv<BOUNDS['curvature_relative'],
      'BoxR_five_point_derivatives':fd[-1]['BoxR_scaled_error']<BOUNDS['FD_BoxR_scaled'],
      'BoxR_derivative_refinement':fd[-1]['BoxR_scaled_error']<fd[0]['BoxR_scaled_error'],
      'Q_and_BoxR_leading_coefficients':max(abs(ratios['Q_leading_ratio'][-1]-1),abs(ratios['BoxR_leading_ratio'][-1]-1))<BOUNDS['leading_asymptotic_relative'],
      'bstar_nonzero_next_order':abs(ratios['bstar_next_order_ratio'][-1]-1)<BOUNDS['next_order_relative'],
      'required_b_is_not_constant':abs(ratios['required_b_shift_ratio'][-1]-1)<BOUNDS['required_b_shift_relative'],
      'symbolic_cancellation_matches_direct_safe_grid':directerr<BOUNDS['stable_vs_direct_relative'],
      'fixed_b0_obstructed':mismatches['b0']>1e20,
      'fixed_b1_obstructed':mismatches['b1']>1e20,
      'fixed_bstar_obstructed':mismatches['bstar']>1e20,
    }
    checks={k:bool(v) for k,v in checks.items()}
    rep={'prefix':pre,'mass':row['mass'],'xi':row['xi'],'rD':rD,
         'independent_geometry_scaled_error':err,'curvature_archive_scaled_error':curv,
         'BoxR_finite_difference_convergence':fd,'stable_direct_scaled_error':directerr,
         'last_w':float(tw[-1]),'last_asymptotic_ratios':{k:float(v[-1]) for k,v in ratios.items()},
         'last_relative_Einstein_trace_defects':mismatches,
         'checks':checks}
    return rep,ar

def make_figures(models,arrays,tag):
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'svg.fonttype':'none'})
    colors=['#0072B2','#D55E00','#009E73'];figures={}
    def save(fig,name):
        for ext in ('png','svg'):
            p=ROOT/f'figures/fullrset_{tag}_{name}.{ext}'
            fig.savefig(p,dpi=190,bbox_inches='tight');figures[p.relative_to(ROOT).as_posix()]=sha(p)
        plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    a=arrays['m1_conformal'];w=a['tail_w']
    for key,lab,col,ls in [('tail_Einstein_trace','Einstein: −R','#333333','-'),
                          ('tail_full_trace_8piG_b0','Anomalie: b=0',colors[0],'--'),
                          ('tail_full_trace_8piG_b1','Anomalie: b=1',colors[1],'-.'),
                          ('tail_full_trace_8piG_bstar','Anomalie: b=4/345',colors[2],':')]:
        axs[0].loglog(w,np.abs(a[key]),label=lab,color=col,ls=ls)
    axs[0].set(xlabel='w = 1/h',ylabel='α × Betrag der Spur',title='M/√α=1: verschiedene Divergenzordnungen');axs[0].legend(fontsize=8)
    for m,col in zip(models,colors):
        a=arrays[m['prefix']]
        axs[1].loglog(a['tail_w'],a['tail_relative_Einstein_defect_bstar'],color=col,label=f"M/√α={m['mass']:g}")
    axs[1].set(xlabel='w = 1/h',ylabel='|R + 8πG Tr T| / |R|',title='Auch die führende Kompensation reicht nicht');axs[1].legend(fontsize=8)
    fig.suptitle('Exakter konformer Spuranomalie-Test der OR-Geometrie · ξ=1/6, α=1')
    save(fig,'trace_obstruction')
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    for m,col in zip(models,colors):
        a=arrays[m['prefix']];w=a['tail_w']
        axs[0].loglog(w,np.abs(a['tail_required_b_shift']),color=col,label=f"M/√α={m['mass']:g}")
        axs[1].loglog(w,np.abs(a['bstar_next_order_ratio']-1),color=col,label=f"M/√α={m['mass']:g}")
    w=arrays[models[0]['prefix']]['tail_w']
    axs[0].loglog(w,BSTAR*w,'k:',label='(4/345) w')
    axs[0].set(xlabel='w = 1/h',ylabel='|b erforderlich − 4/345|',title='Die erforderliche Kopplung ist ortsabhängig');axs[0].legend(fontsize=8)
    axs[1].set(xlabel='w = 1/h',ylabel='|(Q + b* □R) / (c₉ w⁻⁹) − 1|',title='Verifizierter Koeffizient nach Kompensation');axs[1].legend(fontsize=8)
    fig.suptitle('Jede feste endliche effektive □R-Kopplung scheitert an dieser OR-Endpunktfamilie')
    save(fig,'fixed_coupling')
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    for m,col in zip(models,colors):
        a=arrays[m['prefix']];fd=m['BoxR_finite_difference_convergence']
        axs[0].loglog([f['relative_w_step'] for f in fd],[f['BoxR_scaled_error'] for f in fd],'o-',color=col,label=f"M/√α={m['mass']:g}")
        axs[1].loglog(a['tail_w'],np.abs(a['BoxR_leading_ratio']-1),color=col,label=f"M/√α={m['mass']:g}")
    axs[0].set(xlabel='Relativer Ableitungsschritt in w',ylabel='Maximaler skalierter □R-Fehler',title='Unabhängige Fünfpunktableitung');axs[0].legend(fontsize=8)
    axs[1].set(xlabel='w = 1/h',ylabel='|□R / (c₁₀ w⁻¹⁰) − 1|',title='Vollständiger kovarianter □R-Term');axs[1].legend(fontsize=8)
    fig.suptitle('Ableitungs- und Endpunktkontrollen · keine abgeschnittene Modensumme')
    save(fig,'verification')
    return figures

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-tag',default='20261018_final');p.add_argument('--tag',required=True);p.add_argument('--pilot',action='store_true');args=p.parse_args()
    sys.stdout.reconfigure(encoding='utf-8')
    assert all(re.fullmatch('[A-Za-z0-9_-]+',z) for z in (args.input_tag,args.tag))
    out=ROOT/f'work/fullrset_numerics_{args.tag}.json'
    assert not any(out.with_suffix(ext).exists() for ext in ('.json','.npz','.md')),'Choose a fresh output tag'
    start=time.perf_counter();source=ROOT/f'work/fourdim_numerics_{args.input_tag}.json'
    data=json.loads(source.read_text(encoding='utf-8'))
    assert data['all_checks_pass'] and data['provenance']['arrays_sha256']==sha(source.with_suffix('.npz'))
    raw=np.load(source.with_suffix('.npz'),allow_pickle=False);L=data['model']['L_GNhbar_over_16pi']
    symbolic,fun,coeffs=symbolic_functions();models=[];arrays={}
    selected=[m for m in data['models'] if abs(m['xi']-1/6)<1e-15]
    if args.pilot:selected=selected[:1]
    for row in selected:
        rep,ar=one_model(row,raw,L,fun);models.append(rep);arrays[row['prefix']]=ar
        print(json.dumps({'model':row['prefix'],'passed':sum(rep['checks'].values()),'checks':len(rep['checks']),'last_defects':rep['last_relative_Einstein_trace_defects']}),flush=True)
    checks={**symbolic,**{m['prefix']+'_'+k:v for m in models for k,v in m['checks'].items()}}
    fullraw={pre+'_'+key:value for pre,ar in arrays.items() for key,value in ar.items()}
    checks['all_arrays_finite']=all(np.isfinite(v).all() for v in fullraw.values())
    np.savez_compressed(out.with_suffix('.npz'),**fullraw)
    figs=make_figures(models,arrays,args.tag) if not args.pilot else {}
    report={'scope':'Exact necessary conformal trace test on the approximate OR geometries; not a full renormalized mode sum or a newly solved backreaction geometry',
      'excluded_claims':['No conclusion for the minimally coupled full trace','No proof that every full quantum geometry is singular or regular','No interchange of infinite spectral sums and endpoint limit'],
      'alpha':4*L/3,'L':L,'xi':1/6,'effective_BoxR_coefficients':[0,1,BSTAR],
      'trace_equation':'-R = alpha/30 * (K-Ricci_ab_Ricci^ab+b*BoxR), Lambda=0',
      'all_fixed_b_obstruction':'For b!=4/345 the w^-10 term is nonzero. For b=4/345 the conformal OR family has nonzero w^-9 term L^4/(32805000 rD^12), whereas Einstein trace is O(w^-5).',
      'symbolic_coefficients':coeffs,'symbolic':symbolic,'acceptance_bounds':BOUNDS,
      'models':models,'checks':checks,'all_checks_pass':all(checks.values()),'raw_array_count':len(fullraw),
      'input_receipt':source.relative_to(ROOT).as_posix(),'input_receipt_sha256':sha(source),
      'input_arrays_sha256':sha(source.with_suffix('.npz')),'script_sha256':sha(__file__),
      'arrays_sha256':sha(out.with_suffix('.npz')),'figures':figs,
      'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform()},
      'elapsed_seconds':time.perf_counter()-start}
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# Exakter konformer Spuranomalie-Test','',f'Maßgeblicher Tag: `{args.tag}`; {sum(checks.values())}/{len(checks)} Prüfungen bestanden; {len(fullraw)} endliche Roharrays.','',
       'Dies ist keine vollständige renormierte Modensumme. Die Rechnung prüft eine exakte notwendige Bedingung jeder vollständigen konformen Quelle auf den drei konformen OR-Geometrien. α=1, L=3/4, M=1,2,5 und ξ=1/6.','',
       'Die statische OR-Endpunktfamilie erfüllt für keinen festen endlichen effektiven □R-Koeffizienten die vollständige konforme Einstein-Spurgleichung. Bei b≠4/345 bleibt ein Term w⁻¹⁰. Beim einzigen führend kompensierenden b=4/345 bleibt ein strikt positiver Term L⁴/(32805000 rD¹²) w⁻⁹. Die Einstein-Spur divergiert nur wie w⁻⁵.','',
       'Dieser Ausschluss betrifft die konkrete exakte OR-Familie einschließlich ihrer nächstführenden Koeffizienten. Er bestimmt weder die korrekte vollständige Rückwirkungsgeometrie noch das Schicksal ihrer Singularitäten. Für minimale Kopplung enthält die exakte Spur weitere zustandsabhängige Beiträge; darauf wird der konforme Ausschluss nicht übertragen.','',
       'Die Rechnung verwendet vollständige Krümmungsinvarianten und den kovarianten □R-Term ohne Ordnungsreduktion. Positive Endpunktabstände werden aus dem gebundenen OR-Archiv übernommen; bei b=4/345 wird die führende Aufhebung symbolisch vollzogen, um Gleitkomma-Auslöschung zu vermeiden.','',
       '```json',json.dumps({'symbolic_coefficients':coeffs,'models':models,'script_sha256':report['script_sha256'],'arrays_sha256':report['arrays_sha256']},ensure_ascii=False,indent=2),'```']
    out.with_suffix('.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'passed':sum(checks.values()),'checks':len(checks),'all_checks_pass':report['all_checks_pass'],'seconds':report['elapsed_seconds']}),flush=True)
    if not report['all_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
