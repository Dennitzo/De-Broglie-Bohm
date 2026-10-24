"""Self-consistent homogeneous Einstein-Dirac interior, not a renormalized RSET.

Two equally occupied angular partners; zero longitudinal momentum.
Regular x=sqrt(B) chart integrates the spinor, A/(-Bdot), and proper time
all the way to x=0. Metric/curvature are evaluated only at positive B.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
C=1.;B0=1.;A0=1.;MU=1.;KAPPA=1.
GAMMAS=[0.,.1,.3,.6]
THETA=np.pi/3
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bloch(p):
    z=p[0].conj()*p[1]
    return np.array([2*z.real,2*z.imag,abs(p[0])**2-abs(p[1])**2])

def solve(gamma,mu=MU,kappa=KAPPA,C=C,B0=B0,A0=A0,rtol=2e-12,atol=3e-14):
    x0=np.sqrt(B0);q0=np.sqrt(2*C/B0-1);E=np.sqrt(mu*mu+(kappa/B0)**2)
    p0=np.array([np.sqrt((E+mu)/(2*E)),np.sign(kappa)*np.sqrt((E-mu)/(2*E))],complex)
    y0=np.r_[p0,A0/q0,0.,0.,0.,0.].astype(complex)
    def fun(x,y):
        p=y[:2];S=bloch(p);w=np.sqrt(2*C-x*x)
        dp=2j/w*np.array([mu*x*x*p[0]+kappa*p[1],kappa*p[0]-mu*x*x*p[1]])
        f=y[2].real
        return np.r_[dp,gamma*(kappa*S[0]+mu*x*x*S[2])/w**3,-2*x*x/w,
          -2*S[1]*x**3/(f*w*w),-2*S[0]/w,
          -2*x*x/w*np.sqrt(S[2]**2+S[0]**2*np.cos(THETA)**2)]
    sol=solve_ivp(fun,(x0,0.),y0,method='DOP853',rtol=rtol,atol=atol,dense_output=True,max_step=.025)
    assert sol.success
    B=np.geomspace(B0,1e-12*B0,801);x=np.sqrt(B);y=sol.sol(x)
    p=y[:2];f=y[2].real;t=y[3].real;S=bloch(p);w=np.sqrt(2*C-B);q=w/x
    f1=gamma*(kappa*S[0]+mu*B*S[2])/w**3
    f2=gamma*(2*mu*x*S[2]/w**3+3*x*(kappa*S[0]+mu*B*S[2])/w**5)
    assert np.all(f>0),'This run leaves the positive longitudinal metric sector.'
    A=f*q;ux=x*f1/f;wx=-B/w**2
    a1=ux+wx-1 # x A_x / A
    a2=B*f2/f+2*ux*(wx-1)-2*C*B/w**4-2*wx+2 # x² A_xx / A
    KA=(w*w*a2-(B+2*w*w)*a1)/4 # B³ A_tt / A
    KC=w*w*a1/2 # B³ H_A H_B
    Kscaled=4*KA**2+8*KC**2+24*C*C
    energy=mu*S[2]+kappa*S[0]/B
    rho=energy/(A*B*B);pt=kappa*S[0]/(2*A*B**3)
    scaled_source=gamma*x*(kappa*S[0]+mu*B*S[2])/(f*w)
    scaled_pt=gamma*kappa*S[0]*x/(2*f*w)
    constraint=(2*KC+2*C)-scaled_source
    angular=-(KA-C+KC)-scaled_pt
    R=gamma*mu*S[2]/(A*B*B)
    ricci2=gamma**2*(rho*rho+2*pt*pt)
    # Independent one-dimensional quadrature checks all integrated guidance
    # coordinates on the same dense spinor/geometry solution.
    guide_quad=np.array([quad(lambda xx:float(fun(xx,sol.sol(xx))[k].real),x0,0,epsabs=1e-12,epsrel=1e-12)[0] for k in [4,5,6]])
    result={'B':B,'x':x,'psi':p.T,'S':S.T,'f':f,'time':t,'A':A,'Bdot':-q,
      'bohm_z':y[4].real,'bohm_phi':y[5].real,'bohm_tau':y[6].real,
      'HA':-w*a1/(2*x**3),'HB':-w/x**3,'rho':rho,'pt':pt,
      'K':Kscaled/B**6,'B6K':Kscaled,'R':R,'ricci2':ricci2,
      'weyl2':Kscaled/B**6-2*ricci2+R*R/3,
      'constraint_scaled':constraint,'angular_einstein_scaled':angular,
      'gamma':gamma,'mu':mu,'kappa':kappa,'C':C,'B0':B0,'A0':A0,
      'endpoint_psi':sol.y[:2,-1],'endpoint_f':float(sol.y[2,-1].real),'endpoint_time':float(sol.y[3,-1].real),
      'endpoint_guidance':sol.y[4:7,-1].real,'guidance_quadrature':guide_quad,
      'norm_error':float(np.max(abs(np.sum(abs(p)**2,axis=0)-1))),
      'imaginary_auxiliary_error':float(np.max(abs(y[2:].imag))),
      'accepted_steps':len(sol.t)-1,'rtol':rtol,'atol':atol}
    return result

def terminal_time(C,B0):
    th=np.arcsin(np.sqrt(B0/(2*C)))
    return 2*C*(th-np.sin(th)*np.cos(th))

def describe(a):
    C=a['C'];B0=a['B0'];q0=np.sqrt(2*C/B0-1)
    L=(abs(a['mu'])*B0**1.5/3+abs(a['kappa'])*np.sqrt(B0))/(2*C-B0)**1.5
    return {'gamma':a['gamma'],'mu':a['mu'],'kappa':a['kappa'],
      'endpoint_time':a['endpoint_time'],'exact_endpoint_time':terminal_time(C,B0),
      'endpoint_f':a['endpoint_f'],'analytic_f_lower':a['A0']/q0-a['gamma']*L,
      'min_f':float(np.min(a['f'])),'spinor_norm_error':a['norm_error'],
      'constraint_scaled_max':float(np.max(abs(a['constraint_scaled']))),
      'angular_einstein_scaled_max':float(np.max(abs(a['angular_einstein_scaled']))),
      'last_B6K':float(a['B6K'][-1]),'min_B6K':float(np.min(a['B6K'])),
      'endpoint_bloch':bloch(a['endpoint_psi']).tolist(),
      'endpoint_z_phi_tau':a['endpoint_guidance'].tolist(),
      'guidance_quadrature_error':float(np.max(abs(a['endpoint_guidance']-a['guidance_quadrature']))),
      'min_energy_density':float(np.min(a['rho'])),'min_angular_pressure':float(np.min(a['pt'])),
      'accepted_steps':a['accepted_steps']}

def make_plots(runs,tag):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
    files={};colors=plt.cm.viridis(np.linspace(.05,.9,len(runs)))
    def save(fig,name):
        for ext in ['png','svg']:
            p=ROOT/f'figures/diracback_{tag}_{name}.{ext}';assert not p.exists();fig.savefig(p,bbox_inches='tight');files[p.relative_to(ROOT).as_posix()]=sha(p)
        plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    a=runs[0];axs[0].plot(a['time'],a['B'],color='black');axs[0].axvline(a['endpoint_time'],color='gray',ls=':')
    axs[0].set(xlabel='Eigenzeit t / C',ylabel='Kugelradius B / C',title='Identischer Kollapsradius für alle γ')
    for a,c in zip(runs,colors):
        lab=f'γ = {a["gamma"]:g}';axs[1].semilogx(a['B'],a['f'],color=c,label=lab)
        axs[2].loglog(a['B'],a['A'],color=c,label=lab)
    for ax in axs[1:]:ax.invert_xaxis();ax.set_xlabel('Kugelradius B / C');ax.legend(fontsize=8)
    axs[1].set(ylabel='A / √(2C/B − 1)',title='Selbstkonsistente Änderung von A')
    axs[2].set(ylabel='Längenskala A',title='Längsdehnung bis zum Endpunkt')
    fig.suptitle('Homogener Einstein–Dirac-Innenraum: μC = κ = 1, B₀/C = A₀ = 1')
    save(fig,'geometry')
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for a,c in zip(runs,colors):
        lab=f'γ = {a["gamma"]:g}';B=a['B']
        axs[0,0].loglog(B,a['K'],color=c,label=lab);axs[0,1].semilogx(B,a['B6K'],color=c,label=lab)
        if a['gamma']:
            axs[1,0].loglog(B,a['gamma']*a['rho'],color=c,label=lab)
            axs[1,1].semilogx(B,2*a['ricci2']/a['K'],color=c,label=lab)
    axs[0,0].loglog(runs[0]['B'],24/runs[0]['B']**6,'k:',label='Beweis: K ≥ 24C²/B⁶')
    axs[0,1].axhline(48,color='black',ls=':',label='Grenzwert 48C²')
    titles=['Krümmungsdivergenz','Skalierte Krümmung','Rückgekoppelte Energiedichte','Ricci-Beitrag relativ zur Gesamtkrümmung']
    labels=['K C⁴','B⁶K / C²','γρ C²','2 RₐᵦRᵃᵇ / K']
    for ax,title,label in zip(axs.flat,titles,labels):
        ax.invert_xaxis();ax.set(xlabel='Kugelradius B / C',ylabel=label,title=title);ax.legend(fontsize=8)
    save(fig,'curvature')
    fig=plt.figure(figsize=(11,4.7),layout='constrained');ax=fig.add_subplot(121,projection='3d')
    a=runs[-1];s=a['S'];uu,vv=np.meshgrid(np.linspace(0,2*np.pi,40),np.linspace(0,np.pi,24))
    ax.plot_wireframe(np.cos(uu)*np.sin(vv),np.sin(uu)*np.sin(vv),np.cos(vv),color='gray',alpha=.12,lw=.4)
    ax.plot(*s.T,color='#186e9c',lw=2);ax.scatter(*s[0],color='green',s=35,label='Anfang B = 1')
    ax.scatter(*bloch(a['endpoint_psi']),color='red',s=35,label='Grenzwert B → 0')
    ax.set(xlabel='S₁',ylabel='S₂',zlabel='S₃',title='Reguläre innere Spinorentwicklung');ax.set_box_aspect((1,1,1));ax.legend(fontsize=8)
    ax2=fig.add_subplot(122)
    for k in range(3):ax2.semilogx(a['B'],a['S'][:,k],label=f'S{k+1}')
    ax2.invert_xaxis();ax2.set(xlabel='Kugelradius B / C',ylabel='Blochkomponente',title='Gemeinsamer Spinverlauf für alle γ');ax2.legend()
    fig.suptitle('Blochraum der reduzierten Diracamplitude; keine räumliche Teilchenbahn')
    save(fig,'spin')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for a,c in zip(runs,colors):axs[0].plot(a['time'],a['bohm_z'],color=c,label=f'γ = {a["gamma"]:g}')
    a=runs[-1];axs[1].plot(a['time'],a['bohm_phi'],color='#186e9c');axs[1].plot(a['time'],-a['bohm_phi'],'--',color='#a95727')
    axs[2].plot(a['time'],a['bohm_tau'],color='#186e9c');axs[2].plot(a['time'],a['time'],'k:',label='τ = t')
    axs[0].set(xlabel='Geometrische Eigenzeit t / C',ylabel='z / C',title='Längsführung auf der rückgekoppelten Metrik');axs[0].legend(fontsize=8)
    axs[1].set(xlabel='Geometrische Eigenzeit t / C',ylabel='Δφ / rad',title='Entgegengesetzte Winkelströme der Partner')
    axs[2].set(xlabel='Geometrische Eigenzeit t / C',ylabel='Bahn-Eigenzeit τ / C',title='Kausale Führung bei θ = π/3');axs[2].legend(fontsize=8)
    fig.suptitle('Bohmsche Einzelmodenführung; Anfang z = φ = τ = 0')
    save(fig,'guidance')
    return files

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--tag',required=True);args=ap.parse_args()
    assert args.tag.replace('_','').replace('-','').isalnum()
    out=ROOT/f'work/diracback_numerics_{args.tag}.json';assert not out.exists() and not out.with_suffix('.npz').exists()
    start=time.perf_counter();runs=[solve(g) for g in GAMMAS];ref=solve(.6,rtol=2e-13,atol=3e-15)
    conv={'max_spinor':float(np.max(abs(runs[-1]['psi']-ref['psi']))),
      'max_f':float(np.max(abs(runs[-1]['f']-ref['f']))),
      'max_relative_B6K':float(np.max(abs(runs[-1]['B6K']/ref['B6K']-1)))}
    ml=solve(.6,mu=0);q=np.sqrt(2/ml['B']-1);exactA=(1-.6/2)*q+.6/2
    mle=float(np.max(abs(ml['A']/exactA-1)));vacerr=float(np.max(abs(runs[0]['B6K']-48)))
    summaries=[describe(a) for a in runs]
    checks={
      'all_norms_conserved':all(a['norm_error']<2e-11 for a in runs),
      'all_analytic_positive_metric_bounds':all(row['analytic_f_lower']>0 for row in summaries),
      'all_actual_f_above_proven_bounds':all(row['endpoint_f']>=row['analytic_f_lower'] for row in summaries),
      'proper_endpoint_time':all(abs(a['endpoint_time']-terminal_time(1,1))<1e-12 for a in runs),
      'independent_angular_einstein_equation':all(row['angular_einstein_scaled_max']<1e-11 for row in summaries),
      'hamiltonian_constraint':all(row['constraint_scaled_max']<1e-11 for row in summaries),
      'kretschmann_lower_bound':all(np.min(a['B6K'])>=24-1e-11 for a in runs),
      'endpoint_asymptotic_coefficient':all(abs(a['B6K'][-1]-48)<1e-4 for a in runs),
      'vacuum_schwarzschild':vacerr<1e-11,
      'exact_massless_backreaction':mle<1e-10,
      'solver_tolerance_convergence':max(conv.values())<1e-10,
      'guidance_independent_quadrature':all(row['guidance_quadrature_error']<1e-11 for row in summaries),
      'full_bohm_current_causal':all(np.max(a['S'][:,1]**2+a['S'][:,0]**2*np.sin(THETA)**2)<1+1e-12 for a in runs),
      'all_raw_arrays_finite':all(np.isfinite(v).all() for a in runs for v in a.values() if isinstance(v,np.ndarray))}
    arrays={}
    for i,a in enumerate(runs):
        arrays.update({f'g{i}_{k}':v for k,v in a.items() if isinstance(v,np.ndarray)})
    arrays.update({'massless_B':ml['B'],'massless_A':ml['A'],'massless_exact_A':exactA})
    np.savez_compressed(out.with_suffix('.npz'),**arrays)
    figures=make_plots(runs,args.tag)
    result={'schema':'homogeneous-einstein-dirac-v1','source_sha256':sha(__file__),
      'scope':'Classical two-field homogeneous KS interior; zero longitudinal momentum; no absolute renormalized RSET or exterior matching',
      'parameters':{'C':C,'B0':B0,'A0':A0,'mu':MU,'kappa':KAPPA,'gammas':GAMMAS,'theta':THETA},
      'regular_chart':'x=sqrt(B), f=A/sqrt(2C/B-1); smooth integration to x=0',
      'runs':summaries,'convergence':conv,'vacuum_scaled_curvature_error':vacerr,'massless_relative_metric_error':mle,
      'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':bool(all(checks.values())),
      'arrays_file':out.with_suffix('.npz').relative_to(ROOT).as_posix(),'arrays_sha256':sha(out.with_suffix('.npz')),
      'figures':figures,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'seconds':time.perf_counter()-start}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:result[k] for k in ['runs','convergence','checks','seconds']},indent=2))
    assert result['all_checks_pass']

if __name__=='__main__':main()
