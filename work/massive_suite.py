"""Verification and scientific figures for massive_bohm.py; no PDF generation."""
from massive_bohm import *
from scipy.interpolate import RectBivariateSpline


def angular_quantiles(run,curves):
    rates=np.full_like(curves,np.nan)
    for i,psi in enumerate(run['psi']):
        valid=np.isfinite(curves[:,i]);rr=curves[valid,i]
        p=CubicSpline(run['r'],psi,axis=0)(rr)
        rates[valid,i]=-2*np.imag(p[:,0].conj()*p[:,1])/(rr*np.sum(abs(p)**2,axis=1))
    angle=np.full_like(curves,np.nan)
    for i in range(len(curves)):
        valid=np.isfinite(rates[i]);angle[i,valid]=cumulative_simpson(rates[i,valid],x=run['times'][valid],initial=0)
    return rates,angle


def direct_guidance(run,starts,thetas):
    """Separate ODE of j/q with cubic spinor interpolation, not a CDF inversion.

    Proper time uses the full 3+1 current, including azimuthal velocity.
    """
    # Not-a-knot splines reproduce the initial data without imposing a
    # fictitious even time extension. A short polynomial radial extension is
    # used only by ODE stages beyond the terminating cutoff event.
    dr=run['r'][1]-run['r'][0];extra=run['r'][0]-dr*np.arange(10,0,-1)
    rr=np.r_[extra,run['r']]
    pp=np.concatenate([CubicSpline(run['r'][:6],run['psi'][:,:6,:],axis=1)(extra),run['psi']],axis=1)
    interp=[(RectBivariateSpline(run['times'],rr,pp[:,:,i].real),RectBivariateSpline(run['times'],rr,pp[:,:,i].imag)) for i in range(2)]
    def evaluate(t,r):return np.array([a.ev(t,r)+1j*b.ev(t,r) for a,b in interp])
    for it,ir in [(0,1000),(1,1100),(100,600),(300,100),(400,1)]:
        assert max(abs(evaluate(run['times'][it],run['r'][ir])-run['psi'][it,ir]))<2e-13
    results=[]
    for r0,theta in zip(starts,thetas):
        def rhs(t,y):
            r,phi,tau=y;p=evaluate(t,r)
            q=np.sum(abs(p)**2);z=p[0].conj()*p[1]
            return [-shift(r)+(abs(p[0])**2-abs(p[1])**2)/q,
              -2*z.imag/(r*q),2*np.sqrt(z.real**2+z.imag**2*np.cos(theta)**2)/q]
        def event(t,y):return y[0]-LEFT
        event.direction=-1;event.terminal=True
        def horizon(t,y):return y[0]-2*M
        horizon.direction=-1;horizon.terminal=False
        sol=solve_ivp(rhs,(0,run['times'][-1]),[r0,0,0],method='DOP853',events=[event,horizon],dense_output=True,rtol=2e-10,atol=2e-11,max_step=.02)
        assert sol.success
        valid=run['times']<=sol.t[-1];data=np.full((3,len(run['times'])),np.nan);data[:,valid]=sol.sol(run['times'][valid])
        results.append({'r':data[0],'phi':data[1],'proper_time':data[2],
          'cutoff_time':float(sol.t_events[0][0]) if len(sol.t_events[0]) else None,
          'horizon_time':float(sol.t_events[1][0]) if len(sol.t_events[1]) else None,
          'cutoff_phi':float(sol.y[1,-1]),'cutoff_proper_time':float(sol.y[2,-1])})
    return results


def flat_reference(r,t,mu):
    # Finite propagation speed keeps this packet away from the periodic seam.
    p=initial(r,mu);k=2*np.pi*np.fft.fftfreq(len(r),r[1]-r[0]);E=np.sqrt(k*k+mu*mu)
    f=np.fft.fft(p,axis=0);Hf=np.column_stack([k*f[:,0]+mu*f[:,1],mu*f[:,0]-k*f[:,1]])
    return np.fft.ifft(np.cos(E*t)[:,None]*f-1j*(np.sinc(E*t/np.pi)*t)[:,None]*Hf,axis=0)


def summary(run):
    return {'n':run['n'],'dt':run['dt'],'mu':run['mu'],'M':run['M'],'kappa':run['kappa'],
      'initial_norm':float(run['norm'][0]),'final_norm':float(run['norm'][-1]),
      'inner_outflow':float(run['absorbed'][-1]),'outer_outflow':float(run['outer_loss'][-1]),
      'max_norm_balance_error':float(max(abs(run['norm_balance']-run['norm'][0])))}


def plots(tag,run,curves,direct,geo,labels,thetas,massruns,conv):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
    files={};colors=plt.cm.viridis(np.linspace(.08,.9,len(labels)))
    def save(fig,name):
        for ext in ['png','svg']:
            p=ROOT/f'figures/massive_{tag}_{name}.{ext}';assert not p.exists();fig.savefig(p,bbox_inches='tight');files[p.relative_to(ROOT).as_posix()]=sha(p)
        plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    q=np.sum(abs(run['psi'])**2,axis=2);keep=run['r']<=18
    im=axs[0].pcolormesh(run['times'],run['r'][keep],q[:,keep].T,shading='auto',cmap='magma',vmin=0,vmax=.65,rasterized=True)
    fig.colorbar(im,ax=axs[0],label='Radiale Dichte q')
    axs[0].axhline(2,color='cyan',ls='--',lw=1,label='Horizont r = 2M')
    for row,c in zip(curves,colors):axs[0].plot(run['times'],row,color=c,lw=1)
    axs[0].set(xlabel='PG-Zeit T / M',ylabel='Radius r / M',title='Massiver 3+1-Dirac-Zustand, κ = 1');axs[0].legend(loc='upper right',fontsize=8)
    for lab,row,cl,c in zip(labels,curves,geo,colors):
        axs[1].plot(run['times'],row,color=c,label=f'Bohm ζ = {lab:g}')
        axs[1].plot(run['times'],cl['r'],color=c,ls=':',lw=1.2)
    axs[1].plot([],[],':',color='gray',label='Geodäte, gleiche Startgeschwindigkeit')
    axs[1].axhline(2,color='gray',ls='--',lw=1);axs[1].axhline(LEFT,color='black',ls=':',lw=1)
    axs[1].set(xlabel='PG-Zeit T / M',ylabel='Radius r / M',title='Quantilbahnen und klassische Vergleichsbahnen',ylim=(0,14));axs[1].legend(fontsize=8)
    save(fig,'transport')
    fig=plt.figure(figsize=(12,5.7),layout='constrained');az0=np.array([.15,1.4,2.6,3.8,5.])
    for panel,limit in [(1,13),(2,3)]:
        ax=fig.add_subplot(1,2,panel,projection='3d')
        uu,ww=np.meshgrid(np.linspace(0,2*np.pi,45),np.linspace(0,np.pi,25))
        ax.plot_surface(2*np.cos(uu)*np.sin(ww),2*np.sin(uu)*np.sin(ww),2*np.cos(ww),color='#729fbc',alpha=.2,linewidth=0)
        if panel==2:ax.plot_surface(LEFT*np.cos(uu)*np.sin(ww),LEFT*np.sin(uu)*np.sin(ww),LEFT*np.cos(ww),color='gray',alpha=.35,linewidth=0)
        for i,(row,c,theta) in enumerate(zip(direct,colors,thetas)):
            valid=np.isfinite(row['r']);r=np.r_[row['r'][valid],LEFT];ph=np.r_[row['phi'][valid],row['cutoff_phi']]+az0[i]
            x=r*np.sin(theta)*np.cos(ph);y=r*np.sin(theta)*np.sin(ph);z=r*np.cos(theta)
            visible=(abs(x)<=limit)&(abs(y)<=limit)&(abs(z)<=limit)
            ax.plot(x[visible],y[visible],z[visible],color=c,lw=1.8,label=f'ζ = {labels[i]:g}')
            if panel==1:ax.scatter(x[0],y[0],z[0],color=c,s=15)
        ax.set(xlim=(-limit,limit),ylim=(-limit,limit),zlim=(-limit,limit),xlabel='x / M',ylabel='y / M',zlabel='z / M')
        ax.set_box_aspect((1,1,1));ax.view_init(elev=23,azim=-65)
        ax.set_title('Bohmsche Bahnen auf PG-Raumschnitten' if panel==1 else 'Horizontpassage; Abbruch bei r = 0,5M',fontsize=11)
        if panel==1:ax.legend(loc='upper left',fontsize=8)
    fig.suptitle('Massiver Dirac-Testzustand: μM = 1, κ = 1, mⱼ = 1/2',fontsize=13)
    save(fig,'trajectories3d')
    fig,axs=plt.subplots(2,2,figsize=(11,7.2),layout='constrained')
    axs[0,0].plot(run['times'],run['norm'],label='Wahrscheinlichkeit r ≥ 0,5M')
    axs[0,0].plot(run['times'],run['absorbed'],label='Abfluss bei r = 0,5M')
    axs[0,0].plot(run['times'],run['norm_balance'],'k--',label='Summe');axs[0,0].legend(fontsize=8);axs[0,0].set(ylabel='Wahrscheinlichkeit',xlabel='T / M',title='Unabhängige Volumen- und Flussbilanz')
    for row,c,lab in zip(direct,colors,labels):axs[0,1].plot(run['times'],row['phi'],color=c,label=f'ζ = {lab:g}')
    axs[0,1].set(xlabel='T / M',ylabel='Δφ / rad',title='Winkelbewegung aus dem 4D-Strom');axs[0,1].legend(fontsize=8)
    for row in massruns:axs[1,0].plot(row['times'],row['median'],label=f'μM = {row["mu"]:g}')
    axs[1,0].axhline(2,color='gray',ls='--');axs[1,0].set(xlabel='T / M',ylabel='Medianradius r / M',title='Festes p₀M = −2; Startpolarisation variiert');axs[1,0].legend(fontsize=8)
    axs[1,1].loglog([.0025,.00125],[conv['time_coarse_medium']['max_L2_spinor'],conv['time_medium_fine']['max_L2_spinor']],'o-',label='Benachbarte Zeitgitter')
    axs[1,1].loglog([.0025,.00125],np.array([1,.25])*conv['time_coarse_medium']['max_L2_spinor'],'k:',label='Quadratische Referenz')
    axs[1,1].set(xlabel='Größeres ΔT / M',ylabel='max. L²-Spinorunterschied',title='Strang-Konvergenz');axs[1,1].legend(fontsize=8)
    save(fig,'verification')
    return files


def run_suite(tag,out,begin):
    assert not out.with_suffix('.npz').exists()
    print('Evolve kappa=1 refinements',flush=True)
    a=evolve(n=3951,dt=.0025,kappa=1,save_every=.02)
    b=evolve(n=3951,dt=.00125,kappa=1,save_every=.02)
    run=evolve(n=3951,dt=.000625,kappa=1,save_every=.02)
    conv={'time_coarse_medium':compare(a,b),'time_medium_fine':compare(b,run)}
    coarse_summary=[summary(a),summary(b)];del a,b
    fine=evolve(n=7901,dt=.000625,kappa=1,save_every=.02)
    conv['space_refinement']=compare(run,fine);del fine
    labels=np.array([.1,.25,.5,.75,.9]);thetas=np.array([.65,.95,np.pi/2,2.1,2.45])
    curves,cdf=quantile_curves(run,labels);rates,angles=angular_quantiles(run,curves)
    print('Integrate direct 3D guidance and classical geodesics',flush=True)
    direct=direct_guidance(run,curves[:,0],thetas);geo=geodesics(curves[:,0],run['times'])
    rerr=[];pherr=[]
    for i,d in enumerate(direct):
        valid=np.isfinite(curves[i])&np.isfinite(d['r'])&(curves[i]>.65)
        rerr.append(float(max(abs(curves[i,valid]-d['r'][valid]))))
        pherr.append(float(max(abs(angles[i,valid]-d['phi'][valid]))))
    print('Benchmarks and mass comparison',flush=True)
    zero=evolve(n=3951,dt=.00125,kappa=0,save_every=.02);zcurve,_=quantile_curves(zero,labels)
    massruns=[{'mu':MU,'times':run['times'],'median':curves[2],'summary':summary(run)}]
    for mu in [.5,2.]:
        mr=evolve(n=3951,dt=.00125,mu=mu,kappa=1,save_every=.02);mc,_=quantile_curves(mr,[.5])
        massruns.append({'mu':mu,'times':mr['times'],'median':mc[0],'summary':summary(mr)})
    massruns.sort(key=lambda x:x['mu'])
    ml=evolve(n=1976,dt=.01,mu=0,kappa=0,save_every=2.);mlerr=[]
    for t,p in zip(ml['times'],ml['psi']):mlerr.append(float(np.sqrt(simpson(np.sum(abs(p-massless_reference(ml['r'],t))**2,axis=1),x=ml['r']))))
    fl=evolve(n=3951,dt=.00125,mu=1,M=0,kappa=0,end=2,save_every=.1)
    flref=flat_reference(fl['r'][:-1],2.,1.)
    flerr=float(np.sqrt((fl['r'][1]-fl['r'][0])*np.sum(abs(fl['psi'][-1,:-1]-flref)**2)))
    table=[]
    for i,lab in enumerate(labels):
        table.append({'quantile':float(lab),'initial_r':float(curves[i,0]),'theta':float(thetas[i]),
          'horizon_time':direct[i]['horizon_time'],'cutoff_time':direct[i]['cutoff_time'],
          'cutoff_delta_phi':direct[i]['cutoff_phi'],'cutoff_proper_time':direct[i]['cutoff_proper_time'],
          'classical_horizon_time':crossing(run['times'],geo[i]['r'],2),
          'classical_cutoff_time':geo[i]['cutoff_time'],'radial_guidance_vs_quantile':rerr[i],
          'angular_guidance_vs_quantile':pherr[i]})
    q=np.sum(abs(run['psi'])**2,axis=2);z=run['psi'][:,:,0].conj()*run['psi'][:,:,1];supported=q>1e-10
    speed2=((abs(run['psi'][:,:,0])**2-abs(run['psi'][:,:,1])**2)**2+4*z.imag**2)[supported]/q[supported]**2
    checks={
      'normalized_initial_data':abs(run['norm'][0]-1)<1e-12,
      'volume_boundary_balance':max(abs(run['norm_balance']-1))<4e-7,
      'second_order_time_convergence':3.8<conv['time_coarse_medium']['max_L2_spinor']/conv['time_medium_fine']['max_L2_spinor']<4.2,
      'space_convergence':conv['space_refinement']['max_L2_spinor']<3e-8,
      'direct_radial_guidance':max(rerr)<5e-5,
      'direct_angular_guidance':max(pherr)<5e-5,
      'full_current_causal':max(speed2)<1+1e-12,
      'classical_energy_conserved':max(g['energy_drift'] for g in geo)<1e-9,
      'curved_massless_exact_benchmark':max(mlerr)<2e-7,
      'flat_massive_fourier_benchmark':flerr<3e-6,
      'no_outer_boundary_signal':max(abs(run['outer_loss']))<1e-20,
      'inward_inner_flux':max(run['boundary_current'][:,0])<=1e-14,
      'all_selected_quantiles_cross_horizon_and_cutoff':all(row['cutoff_time'] is not None and row['horizon_time'] is not None for row in table)}
    print(json.dumps({'checks':{k:bool(v) for k,v in checks.items()},'convergence':conv,'benchmarks':[max(mlerr),flerr],'radial_errors':rerr,'angular_errors':pherr},indent=2),flush=True)
    arrays={k:v for k,v in run.items() if isinstance(v,np.ndarray)}
    arrays.update({'quantile_labels':labels,'quantile_r':curves,'quantile_phi':angles,'quantile_phi_rate':rates,'theta':thetas,
      'direct_r':np.array([d['r'] for d in direct]),'direct_phi':np.array([d['phi'] for d in direct]),
      'direct_tau':np.array([d['proper_time'] for d in direct]),'classical_r':np.array([g['r'] for g in geo]),
      'kappa0_r':zero['r'],'kappa0_times':zero['times'],'kappa0_psi':zero['psi'][::100],
      'kappa0_profile_times':zero['times'][::100],'kappa0_norm':zero['norm'],
      'kappa0_absorbed':zero['absorbed'],'kappa0_quantile_r':zcurve})
    for mr in massruns:arrays[f'mass_{mr["mu"]:g}_median']=mr['median']
    np.savez_compressed(out.with_suffix('.npz'),**arrays)
    figures=plots(tag,run,curves,direct,geo,labels,thetas,massruns,conv)
    result={'schema':'massive-bohm-pg-v1','source_sha256':sha(ROOT/'work/massive_bohm.py'),'suite_sha256':sha(__file__),
      'model':'Fixed Schwarzschild 3+1 Dirac test-spinor; kappa=1, j=1/2, mj=1/2; no backreaction',
      'parameters':{'M':M,'mu':MU,'p0':P0,'center':RC,'halfwidth':WIDTH,'inner_cutoff':LEFT,'outer_cutoff':RIGHT,'end_time':END},
      'method':'Exact characteristic transport, degree-5 interpolation and exact local rotation; Strang splitting',
      'main':summary(run),'time_refinements':coarse_summary,'convergence':conv,'trajectory_table':table,
      'kappa0':summary(zero),'mass_comparison':[{'mu':x['mu'],'horizon_time':crossing(x['times'],x['median'],2),'summary':x['summary']} for x in massruns],
      'benchmarks':{'massless_curved_max_L2':max(mlerr),'massive_flat_final_L2':flerr,
        'max_classical_energy_drift':max(g['energy_drift'] for g in geo),'max_supported_equatorial_speed_squared':float(max(speed2))},
      'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':bool(all(checks.values())),
      'arrays_file':out.with_suffix('.npz').relative_to(ROOT).as_posix(),'arrays_sha256':sha(out.with_suffix('.npz')),
      'nan_policy':'Trajectory arrays are NaN after the inner-cutoff crossing; wavefunction and volume/flux arrays must be finite.',
      'figures':figures,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'seconds':time.perf_counter()-begin}
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'output':str(out),'seconds':result['seconds'],'all_checks_pass':result['all_checks_pass'],'main':result['main']},indent=2))
    assert result['all_checks_pass'],'Inspect failed checks; do not publish this run.'
