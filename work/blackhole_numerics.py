"""Massless real scalar Regge--Wheeler scattering on fixed Schwarzschild m=1.

No imports of older models. Units G=hbar=c=1. This is a stationary test-field
Hawking-flux calculation, not a self-consistent evaporating black-hole metric.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import scipy
from scipy.integrate import solve_ivp, simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.special import lambertw

HORIZON_EPS=1e-5
HORIZON_ORDER=3
OUTER_ORDER=12
OUTER_WAVELENGTH=20.
OUTER_MINIMUM=60.
RTOL=2e-10
L_MAX=5
FREQUENCY_STEP=.0025
FREQUENCIES=np.arange(1,401)*FREQUENCY_STEP
BOUNDS={
    "probability_range_tolerance":2e-7,
    "unitarity_absolute":2e-7,
    "outer_basis_flux_absolute":2e-8,
    "horizon_flux_absolute":2e-8,
    "well_conditioned_Wronskian_relative":2e-7,
    "scaled_Wronskian":2e-8,
    "outer_radius_relative":2e-7,
    "horizon_epsilon_relative":2e-7,
    "tolerance_relative":2e-7,
    "matching_order_relative":2e-7,
    "horizon_order_relative":2e-7,
    "low_frequency_ratio_error":1e-3,
    "power_grid_relative":2e-5,
    "power_interpolation_relative":2e-5,
    "last_partial_wave_relative":2e-7,
    "retained_frequency_tail_relative":1e-6,
}


def tortoise(r):
    return r+2*math.log(r/2-1)


def radius(x):
    if x<600.:
        return 2*(1+float(lambertw(math.exp(x/2-1)).real))
    rr=x-2*math.log(x/2-1)
    for _ in range(3):
        rr-=(rr+2*math.log(rr/2-1)-x)*(1-2/rr)
    return rr


def potential(r,ell):
    return (1-2/r)*(ell*(ell+1)/r**2+2/r**3)


def horizon_series(omega,ell,eps,order):
    """psi=e^(-i omega x) sum a_n ((r-2)/2)^n; a_0=1."""
    angular=ell*(ell+1)
    coeff=[1.+0j]
    for n in range(order):
        value=(2*n*n-n-12j*omega*n-angular-1)*coeff[n]
        if n>=1:value+=((n-1)*(n-2)-12j*omega*(n-1)-angular)*coeff[n-1]
        if n>=2:value+=(-4j*omega*(n-2))*coeff[n-2]
        coeff.append(-value/((n+1)*(n+1-4j*omega)))
    s=sum(c*eps**n for n,c in enumerate(coeff))
    ds=sum(n*coeff[n]*eps**(n-1) for n in range(1,len(coeff)))
    f=eps/(1+eps)
    # Common horizon phase is set to one; this does not alter any flux.
    return s,-1j*omega*s+f*ds/2


def outer_series(omega,ell,r,order):
    """Outgoing Jost basis e^(i omega x) sum c_n/r^n; c_0=1."""
    angular=ell*(ell+1)
    coeff=[1.+0j]
    for n in range(order):
        value=(n*(n+1)-angular)*coeff[n]
        if n>=1:value-=2*n*n*coeff[n-1]
        coeff.append(value/(2j*omega*(n+1)))
    s=sum(c/r**n for n,c in enumerate(coeff))
    ds=sum(-n*coeff[n]/r**(n+1) for n in range(1,len(coeff)))
    phase=np.exp(1j*omega*tortoise(r))
    outgoing=phase*s
    deriv=phase*(1j*omega*s+(1-2/r)*ds)
    return outgoing,deriv


def scatter(omega,ell,rtol=RTOL,eps=HORIZON_EPS,horizon_order=HORIZON_ORDER,
            outer_scale=OUTER_WAVELENGTH,outer_min=OUTER_MINIMUM,outer_order=OUTER_ORDER):
    rh=2*(1+eps); rout=max(outer_min,outer_scale/omega)
    xh=tortoise(rh);xo=tortoise(rout)
    initial=np.array(horizon_series(omega,ell,eps,horizon_order),dtype=complex)
    def rhs(x,y):
        r=radius(x)
        return np.array([y[1],(potential(r,ell)-omega*omega)*y[0]])
    sol=solve_ivp(rhs,(xh,xo),initial,method="DOP853",rtol=rtol,atol=rtol*1e-3)
    if not sol.success:raise RuntimeError(sol.message)
    psi,derivative=sol.y[:,-1]
    outgoing,dp=outer_series(omega,ell,rout,outer_order)
    incoming=outgoing.conjugate();dm=dp.conjugate()
    determinant=incoming*dp-outgoing*dm
    ain=(psi*dp-derivative*outgoing)/determinant
    aout=(incoming*derivative-dm*psi)/determinant
    transmission=1/abs(ain)**2
    reflection=abs(aout/ain)**2
    horizon_flux=float(np.imag(initial[0].conjugate()*initial[1]))
    endpoint_flux=float(np.imag(psi.conjugate()*derivative))
    outer_flux=float(np.imag(outgoing.conjugate()*dp))
    scaled_wronskian=abs(endpoint_flux-horizon_flux)/(omega+abs(psi)*abs(derivative))
    pathflux=np.imag(sol.y[0].conjugate()*sol.y[1])
    pathscale=omega+np.abs(sol.y[0])*np.abs(sol.y[1])
    return {"omega":float(omega),"ell":int(ell),"transmission":float(transmission),
            "reflection":float(reflection),"unitarity_residual":float(abs(reflection+transmission-1)),
            "horizon_flux_over_minus_omega":float(-horizon_flux/omega),
            "outer_flux_over_omega":float(outer_flux/omega),
            "relative_endpoint_wronskian_error":float(abs(endpoint_flux-horizon_flux)/omega),
            "scaled_endpoint_wronskian_error":float(scaled_wronskian),
            "max_relative_path_wronskian_error":float(np.max(abs(pathflux-horizon_flux))/omega),
            "max_scaled_path_wronskian_error":float(np.max(abs(pathflux-horizon_flux)/pathscale)),
            "ain_real":float(ain.real),"ain_imag":float(ain.imag),
            "aout_real":float(aout.real),"aout_imag":float(aout.imag),
            "r_outer":rout,"x_horizon":xh,"nfev":sol.nfev}


def scatter_task(task):
    ell,omega=task
    return scatter(omega,ell)


def scan(workers):
    jobs=[(ell,float(w)) for ell in range(L_MAX+1) for w in FREQUENCIES]
    if workers==1:
        return list(map(scatter_task,jobs))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(scatter_task,jobs,chunksize=12))


def convergence_controls():
    witnesses=[(0,.01),(0,.1),(1,.1),(2,.2),(3,.5),(5,1.)]
    controls=[]
    for ell,omega in witnesses:
        baseline=scatter(omega,ell)
        target=baseline["transmission"]
        row={"ell":ell,"omega":omega,"baseline_transmission":target}
        specifications={
            "outer_radius":[{"outer_scale":s,"outer_min":3*s} for s in (10.,20.,40.)],
            "horizon_epsilon":[{"eps":e} for e in (1e-3,1e-4,1e-5,1e-6)],
            "ODE_tolerance":[{"rtol":t} for t in (2e-8,2e-10,2e-12)],
            "outer_matching_order":[{"outer_order":n} for n in (6,9,12,15)],
            "horizon_series_order":[{"horizon_order":n} for n in (0,1,3,4)]}
        for kind,settings in specifications.items():
            values=[]
            for params in settings:
                result=scatter(omega,ell,**params)
                values.append({"parameters":params,"transmission":result["transmission"],
                               "relative_difference_from_baseline":abs(result["transmission"]-target)/target})
            row[kind]=values
        controls.append(row)
    return controls


def build_results(rows,controls):
    frequency=np.r_[0.,FREQUENCIES]
    shape=(L_MAX+1,len(frequency))
    gamma=np.zeros(shape);reflection=np.ones(shape)
    fields=["unitarity_residual","horizon_flux_over_minus_omega","outer_flux_over_omega",
            "relative_endpoint_wronskian_error","scaled_endpoint_wronskian_error",
            "max_relative_path_wronskian_error","max_scaled_path_wronskian_error",
            "ain_real","ain_imag","aout_real","aout_imag","r_outer","nfev"]
    raw={key:np.zeros(shape) for key in fields}
    for row in rows:
        i=row["ell"];j=int(round(row["omega"]/FREQUENCY_STEP))
        gamma[i,j]=row["transmission"];reflection[i,j]=row["reflection"]
        for key in fields:raw[key][i,j]=row[key]
    thermal=np.zeros_like(frequency)
    thermal[1:]=frequency[1:]/np.expm1(8*np.pi*frequency[1:])/(2*np.pi)
    partial=(2*np.arange(L_MAX+1)+1)[:,None]*gamma*thermal
    spectrum=partial.sum(axis=0)
    powers=simpson(partial,x=frequency,axis=1)
    total=float(powers.sum())
    quadrature=[]
    for stride in (4,2,1):
        grid=frequency[::stride];values=partial[:,::stride]
        pp=simpson(values,x=grid,axis=1)
        quadrature.append({"frequency_step":float(grid[1]-grid[0]),"sample_points":len(grid),
            "partial_powers":pp.tolist(),"total_power":float(pp.sum()),
            "relative_difference_from_finest":float(abs(pp.sum()-total)/total)})
    interpolated=[]
    for stride in (4,2,1):
        grid=frequency[::stride];interp=PchipInterpolator(grid,gamma[:,::stride],axis=1)
        degeneracy=2*np.arange(L_MAX+1)+1
        def fn(w):
            if w==0:return 0.
            return float(np.dot(degeneracy,interp(w)))*w/np.expm1(8*np.pi*w)/(2*np.pi)
        result,error=quad(fn,0,1,points=grid[1:-1],epsabs=2e-13,epsrel=2e-10,limit=800)
        interpolated.append({"frequency_step":float(grid[1]-grid[0]),"power_PCHIP_Gamma_then_quad":result,
                             "quad_internal_error_estimate":error,"relative_difference_from_Simpson_finest":abs(result-total)/total})
    tail=(L_MAX+1)**2/(2*np.pi)*np.exp(-8*np.pi)/(1-np.exp(-8*np.pi))*(1/(8*np.pi)+1/(8*np.pi)**2)
    low=[]
    for w in (.0001,.0002,.0005,.001,.002):
        result=scatter(w,0,rtol=2e-12)
        result["ratio_to_16omega_squared"]=result["transmission"]/(16*w*w)
        result["absorption_cross_section_over_horizon_area"]=result["ratio_to_16omega_squared"]
        low.append(result)
    high_condition=[r for r in rows if r["transmission"]>=1e-4]
    convergence_max={
        "outer_radius":max(c["outer_radius"][-1]["relative_difference_from_baseline"] for c in controls),
        "horizon_epsilon":max(c["horizon_epsilon"][-1]["relative_difference_from_baseline"] for c in controls),
        "ODE_tolerance":max(c["ODE_tolerance"][-1]["relative_difference_from_baseline"] for c in controls),
        "outer_matching_order":max(c["outer_matching_order"][-1]["relative_difference_from_baseline"] for c in controls),
        "horizon_series_order":max(c["horizon_series_order"][-1]["relative_difference_from_baseline"] for c in controls)}
    maxima={"unitarity_absolute":max(r["unitarity_residual"] for r in rows),
        "outer_basis_flux_absolute":max(abs(r["outer_flux_over_omega"]-1) for r in rows),
        "horizon_flux_absolute":max(abs(r["horizon_flux_over_minus_omega"]-1) for r in rows),
        "well_conditioned_Wronskian_relative":max(r["max_relative_path_wronskian_error"] for r in high_condition),
        "scaled_Wronskian":max(r["max_scaled_path_wronskian_error"] for r in rows),
        "max_unconditioned_Wronskian_relative":max(r["max_relative_path_wronskian_error"] for r in rows),
        "well_conditioned_samples":len(high_condition),"total_scattering_samples":len(rows)}
    checks={
        "positive_transmission":bool(np.all(gamma[:,1:]>0)),
        "bounded_probabilities":bool(np.max(gamma)<=1+BOUNDS["probability_range_tolerance"] and np.max(reflection)<=1+BOUNDS["probability_range_tolerance"]),
        "flux_unitarity":maxima["unitarity_absolute"]<BOUNDS["unitarity_absolute"],
        "outer_Jost_flux_normalization":maxima["outer_basis_flux_absolute"]<BOUNDS["outer_basis_flux_absolute"],
        "horizon_Frobenius_flux_normalization":maxima["horizon_flux_absolute"]<BOUNDS["horizon_flux_absolute"],
        "Wronskian_on_transmitting_cases":maxima["well_conditioned_Wronskian_relative"]<BOUNDS["well_conditioned_Wronskian_relative"],
        "scaled_Wronskian_all_cases":maxima["scaled_Wronskian"]<BOUNDS["scaled_Wronskian"],
        "outer_radius_convergence":convergence_max["outer_radius"]<BOUNDS["outer_radius_relative"],
        "horizon_distance_convergence":convergence_max["horizon_epsilon"]<BOUNDS["horizon_epsilon_relative"],
        "ODE_tolerance_convergence":convergence_max["ODE_tolerance"]<BOUNDS["tolerance_relative"],
        "outer_series_convergence":convergence_max["outer_matching_order"]<BOUNDS["matching_order_relative"],
        "horizon_series_convergence":convergence_max["horizon_series_order"]<BOUNDS["horizon_order_relative"],
        "low_frequency_area_limit":abs(low[0]["ratio_to_16omega_squared"]-1)<BOUNDS["low_frequency_ratio_error"] and abs(low[0]["ratio_to_16omega_squared"]-1)<abs(low[-1]["ratio_to_16omega_squared"]-1),
        "frequency_grid_convergence":quadrature[-2]["relative_difference_from_finest"]<BOUNDS["power_grid_relative"],
        "independent_interpolation_integral":interpolated[-1]["relative_difference_from_Simpson_finest"]<BOUNDS["power_interpolation_relative"],
        "partial_wave_cutoff_convergence":float(powers[-1]/total)<BOUNDS["last_partial_wave_relative"],
        "retained_mode_frequency_tail_bound":tail/total<BOUNDS["retained_frequency_tail_relative"],
    }
    arrays={"frequency":frequency,"ell":np.arange(L_MAX+1),"transmission":gamma,"reflection":reflection,
        "frequency_was_numerically_scattered":frequency>0,"partial_spectrum":partial,"total_spectrum":spectrum,
        "partial_powers":powers,"cumulative_powers":np.cumsum(powers),"low_frequency":np.array([r["omega"] for r in low]),
        "low_transmission":np.array([r["transmission"] for r in low]),
        "low_ratio_to_area_limit":np.array([r["ratio_to_16omega_squared"] for r in low])}
    arrays.update(raw)
    xx=np.linspace(-20.,50.,1001)
    rr=np.array([radius(x) for x in xx])
    arrays["potential_radius"]=rr
    arrays["potential_tortoise"]=xx
    arrays["potential"]=np.array([potential(rr,l) for l in range(L_MAX+1)])
    report={"model":{"mass":1.,"field":"One real massless scalar in four-dimensional Schwarzschild; xi=0 and xi=1/6 have the same radial equation because R=0",
        "potential":"(1-2/r)[ell(ell+1)/r^2+2/r^3]","Hawking_temperature":1/(8*np.pi),
        "horizon_epsilon":HORIZON_EPS,"horizon_Frobenius_order":HORIZON_ORDER,"outer_Jost_order":OUTER_ORDER,
        "outer_radius_rule":"max(60,20/omega)","rtol":RTOL,"atol":RTOL*1e-3,
        "frequency_min":float(FREQUENCIES[0]),"frequency_max":float(FREQUENCIES[-1]),"frequency_step":FREQUENCY_STEP,
        "l_max":L_MAX,"m_squared_times_power":total,
        "scope":"Stationary fixed-background scattering and adiabatic outgoing Hawking flux; no self-consistent evaporation metric or singularity-resolution result"},
        "fixed_acceptance_bounds":BOUNDS,"flux_checks":maxima,"boundary_convergence_max_relative":convergence_max,
        "boundary_convergence_witnesses":controls,"low_frequency_area_test":low,
        "frequency_quadrature":quadrature,"independent_PCHIP_quadrature":interpolated,
        "partial_powers":powers.tolist(),"cumulative_powers":np.cumsum(powers).tolist(),
        "last_partial_wave_fraction":float(powers[-1]/total),
        "frequency_tail_upper_bound_retained_l_only":tail,
        "frequency_tail_fraction_retained_l_only":tail/total,
        "cutoff_scope":"The upper-frequency bound uses Gamma<=1 only for ell<=5. The ell cutoff is a numerical convergence study, not an analytic bound on the infinite angular sum.",
        "Wronskian_scope":"Absolute relative Wronskians are checked only where Gamma>=1e-4. At tiny transmission, cancellation of huge almost-real amplitudes destroys the unscaled flux diagnostic; the scaled residual and independent boundary/tolerance convergence are recorded instead.",
        "checks":{k:bool(v) for k,v in checks.items()},"all_checks_pass":all(checks.values())}
    return arrays,report


def plot_results(root,stem,arrays,report):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6))
    fig.subplots_adjust(left=.08,right=.98,bottom=.17,top=.82,wspace=.27)
    for ell in range(L_MAX+1):
        axes[0].plot(arrays["potential_tortoise"],arrays["potential"][ell],label=f"l={ell}")
        axes[1].plot(arrays["frequency"],arrays["transmission"][ell],label=f"l={ell}")
    axes[0].set(xlim=(-12,35),xlabel="Tortoisenkoordinate r*/m",ylabel="m² V_l",title="Vierdimensionale Regge–Wheeler-Barriere")
    axes[1].set(xlabel="mω",ylabel="Transmissionswahrscheinlichkeit Γ_l",title="Berechnete Graukörperfaktoren",ylim=(-.025,1.025))
    for ax in axes:ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.suptitle("Reelles massloses Skalarfeld auf festem Schwarzschild-Hintergrund",fontsize=13)
    fig.savefig(root/"figures"/(stem+"_barrier.png"),dpi=165);fig.savefig(root/"figures"/(stem+"_barrier.svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6))
    fig.subplots_adjust(left=.08,right=.98,bottom=.17,top=.82,wspace=.28)
    for ell in range(4):
        axes[0].plot(arrays["frequency"],arrays["partial_spectrum"][ell],label=f"l={ell}")
    axes[0].plot(arrays["frequency"],arrays["total_spectrum"],color="black",ls="--",label="Summe l=0…5")
    axes[0].set(xlim=(0,.6),xlabel="mω",ylabel="m dP/dω",title="Hawking-Spektrum mit Barrierenstreuung")
    axes[1].bar(np.arange(L_MAX+1),arrays["partial_powers"],color="#286b92")
    axes[1].set_yscale("log")
    axes[1].set(xlabel="Drehimpuls l",ylabel="m² P_l",title=f"Gesamtleistung m²P = {report['model']['m_squared_times_power']:.9g}",xticks=np.arange(L_MAX+1))
    for ax in axes:ax.grid(alpha=.22)
    axes[0].legend(fontsize=8)
    fig.suptitle("Ein reelles Skalarfeld: P=(1/2π) Σ(2l+1) ∫ ωΓ_l/(exp(8πmω)−1) dω",fontsize=12)
    fig.savefig(root/"figures"/(stem+"_spectrum.png"),dpi=165);fig.savefig(root/"figures"/(stem+"_spectrum.svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6))
    fig.subplots_adjust(left=.09,right=.98,bottom=.18,top=.82,wspace=.30)
    axes[0].loglog(arrays["low_frequency"],arrays["low_transmission"],"o-",label="Numerische Γ₀")
    axes[0].loglog(arrays["low_frequency"],16*arrays["low_frequency"]**2,"--",label="16(mω)²")
    axes[0].set(xlabel="mω",ylabel="Γ₀",title="Niederfrequenzgrenze: Horizontfläche")
    rows=report["frequency_quadrature"]
    steps=[r["frequency_step"] for r in rows]
    delta=[abs(r["total_power"]-rows[-1]["total_power"])/rows[-1]["total_power"] for r in rows]
    axes[1].loglog(steps[:-1],delta[:-1],"o-",label="Simpson gegen feinstes Gitter")
    interp=report["independent_PCHIP_quadrature"]
    axes[1].loglog(steps,[r["relative_difference_from_Simpson_finest"] for r in interp],"s-",label="PCHIP(Γ) + unabhängige Quadratur")
    axes[1].set(xlabel="Frequenzschritt Δ(mω)",ylabel="Relativer Leistungsunterschied",title="Quadratur- und Interpolationskontrolle")
    for ax in axes:ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.suptitle("Unabhängige Grenzfall- und Integrationskontrollen",fontsize=13)
    fig.savefig(root/"figures"/(stem+"_validation.png"),dpi=165);fig.savefig(root/"figures"/(stem+"_validation.svg"),metadata={"Date":None});plt.close(fig)
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--pilot",action="store_true")
    parser.add_argument("--tag",default="20261018");parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    if not args.tag.replace("_","").replace("-","").isalnum():parser.error("Invalid tag")
    if not 1<=args.workers<=8:parser.error("Choose 1 through 8 workers")
    start=time.perf_counter()
    if args.pilot:
        for ell,omega in [(0,.0002),(0,.001),(0,.01),(0,.1),(1,.1),(2,.2),(5,.01),(5,.5),(5,1.)]:
            result=scatter(omega,ell)
            if ell==0:result["ratio_to_16omega_squared"]=result["transmission"]/(16*omega*omega)
            print(json.dumps(result),flush=True)
        print("elapsed",time.perf_counter()-start,flush=True);return
    root=Path(__file__).resolve().parents[1]
    datastem="blackhole_numerics_"+args.tag;figstem="blackhole_scattering_"+args.tag
    paths=[root/"work"/(datastem+s) for s in (".json",".npz",".md")]
    paths += [root/"figures"/(figstem+s+e) for s in ("_barrier","_spectrum","_validation") for e in (".png",".svg")]
    if any(p.exists() for p in paths):parser.error("Historical outputs exist; use a new tag")
    print(f"Scattering {(L_MAX+1)*len(FREQUENCIES)} cases on {args.workers} workers...",flush=True)
    rows=scan(args.workers)
    print(f"Scattering finished after {time.perf_counter()-start:.1f}s; checking boundaries and quadratures...",flush=True)
    controls=convergence_controls()
    arrays,report=build_results(rows,controls)
    mplversion=plot_results(root,figstem,arrays,report)
    np.savez_compressed(paths[1],**arrays)
    report["provenance"]={"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":mplversion,"platform":platform.platform(),
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"arrays_sha256":hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        "command":f"python work/blackhole_numerics.py --tag {args.tag} --workers {args.workers}","elapsed_seconds":time.perf_counter()-start}
    paths[0].write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    passed=sum(report["checks"].values());total=len(report["checks"])
    reporttext=f"""# Vierdimensionale Skalarstreuung am Schwarzschild-Loch

Lauf `{args.tag}`: **{passed}/{total}** festgelegte Kontrollen bestanden.
Ein reelles massloses Skalarfeld, m=1, l=0…5, 2400 direkte Streulösungen.
Die berechnete dimensionslose Leistung beträgt **m²P={report['model']['m_squared_times_power']:.12g}**.

Die Horizontlösung wird durch eine Frobeniusreihe mit einlaufender Phase
normiert und außen an die asymptotischen Jostlösungen angepasst. Getrennte
Kontrollen variieren beide Ränder, beide Reihenordnungen und die ODE-Toleranz.
Drei verschachtelte Frequenzgitter und eine unabhängige PCHIP-Quadratur
prüfen die Flussintegration. Die Frequenznullstelle ist ein analytischer
Grenzwert und keine ausgewertete Streulösung.

Bei sehr kleiner Transmission ist ein unskalierter Wronskian aus den großen
reflektierten Amplituden durch Auslöschung schlecht konditioniert. Deshalb
wird er nur für Γ>=1e−4 als scharfer Test verwendet; alle ungekürzten Rohwerte
und zusätzlich die skalierten Residuen bleiben im Datensatz.

Der obere Frequenzrest ist nur für die behaltenen l<=5-Moden analytisch
abgeschätzt. Der Abschluss der unendlichen Drehimpulssumme wird durch
numerische Konvergenz untersucht, nicht durch einen behaupteten strikten
Restbeweis. Die Rechnung beschreibt einen festen Schwarzschild-Hintergrund
und seinen adiabatisch interpretierten Hawkingfluss; sie löst keine
selbstkonsistente Verdampfung und keine innere Raumzeitsingularität.

QuelleSHA256: `{report['provenance']['source_sha256']}`.
RoharraySHA256: `{report['provenance']['arrays_sha256']}`.
Reproduktion mit neuem Tag: `python work/blackhole_numerics.py --tag NEUER_TAG --workers 4`.
"""
    paths[2].write_text(reporttext,encoding="utf-8")
    print(json.dumps({"power":report["model"]["m_squared_times_power"],"partial_powers":report["partial_powers"],
        "flux_checks":report["flux_checks"],"boundary_convergence":report["boundary_convergence_max_relative"],
        "frequency_quadrature":report["frequency_quadrature"],"PCHIP":report["independent_PCHIP_quadrature"],
        "checks":report["checks"],"passed":passed,"total":total,"files":[str(p.relative_to(root)) for p in paths]},indent=2),flush=True)
    if not report["all_checks_pass"]:raise SystemExit(1)


if __name__=="__main__":main()
