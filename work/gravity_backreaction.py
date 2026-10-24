"""Closed finite-mode semiclassical FLRW with full Gaussian operator source.

python work/gravity_backreaction.py --tag 20261017

There are three complete cubic axis orbits, k=1,2,3. Each orbit consists of
three directions and both real sine/cosine modes, all with the same Gaussian
covariance. These are not complete Euclidean energy shells. The physical
normalization is six real modes per orbit divided by the comoving box volume.
H is independently evolved; the Friedmann square root is used only at t=0.
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
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

K = np.array([1., 2., 3.])
BOX_L = 2*np.pi
VOLUME = BOX_L**3
WEIGHTS = np.full(3,6/VOLUME)
WIDTH_FACTORS = np.array([.7,1.4,.9])
SIGMA0 = WIDTH_FACTORS/np.sqrt(2*K)
V0 = np.array([.1,-.2,.15])*SIGMA0
CGRAV = 1.0
LAMBDAS = (0.0,.03)
T_FINAL = 30.0
TIME = np.linspace(0,T_FINAL,601)
PROBE_XI = np.array([[.7,-.6,.9,-1.1,.4,-.8],[-.8,.5,1.1,-.4,.9,.3],[.6,-1.2,.8,-.5,.7,-.9]])
BOUNDS = {
    "friedmann_scaled_error":2e-8,
    "complex_background_absolute_error":2e-7,
    "complex_width_relative_error":2e-7,
    "complex_wronskian_error":2e-7,
    "conformal_background_absolute_error":2e-7,
    "conformal_width_relative_error":2e-7,
    "schrodinger_18_mode_relative":5e-5,
    "operator_continuity_scaled":2e-11,
    "continuity_finite_difference":2e-5,
    "gaussian_flow_identity":2e-10,
    "naive_constraint_integral_error":2e-7,
    "naive_defect_minimum":1e-4,
    "tolerance_convergence_fine":2e-7,
}


def source(a,sigma,v,naive=False):
    omega2=K*K/(a*a)
    qkin=1/(4*a**6*sigma**2)
    rho_B=.5*np.sum(WEIGHTS*(v*v+omega2*sigma*sigma),axis=-1)
    p_B=.5*np.sum(WEIGHTS*(v*v-omega2*sigma*sigma/3),axis=-1)
    delta=.5*np.sum(WEIGHTS*qkin,axis=-1)
    defect=np.sum(WEIGHTS*v*qkin/sigma,axis=-1)
    return (rho_B if naive else rho_B+delta),(p_B if naive else p_B+delta),rho_B,p_B,delta,defect


def initial_state(lam,naive=False):
    rho=source(1.,SIGMA0,V0,naive)[0]
    hubble=np.sqrt(CGRAV*rho+lam/3)
    return np.concatenate([[1.,hubble],SIGMA0,V0,np.zeros(3),[0.]])


def cosmic_rhs(t,y,naive=False):
    a,H=y[:2];s=y[2:5];v=y[5:8]
    rho,p,_,_,_,defect=source(a,s,v,naive)
    vd=-3*H*v-K*K*s/(a*a)+1/(4*a**6*s**3)
    return np.concatenate([[a*H,-1.5*CGRAV*(rho+p)],v,vd,-1/(4*a**3*s*s),[defect]])


def solve_cosmic(lam,rtol,naive=False):
    events=None
    if naive:
        def stop_at_H_zero(t,y):return y[1]
        stop_at_H_zero.terminal=True
        stop_at_H_zero.direction=-1
        events=stop_at_H_zero
    solution=solve_ivp(lambda t,y:cosmic_rhs(t,y,naive),(0,T_FINAL),initial_state(lam,naive),
        method="DOP853",rtol=rtol,atol=rtol*1e-3,dense_output=True,events=events)
    if not solution.success:raise RuntimeError(solution.message)
    return solution


def solve_complex(lam):
    h0=initial_state(lam)[1]
    u0=SIGMA0.astype(complex)
    du0=V0-1j/(2*SIGMA0)
    initial=np.concatenate([[1.,h0],u0,du0])
    def rhs(t,y):
        a,H=y[:2].real;u=y[2:5];du=y[5:8]
        kinetic=abs(du)**2;gradient=K*K*abs(u)**2/(a*a)
        rho=.5*np.sum(WEIGHTS*(kinetic+gradient))
        pressure=.5*np.sum(WEIGHTS*(kinetic-gradient/3))
        return np.concatenate([[a*H,-1.5*CGRAV*(rho+pressure)],du,-3*H*du-K*K*u/(a*a)])
    result=solve_ivp(rhs,(0,T_FINAL),initial,method="DOP853",rtol=2e-13,atol=2e-15,t_eval=TIME)
    if not result.success:raise RuntimeError(result.message)
    return result


def solve_conformal(lam):
    initial=initial_state(lam)
    # State a,H,sigma,w=d_eta sigma,theta,cosmic time.
    initial=np.concatenate([initial[:11],[0.]])
    def rhs(eta,y):
        a,H=y[:2];s=y[2:5];w=y[5:8]
        rho,p,*_=source(a,s,w/a)
        return np.concatenate([[a*a*H,-1.5*CGRAV*a*(rho+p)],w,
             -2*a*H*w-K*K*s+1/(4*a**4*s**3),-1/(4*a*a*s*s),[a]])
    def stop_at_final(eta,y):return y[-1]-T_FINAL
    stop_at_final.terminal=True;stop_at_final.direction=1
    result=solve_ivp(rhs,(0,T_FINAL),initial,method="DOP853",rtol=2e-13,atol=2e-15,
                     dense_output=True,events=stop_at_final)
    if not result.success or len(result.t_events[0])!=1:raise RuntimeError("Conformal integration did not reach final cosmic time")
    endpoint=result.t_events[0][0]
    etas=np.array([0. if t==0 else endpoint if t==T_FINAL else brentq(lambda e:result.sol(e)[-1]-t,0,endpoint,xtol=2e-13) for t in TIME])
    values=result.sol(etas)
    return result,etas,values


def log_wavefunction(sol,t,q):
    y=sol.sol(t);a=y[0];s=y[2:5];v=y[5:8];theta=y[8:11]
    phase_coefficient=a**3*v/s
    return np.sum(-1.5*np.log(2*np.pi*s*s)-np.sum(q*q,axis=1)/(4*s*s)
                  +1j*(.5*phase_coefficient*np.sum(q*q,axis=1)+6*theta))


def finite_difference_schrodinger(sol):
    rows=[]
    for h in (.0002,.0001,.00005):
        worst=0.0
        for t in (.25,2.,10.,25.):
            y=sol.sol(t);a=y[0];s=y[2:5]
            q=s[:,None]*PROBE_XI
            logpsi=log_wavefunction(sol,t,q);psi=np.exp(logpsi)
            derivative=(np.exp(log_wavefunction(sol,t+h,q))-np.exp(log_wavefunction(sol,t-h,q)))/(2*h)
            laplacian=0j
            for group in range(3):
                for mode in range(6):
                    delta=np.zeros_like(q);delta[group,mode]=h
                    laplacian+=(np.exp(log_wavefunction(sol,t,q+delta))-2*psi+np.exp(log_wavefunction(sol,t,q-delta)))/h**2
            potential=.5*a*np.sum(K*K*np.sum(q*q,axis=1))*psi
            residual=abs(1j*derivative+laplacian/(2*a**3)-potential)/abs(psi)
            worst=max(worst,float(residual))
        rows.append({"h":h,"times":[.25,2.,10.,25.],"real_modes":18,"max_relative_residual":worst})
    return rows


def energy_checks(sol):
    values=sol.sol(TIME)
    a,H=values[:2];s=values[2:5].T;v=values[5:8].T
    aa=a[:,None];hh=H[:,None]
    omega2=K*K/(aa*aa);qkin=1/(4*aa**6*s*s)
    vd=-3*hh*v-omega2*s+qkin/s
    qkin_dot=qkin*(-6*hh-2*v/s)
    rho,p,rhoB,pB,delta,defect=source(aa,s,v)
    rhodot=np.sum(WEIGHTS*(v*vd+omega2*s*v-hh*omega2*s*s+.5*qkin_dot),axis=1)
    balance=rhodot+3*H*(rho+p)
    delta_dot=.5*np.sum(WEIGHTS*qkin_dot,axis=1)
    correction_balance=delta_dot+6*H*delta
    mean_balance=rhodot-delta_dot+3*H*(rhoB+pB)
    rows=[]
    for h in (.002,.001,.0005):
        worst=0.
        for t in (.25,2.,10.,25.):
            def rho_at(tau):
                y=sol.sol(tau);return source(y[0],y[2:5],y[5:8])[0]
            y=sol.sol(t);r,pp,*_=source(y[0],y[2:5],y[5:8])
            drho=(rho_at(t+h)-rho_at(t-h))/(2*h)
            worst=max(worst,abs(drho+3*y[1]*(r+pp)))
        rows.append({"h":h,"max_absolute_residual":float(worst)})
    return {"max_operator_continuity_scaled":float(np.max(abs(balance))/(1+np.max(abs(rhodot)))),
            "max_bohm_defect_identity":float(np.max(abs(mean_balance-defect))),
            "max_operator_compensation_error":float(np.max(abs(defect+correction_balance))),
            "max_naive_ensemble_source_defect":float(np.max(abs(defect))),"finite_difference":rows},(rho,p,rhoB,pB,delta,defect,balance)


def inspect_case(lam,arrays):
    label="L0" if lam==0 else "L003"
    reference=solve_cosmic(lam,2e-13)
    values=reference.sol(TIME)
    a,H=values[:2];s=values[2:5].T;v=values[5:8].T;theta=values[8:11].T
    energy,(rho,p,rhoB,pB,delta,defect,balance)=energy_checks(reference)
    constraint=H*H-CGRAV*rho-lam/3
    convergence=[]
    for tol in (1e-7,1e-9,1e-11):
        candidate=solve_cosmic(lam,tol)
        yy=candidate.sol(TIME)
        crho=source(yy[0,:,None],yy[2:5].T,yy[5:8].T)[0]
        drift=yy[1]**2-CGRAV*crho-lam/3
        scaled_error=np.max(abs(yy[:8]-values[:8])/(1+abs(values[:8])))
        convergence.append({"rtol":tol,"atol":tol*1e-3,"nfev":candidate.nfev,
            "state_scaled_error_vs_rtol2e13":float(scaled_error),"constraint_max_absolute":float(np.max(abs(drift)))})
    complex_result=solve_complex(lam)
    cu=complex_result.y[2:5].T;cdu=complex_result.y[5:8].T
    ca=complex_result.y[0].real;cH=complex_result.y[1].real
    wronskian=ca[:,None]**3*(cu*cdu.conj()-cu.conj()*cdu)
    complex_info={"nfev":complex_result.nfev,"max_background_absolute_error":float(max(np.max(abs(ca-a)),np.max(abs(cH-H)))),
        "max_width_relative_error":float(np.max(abs(abs(cu)-s)/s)),
        "max_wronskian_error":float(np.max(abs(wronskian-1j))),
        "max_global_phase_error":float(np.max(abs(.5*np.unwrap(np.angle(cu),axis=0)-theta))),
        "method":"Independent coupled complex KG modes plus a,H; rho and p computed directly from |u_dot|^2 and |u|^2"}
    conformal_result,etas,eta_y=solve_conformal(lam)
    conformal_info={"nfev":conformal_result.nfev,"eta_final":float(etas[-1]),
        "max_background_absolute_error":float(max(np.max(abs(eta_y[0]-a)),np.max(abs(eta_y[1]-H)))),
        "max_width_relative_error":float(np.max(abs(eta_y[2:5].T-s)/s)),
        "max_clock_absolute_error":float(np.max(abs(eta_y[-1]-TIME))),
        "method":"Separate conformal-time ODE with w=d_eta sigma, independent cosmic clock and dense-output root inversion"}
    pde=finite_difference_schrodinger(reference)
    logflow_error=0.0
    for ss in s:
        q=ss[:,None]*PROBE_XI
        logp=np.sum(-3*np.log(2*np.pi*ss*ss)-np.sum(q*q,axis=1)/(2*ss*ss))
        q0=SIGMA0[:,None]*PROBE_XI
        logp0=np.sum(-3*np.log(2*np.pi*SIGMA0*SIGMA0)-np.sum(q0*q0,axis=1)/(2*SIGMA0*SIGMA0))
        logjac=6*np.sum(np.log(ss/SIGMA0))
        logflow_error=max(logflow_error,abs(logp+logjac-logp0))
    naive=solve_cosmic(lam,2e-13,True)
    naive_end=float(naive.t[-1]);nt=np.linspace(0,naive_end,601);ny=naive.sol(nt)
    nrho,npres,nrhoB,npB,ndelta,ndefect=source(ny[0,:,None],ny[2:5].T,ny[5:8].T,True)
    nconstraint=ny[1]**2-CGRAV*nrho-lam/3
    consistency=nconstraint+CGRAV*ny[-1]
    naive_info={"time_end":naive_end,"termination":"first H=0" if naive_end<T_FINAL-1e-8 else "fixed t=30",
        "initial_H":float(ny[1,0]),"final_H":float(ny[1,-1]),"initial_rho":float(nrho[0]),
        "max_constraint_absolute":float(np.max(abs(nconstraint))),"final_constraint":float(nconstraint[-1]),
        "constraint_plus_integrated_defect_max_error":float(np.max(abs(consistency))),
        "minimum_sigma":float(np.min(ny[2:5])),
        "interpretation":"Deliberately inconsistent closure using only ensemble means of classical Bohm-path tensors; independently constrained initial H, same a/sigma/v. No physical bounce is inferred."}
    arrays.update({f"a_{label}":a,f"H_{label}":H,f"sigma_{label}":s,f"v_{label}":v,f"theta_{label}":theta,
        f"rho_{label}":rho,f"p_{label}":p,f"rho_B_{label}":rhoB,f"p_B_{label}":pB,
        f"delta_{label}":delta,f"defect_B_{label}":defect,f"constraint_{label}":constraint,
        f"operator_balance_{label}":balance,f"integrated_defect_B_{label}":values[-1],
        f"complex_u_{label}":cu,f"complex_udot_{label}":cdu,f"complex_a_{label}":ca,f"complex_H_{label}":cH,
        f"complex_wronskian_{label}":wronskian,f"conformal_eta_{label}":etas,
        f"conformal_a_{label}":eta_y[0],f"conformal_H_{label}":eta_y[1],
        f"conformal_sigma_{label}":eta_y[2:5].T,f"conformal_w_{label}":eta_y[5:8].T,
        f"naive_time_{label}":nt,f"naive_a_{label}":ny[0],f"naive_H_{label}":ny[1],
        f"naive_sigma_{label}":ny[2:5].T,f"naive_v_{label}":ny[5:8].T,
        f"naive_constraint_{label}":nconstraint,f"naive_integrated_defect_{label}":ny[-1]})
    result={"lambda":lam,"initial_H":float(H[0]),"initial_rho":float(rho[0]),"initial_p":float(p[0]),
        "final_a":float(a[-1]),"final_H":float(H[-1]),"final_rho":float(rho[-1]),"final_p":float(p[-1]),
        "minimum_sigma":float(s.min()),"minimum_H":float(H.min()),"max_rho":float(rho.max()),
        "w_range":[float(np.min(p/rho)),float(np.max(p/rho))],
        "constraint_max_absolute":float(np.max(abs(constraint))),
        "constraint_max_scaled":float(np.max(abs(constraint)/(1+CGRAV*rho+lam/3))),
        "tolerance_convergence":convergence,"energy":energy,"independent_complex_KG":complex_info,
        "conformal_time_comparison":conformal_info,"schrodinger_finite_difference_18_modes":pde,
        "log_gaussian_density_jacobian_error":float(logflow_error),"naive_control":naive_info}
    bb=BOUNDS
    checks={"friedmann_constraint":result["constraint_max_scaled"]<bb["friedmann_scaled_error"],
        "positive_widths":s.min()>0,"expanding_background":H.min()>0,
        "operator_w_bounds":np.min(p/rho)>=-1/3-1e-12 and np.max(p/rho)<=1+1e-12,
        "tolerance_convergence":convergence[-1]["state_scaled_error_vs_rtol2e13"]<bb["tolerance_convergence_fine"] and convergence[-1]["state_scaled_error_vs_rtol2e13"]<convergence[0]["state_scaled_error_vs_rtol2e13"]/100,
        "independent_complex_background":complex_info["max_background_absolute_error"]<bb["complex_background_absolute_error"],
        "independent_complex_width":complex_info["max_width_relative_error"]<bb["complex_width_relative_error"],
        "canonical_wronskian":complex_info["max_wronskian_error"]<bb["complex_wronskian_error"],
        "conformal_background":conformal_info["max_background_absolute_error"]<bb["conformal_background_absolute_error"],
        "conformal_width":conformal_info["max_width_relative_error"]<bb["conformal_width_relative_error"],
        "schrodinger_18_modes":pde[-1]["max_relative_residual"]<bb["schrodinger_18_mode_relative"],
        "operator_source_continuity":energy["max_operator_continuity_scaled"]<bb["operator_continuity_scaled"],
        "source_continuity_finite_difference":energy["finite_difference"][-1]["max_absolute_residual"]<bb["continuity_finite_difference"],
        "density_jacobian":logflow_error<bb["gaussian_flow_identity"],
        "naive_constraint_drift_detected":naive_info["max_constraint_absolute"]>bb["naive_defect_minimum"],
        "naive_drift_integral_identity":naive_info["constraint_plus_integrated_defect_max_error"]<bb["naive_constraint_integral_error"]}
    result["acceptance_checks"]={k:bool(v) for k,v in checks.items()};result["all_checks_pass"]=all(checks.values())
    print(json.dumps({"lambda":lam,"initial_H":result["initial_H"],"final_a":result["final_a"],
        "constraint_max":result["constraint_max_absolute"],"complex":complex_info,"conformal":conformal_info,
        "pde":pde,"naive":naive_info,"checks":result["acceptance_checks"]},indent=2),flush=True)
    return result


def plot_results(root,stem,arrays,metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8.5))
    fig.subplots_adjust(left=.09,right=.97,bottom=.08,top=.89,wspace=.29,hspace=.35)
    colors=("#17628d","#b64832")
    for case,color,label in zip(metrics["cases"],colors,("L0","L003")):
        lam=case["lambda"];legend=f"Λ={lam:g}"
        axes[0,0].plot(TIME,arrays[f"a_{label}"],color=color,label=legend)
        axes[0,1].semilogy(TIME,arrays[f"rho_{label}"],color=color,label=legend)
        axes[1,0].plot(TIME,arrays[f"p_{label}"]/arrays[f"rho_{label}"],color=color,label=legend)
        axes[1,1].plot(arrays[f"naive_time_{label}"],arrays[f"naive_constraint_{label}"],color=color,ls="--",label=f"Naiver Abschluss, {legend}")
        axes[1,1].plot(TIME,arrays[f"constraint_{label}"],color=color,label=f"Operatorquelle, {legend}")
    axes[0,0].set(xlabel="Kosmische Zeit t",ylabel="Skalenfaktor a",title="Selbstkonsistent berechnete Expansion")
    axes[0,1].set(xlabel="Kosmische Zeit t",ylabel="Feldenergiedichte ρ",title="Vollständige endliche Operatorquelle")
    axes[1,0].axhline(-1/3,color="gray",ls=":",lw=1);axes[1,0].axhline(1,color="gray",ls=":",lw=1)
    axes[1,0].set(xlabel="Kosmische Zeit t",ylabel="w = p / ρ",title="Zustandsgleichung des Feldzustands")
    axes[1,1].set(xlabel="Kosmische Zeit t",ylabel="C = H² − cρ − Λ/3",title="Constraint: Operatorquelle und naiver Vergleich")
    axes[1,1].set_yscale("symlog",linthresh=1e-10)
    for ax in axes.flat:ax.grid(alpha=.23);ax.legend(fontsize=8)
    fig.suptitle("Geschlossenes semiklassisches FLRW mit 18 reellen Gaußmoden\n"
                 "Drei kubische Achsenorbits, k=1,2,3; Volumen (2π)³, c=8πG/3=1",fontsize=13)
    fig.savefig(root/"figures"/(stem+".png"),dpi=165)
    fig.savefig(root/"figures"/(stem+".svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    fig.subplots_adjust(left=.08,right=.98,bottom=.17,top=.84,wspace=.3)
    for case,color in zip(metrics["cases"],colors):
        rows=case["tolerance_convergence"]
        axes[0].loglog([r["rtol"] for r in rows],[r["state_scaled_error_vs_rtol2e13"] for r in rows],"o-",color=color,label=f"Λ={case['lambda']:g}")
        rows=case["schrodinger_finite_difference_18_modes"]
        axes[1].loglog([r["h"] for r in rows],[r["max_relative_residual"] for r in rows],"o-",color=color,label=f"Λ={case['lambda']:g}")
    axes[0].set(xlabel="DOP853-Relativtoleranz",ylabel="Maximaler skalierter Zustandsfehler",title="Konvergenz der gekoppelten Entwicklung")
    axes[1].set(xlabel="Differenzenschritt h",ylabel="Relatives Schrödingerresiduum",title="Vollständige Wellenfunktion: 18 Moden")
    for ax in axes:ax.grid(alpha=.23);ax.legend(fontsize=9)
    fig.savefig(root/"figures"/(stem+"_validation.png"),dpi=165)
    fig.savefig(root/"figures"/(stem+"_validation.svg"),metadata={"Date":None});plt.close(fig)
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--tag",default="20261017")
    parser.add_argument("--pilot",action="store_true");args=parser.parse_args()
    if not args.tag.replace("_","").replace("-","").isalnum():parser.error("Invalid tag")
    root=Path(__file__).resolve().parents[1];stem="gravity_backreaction_"+args.tag
    paths=[root/"work"/(stem+e) for e in (".json",".npz",".md")]+[root/"figures"/(stem+e) for e in (".png",".svg","_validation.png","_validation.svg")]
    if not args.pilot and any(p.exists() for p in paths):parser.error("Historical outputs exist; use a new tag")
    start=time.perf_counter();arrays={"time":TIME,"k":K,"weights":WEIGHTS,"initial_sigma":SIGMA0,"initial_v":V0,"probe_xi":PROBE_XI}
    metrics={"model":{"box_length":float(BOX_L),"comoving_volume":float(VOLUME),"real_modes":18,
        "groups":"three complete cubic axis orbits, each three directions times real sine/cosine pair",
        "k":K.tolist(),"weights":WEIGHTS.tolist(),"sigma0_width_factors":WIDTH_FACTORS.tolist(),
        "sigma0":SIGMA0.tolist(),"v0":V0.tolist(),"a0":1.,"c_gravity":CGRAV,"lambdas":list(LAMBDAS),
        "time_interval":[0,T_FINAL],"H_evolution":"Raychaudhuri ODE; square root only for initial constraint",
        "regularization":"fixed finite comoving mode set; no infinite renormalized stress tensor or adiabatic/BD assumption"},
        "fixed_acceptance_bounds":BOUNDS,"cases":[inspect_case(lam,arrays) for lam in LAMBDAS]}
    metrics["all_checks_pass"]=all(case["all_checks_pass"] for case in metrics["cases"])
    if args.pilot:
        print(json.dumps({"all_checks_pass":metrics["all_checks_pass"],"case_summaries":[{k:c[k] for k in ("lambda","initial_rho","initial_p","final_H","final_rho","w_range","energy")} for c in metrics["cases"]]},indent=2),flush=True);return
    version=plot_results(root,stem,arrays,metrics);np.savez_compressed(paths[1],**arrays)
    metrics["provenance"]={"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":version,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"arrays_sha256":hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        "command":f"python work/gravity_backreaction.py --tag {args.tag}","elapsed_seconds":time.perf_counter()-start}
    paths[0].write_text(json.dumps(metrics,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    total=sum(len(c["acceptance_checks"]) for c in metrics["cases"]);passed=sum(sum(c["acceptance_checks"].values()) for c in metrics["cases"])
    report=f"""# Geschlossenes semiklassisches Gauß-FLRW-Modell

Verbindlicher Lauf: **{args.tag}**, **{passed}/{total}** feste Prüfungen bestanden.
Reproduktion mit neuem Tag: `python work/gravity_backreaction.py --tag NEUER_TAG`.

Die Metrik wird gemeinsam mit den Feldbreiten berechnet. Der Hubbleparameter
wird durch die Raychaudhuri-Gleichung entwickelt und nicht nach jedem Schritt
auf eine Friedmann-Wurzel projiziert. Nur H(0)>0 erfüllt den Anfangsconstraint.
Die vollständige Operatorquelle enthält die Impulsvarianz 1/(4a⁶σ²).

Das periodische Volumen ist (2π)³. Für k=1,2,3 enthält jede kubische Achsenbahn
drei Richtungen mit je einer reellen Sinus- und Kosinusmode. Alle sechs Moden
einer Gruppe haben gleiche unabhängige Gaußkovarianzen. Das Gewicht ist 6/V.
Diese Modenauswahl erzeugt einen homogenen, isotropen Erwartungs-Stresstensor;
sie ist keine vollständige euklidische Energieschale und keine volle Rotations-
invarianz sämtlicher höherer Feldkorrelationen. Es handelt sich um ein endliches
reguliertes System, ohne UV-renormierte unendliche Quelle und ohne angenommene
adiabatische oder Bunch-Davies-Zustandspräparation.

Unabhängige Kontrollen umfassen eine direkt gekoppelte komplexe KG-Formulierung
mit dem kanonischen Wronskian, eine getrennte Konformzeitintegration mit eigener
Uhr und die vollständige 18-Moden-Schrödingergleichung durch finite Differenzen
einschließlich integrierter globaler Phase. Die Gauß-Äquivarianz wird mittels
Dichte-Jakobian geprüft. Drei Integrationstoleranzen prüfen die Konvergenz.

Der naive Vergleich lässt die Impulsvarianzquelle weg, behält jedoch dieselbe
Quantenführung. Er beginnt mit denselben a,σ,v und seinem eigenen anfänglich
erfüllten Friedmannconstraint. Danach erzeugt die fehlende Kontinuität die
explizit geprüfte Drift C_B(t)=−c∫D_B dt. Falls H erstmals null wird, endet
dieser Kontrolllauf dort. Daraus wird kein physikalischer Bounce abgeleitet.

QuelleSHA256: `{metrics['provenance']['source_sha256']}`.
RoharraySHA256: `{metrics['provenance']['arrays_sha256']}`.
Alle Zahlen, Zeitreihen und Schwellen: `{paths[0].name}`.
"""
    paths[2].write_text(report,encoding="utf-8")
    print(json.dumps({"all_checks_pass":metrics["all_checks_pass"],"passed":passed,"total":total,
                     "files":[str(p.relative_to(root)) for p in paths]},indent=2),flush=True)
    if not metrics["all_checks_pass"]:raise SystemExit(1)


if __name__=="__main__":main()
