"""Renormalized conformal scalar FLRW: two b=0 branches and a b!=0 bounce.

Reproduce with: python work/gravity_renormalized.py --tag NEW_TAG
All quantities below are dimensionless: tau=t/sqrt(alpha), h=sqrt(alpha) H,
alpha=c_G A, A=1/(960*pi**2), ell=alpha*Lambda/3, beta=b/A.
The effective local R^2 coupling beta is part of the model. Changing a
subtraction convention while compensating gravitational couplings does not
change the physical solution; the examples here change the effective coupling.
No ultraviolet mode cutoff, experimental fit, or generic quantum-gravity claim.
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
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

A_ANOMALY = 1/(960*np.pi**2)
D_RADIATION = .04
ELLS = (0., .02)
FUTURE_TIME = np.linspace(0.,30.,301)
DELTAS = np.logspace(-8,-2,25)
BOUNCE_TIME = np.linspace(-6.,6.,601)
BOUNDS = {
    "branch_constraint_absolute":2e-8,
    "branch_ODE_quadrature_scaled":2e-7,
    "branch_clock_absolute":2e-8,
    "past_ODE_quadrature_scaled":2e-6,
    "tolerance_fine_scaled":2e-7,
    "source_continuity_FD":2e-6,
    "trace_einstein_absolute":2e-9,
    "asymptotic_hprime_coefficient_relative":2e-6,
    "asymptotic_curvature_coefficient_relative":2e-3,
    "asymptotic_log_slope_absolute":.015,
    "deSitter_null_integral_absolute":2e-11,
    "bounce_exact_scaled":2e-7,
    "bounce_Radau_scaled":2e-7,
    "bounce_constraint_absolute":2e-7,
    "bounce_FD_balance_absolute":2e-6,
    "bounce_FD_trace_absolute":2e-5,
    "bounce_source_cancellation":2e-14,
}


def branch_data(loga,ell,branch):
    """Algebraic branch; branch=-1 is Einstein/low branch, +1 high branch."""
    rad=D_RADIATION*np.exp(-4*np.asarray(loga))
    discr=1-4*(ell+rad)
    if np.any(discr < -2e-14):
        raise ValueError("No real b=0 branch below its limiting scale")
    z=np.sqrt(np.maximum(0.,discr))
    # Rationalize the low branch to avoid cancellation in its vacuum limit.
    y=2*(ell+rad)/(1+z) if branch==-1 else (1+z)/2
    h=np.sqrt(y)
    hp=branch*2*rad/z
    return h,hp,rad,z


def branch_source(h,hp,rad):
    rho=h**4+rad
    pressure=-h**4-(4/3)*h*h*hp+rad/3
    trace=-rho+3*pressure
    ricci=6*(hp+2*h*h)
    kretschmann=12*((hp+h*h)**2+h**4)
    return rho,pressure,trace,ricci,kretschmann


def branch_clock(loga,ell,branch):
    return quad(lambda x:1/float(branch_data(x,ell,branch)[0]),0.,float(loga),
                epsabs=2e-12,epsrel=2e-13,limit=150)[0]


def branch_logscale(tau,ell,branch):
    if tau==0:return 0.
    return brentq(lambda x:branch_clock(x,ell,branch)-tau,0.,tau+1.,xtol=2e-13)


def boundary_clock(z,ell,branch):
    """Elapsed time from branch boundary; regular integrand in z=sqrt(discriminant)."""
    B=1-4*ell
    return quad(lambda zz:zz/(2*np.sqrt((1+branch*zz)/2)*(B-zz*zz)),
                0.,float(z),epsabs=2e-14,epsrel=2e-13,limit=150)[0]


def branch_rhs(tau,y):
    loga,h=y
    return np.array([h,-2*D_RADIATION*np.exp(-4*loga)/(1-2*h*h)])


def solve_branch(ell,branch,tol,past=False):
    h0=float(branch_data(0.,ell,branch)[0])
    initial=[0.,h0]
    events=None
    endpoint=30.
    if past:
        B=1-4*ell; astar=(4*D_RADIATION/B)**.25
        tau_star=-boundary_clock(np.sqrt(B-4*D_RADIATION),ell,branch)
        target=np.log(astar)+np.log1p(1e-6)
        def stop(t,y):return y[0]-target
        stop.terminal=True;stop.direction=-1
        events=stop;endpoint=1.01*tau_star
    sol=solve_ivp(branch_rhs,(0.,endpoint),initial,method="DOP853",rtol=tol,
                  atol=tol*1e-3,dense_output=True,events=events)
    if not sol.success:raise RuntimeError(sol.message)
    if past and len(sol.t_events[0])!=1:raise RuntimeError("Past integration failed to reach its nonsingular stop")
    return sol


def branch_FD_balance(ell,branch):
    rows=[]
    for step in (1e-3,5e-4,2.5e-4):
        worst=0.
        for tau in (.25,1.,3.,10.,25.):
            def at(t):
                uu=branch_logscale(t,ell,branch)
                hh,hp,rad,_=branch_data(uu,ell,branch)
                rho,p,*_=branch_source(hh,hp,rad)
                return rho,p,hh
            rho,p,hh=at(tau)
            drho=(at(tau+step)[0]-at(tau-step)[0])/(2*step)
            worst=max(worst,abs(drho+3*hh*(rho+p)))
        rows.append({"step":step,"max_absolute_residual":float(worst)})
    return rows


def inspect_branch(ell,branch,arrays):
    label=("ell0" if ell==0 else "ell002")+("_low" if branch==-1 else "_high")
    reference=solve_branch(ell,branch,2e-13)
    yy=reference.sol(FUTURE_TIME);loga,h=yy
    rad=D_RADIATION*np.exp(-4*loga)
    hp=-2*rad/(1-2*h*h)
    rho,p,trace,ricci,kretschmann=branch_source(h,hp,rad)
    constraint=h*h-rho-ell
    qloga=np.array([branch_logscale(t,ell,branch) for t in FUTURE_TIME])
    qh=np.array(branch_data(qloga,ell,branch)[0])
    qexpected=np.vstack((qloga,qh))
    qerr=float(np.max(abs(yy-qexpected)/(1+abs(qexpected))))
    clock_error=float(max(abs(branch_clock(u,ell,branch)-t) for u,t in zip(loga[::15],FUTURE_TIME[::15])))
    convergence=[]
    for tol in (1e-7,1e-9,1e-11):
        sol=solve_branch(ell,branch,tol);val=sol.sol(FUTURE_TIME)
        err=float(np.max(abs(val-qexpected)/(1+abs(qexpected))))
        cc=val[1]**2-val[1]**4-ell-D_RADIATION*np.exp(-4*val[0])
        convergence.append({"rtol":tol,"atol":tol*1e-3,"nfev":sol.nfev,
                            "scaled_error_vs_independent_quadrature":err,"constraint_max":float(np.max(abs(cc)))})
    B=1-4*ell;astar=(4*D_RADIATION/B)**.25
    tau_star=-boundary_clock(np.sqrt(B-4*D_RADIATION),ell,branch)
    nearz=np.sqrt(B*(-np.expm1(-4*np.log1p(DELTAS))))
    nearh=np.sqrt((1+branch*nearz)/2)
    nearhp=branch*(B-nearz*nearz)/(2*nearz)
    nearrad=(B-nearz*nearz)/4
    near_rho,near_p,near_trace,near_R,near_K=branch_source(nearh,nearhp,nearrad)
    elapsed=np.array([boundary_clock(z,ell,branch) for z in nearz])
    slope_hp=float(np.polyfit(np.log(DELTAS[:9]),np.log(abs(nearhp[:9])),1)[0])
    slope_K=float(np.polyfit(np.log(elapsed[:9]),np.log(near_K[:9]),1)[0])
    hp_coeff=float(abs(abs(nearhp[0])*np.sqrt(DELTAS[0])/(np.sqrt(B)/4)-1))
    k_coeff=float(abs(near_K[0]*DELTAS[0]/(.75*B)-1))
    past=solve_branch(ell,branch,2e-13,past=True)
    pdelta=np.geomspace(1e-5,1/astar-1,101)
    pz=np.sqrt(B*(-np.expm1(-4*np.log1p(pdelta))))
    ptime=tau_star+np.array([boundary_clock(z,ell,branch) for z in pz])
    ptime[-1]=0.
    pexpected=np.vstack((np.log(astar)+np.log1p(pdelta),np.sqrt((1+branch*pz)/2)))
    py=past.sol(ptime)
    past_error=float(np.max(abs(py-pexpected)/(1+abs(pexpected))))
    pcon=py[1]**2-py[1]**4-ell-D_RADIATION*np.exp(-4*py[0])
    fd=branch_FD_balance(ell,branch)
    trace_einstein_error=float(np.max(abs(ricci-(12*ell-3*trace))))
    arrays.update({f"loga_{label}":loga,f"a_{label}":np.exp(loga),f"h_{label}":h,
        f"hprime_{label}":hp,f"rho_{label}":rho,f"pressure_{label}":p,
        f"constraint_{label}":constraint,f"Ricci_{label}":ricci,f"Kretschmann_{label}":kretschmann,
        f"quadrature_loga_{label}":qloga,f"quadrature_h_{label}":qh,
        f"past_time_{label}":ptime,f"past_loga_{label}":py[0],f"past_h_{label}":py[1],
        f"past_quadrature_loga_{label}":pexpected[0],f"past_quadrature_h_{label}":pexpected[1],
        f"near_time_from_boundary_{label}":elapsed,f"near_h_{label}":nearh,
        f"near_hprime_{label}":nearhp,f"near_rho_{label}":near_rho,
        f"near_pressure_{label}":near_p,f"near_Ricci_{label}":near_R,f"near_Kretschmann_{label}":near_K})
    result={"ell":ell,"branch":branch,"label":label,"D":D_RADIATION,
        "a_boundary":float(astar),"h_boundary":float(1/np.sqrt(2)),"tau_boundary":tau_star,
        "h_initial":float(h[0]),"a_final":float(np.exp(loga[-1])),"h_final":float(h[-1]),
        "deSitter_h_limit":float(np.sqrt(2*ell/(1+np.sqrt(1-4*ell))) if branch==-1 else np.sqrt((1+np.sqrt(1-4*ell))/2)),
        "future_constraint_max_absolute":float(np.max(abs(constraint))),
        "quadrature_comparison_scaled_error":qerr,"quadrature_clock_error":clock_error,
        "past_comparison_scaled_error":past_error,"past_constraint_max_absolute":float(np.max(abs(pcon))),
        "past_stop_delta":1e-6,"near_delta_smallest":float(DELTAS[0]),
        "near_hprime_log_slope":slope_hp,"near_K_vs_elapsed_log_slope":slope_K,
        "near_hprime_coefficient_relative_error":hp_coeff,"near_K_coefficient_relative_error":k_coeff,
        "near_K_at_smallest_delta":float(near_K[0]),"near_Ricci_at_smallest_delta":float(near_R[0]),
        "near_pressure_at_smallest_delta":float(near_p[0]),
        "trace_Einstein_max_error":trace_einstein_error,"source_continuity_finite_difference":fd,
        "tolerance_convergence":convergence}
    checks={"future_constraint":result["future_constraint_max_absolute"]<BOUNDS["branch_constraint_absolute"],
        "ODE_vs_quadrature":qerr<BOUNDS["branch_ODE_quadrature_scaled"],
        "quadrature_clock":clock_error<BOUNDS["branch_clock_absolute"],
        "past_ODE_vs_quadrature":past_error<BOUNDS["past_ODE_quadrature_scaled"],
        "past_constraint":result["past_constraint_max_absolute"]<BOUNDS["branch_constraint_absolute"],
        "tolerance_convergence":convergence[-1]["scaled_error_vs_independent_quadrature"]<BOUNDS["tolerance_fine_scaled"] and convergence[-1]["scaled_error_vs_independent_quadrature"]<convergence[0]["scaled_error_vs_independent_quadrature"]/30,
        "source_continuity_FD":fd[-1]["max_absolute_residual"]<BOUNDS["source_continuity_FD"],
        "trace_Einstein":trace_einstein_error<BOUNDS["trace_einstein_absolute"],
        "asymptotic_hprime_coefficient":hp_coeff<BOUNDS["asymptotic_hprime_coefficient_relative"],
        "asymptotic_K_coefficient":k_coeff<BOUNDS["asymptotic_curvature_coefficient_relative"],
        "asymptotic_hprime_exponent":abs(slope_hp+.5)<BOUNDS["asymptotic_log_slope_absolute"],
        "asymptotic_K_exponent":abs(slope_K+1)<BOUNDS["asymptotic_log_slope_absolute"]}
    result["checks"]={k:bool(v) for k,v in checks.items()}
    print(json.dumps({k:result[k] for k in ("label","tau_boundary","a_boundary","quadrature_comparison_scaled_error","past_comparison_scaled_error","near_hprime_log_slope","near_K_vs_elapsed_log_slope","checks")},indent=2),flush=True)
    return result


def inspect_deSitter():
    rows=[]
    for ell in ELLS:
        for branch in (-1,1):
            y=2*ell/(1+np.sqrt(1-4*ell)) if branch==-1 else (1+np.sqrt(1-4*ell))/2
            h=np.sqrt(y)
            item={"ell":ell,"branch":branch,"h":float(h),"D":0.,
                  "constraint_absolute":float(abs(h*h-h**4-ell)),
                  "Ricci":float(12*h*h),"Kretschmann":float(24*h**4)}
            if h>0:
                exact=1/h
                numeric=quad(lambda t:np.exp(h*t),-100/h,0.,epsabs=2e-13,epsrel=2e-13)[0]
                item.update({"past_null_affine_length_a0_one":float(exact),
                             "numeric_integral_from_minus100_over_h":float(numeric),
                             "finite_integral_error":float(abs(numeric-exact)),
                             "comoving_proper_past_time":"infinite",
                             "geometric_statement":"Past null incompleteness of the flat expanding de Sitter patch; bounded curvature. No general extension claim."})
            else:
                item.update({"past_null_affine_length_a0_one":"infinite","finite_integral_error":0.,
                             "comoving_proper_past_time":"infinite","geometric_statement":"Minkowski lower branch; not de Sitter and no past null incompleteness."})
            rows.append(item)
    return {"cases":rows,"checks":{"constant_constraint":all(r["constraint_absolute"]<2e-14 for r in rows),
                "null_affine_quadrature":all(r["finite_integral_error"]<BOUNDS["deSitter_null_integral_absolute"] for r in rows)}}


def bounce_parameters(u):
    return (2*u-1)/(6*u),u-u*u,u*(4*u-3)/2


def bounce_exact(t,u=1.):
    tt=np.sqrt(u)*np.asarray(t);a=np.cosh(tt);h=np.sqrt(u)*np.tanh(tt);z=1/a**2
    return np.array([a,h,6*u*(2-z),12*u*z*h])


def bounce_rhs(t,y,u=1.):
    beta,ell,D=bounce_parameters(u)
    a,h,r,s=y
    hp=r/6-2*h*h
    sp=(r-12*ell+12*h**4-2*h*h*r)/(3*beta)-3*h*s
    return np.array([a*h,hp,s,sp])


def solve_bounce(tol,method="DOP853",u=1.):
    initial=[1.,0.,6*u,0.]
    left=solve_ivp(lambda t,y:bounce_rhs(t,y,u),(0.,-6.),initial,method=method,rtol=tol,atol=tol*1e-3,dense_output=True)
    right=solve_ivp(lambda t,y:bounce_rhs(t,y,u),(0.,6.),initial,method=method,rtol=tol,atol=tol*1e-3,dense_output=True)
    if not left.success or not right.success:raise RuntimeError("Trace ODE bounce integration failed")
    def sample(t):
        xx=np.atleast_1d(t);val=np.empty((4,len(xx)));mask=xx<0
        if np.any(mask):val[:,mask]=left.sol(xx[mask])
        if np.any(~mask):val[:,~mask]=right.sol(xx[~mask])
        return val[:,0] if np.ndim(t)==0 else val
    return sample,left.nfev+right.nfev


def bounce_source(y,u=1.):
    beta,ell,D=bounce_parameters(u)
    a,h,r,s=y;hp=r/6-2*h*h
    sp=(r-12*ell+12*h**4-2*h*h*r)/(3*beta)-3*h*s
    curvature=h**4+beta*(h*s+h*h*r-r*r/12)
    radiation=D/a**4
    rho=curvature+radiation
    trace=-4*h*h*(hp+h*h)-beta*(sp+3*h*s)
    pressure=(rho+trace)/3
    constraint=h*h-rho-ell
    kretschmann=12*((hp+h*h)**2+h**4)
    return rho,pressure,trace,constraint,curvature,radiation,kretschmann


def inspect_bounce(arrays,u=1.):
    beta,ell,D=bounce_parameters(u)
    prefix="bounce" if u==1 else "vacuum_bounce"
    exact=bounce_exact(BOUNCE_TIME,u)
    sample,nfev=solve_bounce(2e-13,u=u)
    y=sample(BOUNCE_TIME)
    rho,p,trace,constraint,curv,rad,kretschmann=bounce_source(y,u)
    error=float(np.max(abs(y-exact)/(1+abs(exact))))
    radau,nrad=solve_bounce(2e-11,"Radau",u=u)
    ry=radau(BOUNCE_TIME)
    radau_error=float(np.max(abs(ry-exact)/(1+abs(exact))))
    conv=[]
    for tol in (1e-7,1e-9,1e-11):
        solver,n=solve_bounce(tol,u=u);val=solver(BOUNCE_TIME)
        conv.append({"rtol":tol,"nfev":n,"scaled_error_vs_exact":float(np.max(abs(val-exact)/(1+abs(exact))))})
    fd=[]
    for dt in (1e-3,5e-4,2.5e-4):
        continuity=0.;trace_error=0.;constraint_identity=0.
        for t in (-4.,-1.,0.,1.,4.):
            yp=bounce_exact(t+dt,u);ym=bounce_exact(t-dt,u);yc=bounce_exact(t,u)
            rc,pc,tc,cc,*_=bounce_source(yc,u)
            rp,_,_,cp,*_=bounce_source(yp,u);rm,_,_,cm,*_=bounce_source(ym,u)
            continuity=max(continuity,abs((rp-rm)/(2*dt)+3*yc[1]*(rc+pc)))
            rsecond=(yp[2]-2*yc[2]+ym[2])/dt**2
            rfirst=(yp[2]-ym[2])/(2*dt)
            hp=u/np.cosh(np.sqrt(u)*t)**2
            fdtrace=-4*yc[1]**2*(hp+yc[1]**2)-beta*(rsecond+3*yc[1]*rfirst)
            trace_error=max(trace_error,abs(fdtrace-tc))
            constraint_identity=max(constraint_identity,abs((cp-cm)/(2*dt)+4*yc[1]*cc))
        fd.append({"step":dt,"continuity_residual":float(continuity),"trace_residual":float(trace_error),"constraint_evolution_residual":float(constraint_identity)})
    center=bounce_source(bounce_exact(0.,u),u)
    arrays.update({f"{prefix}_time":BOUNCE_TIME,f"{prefix}_a":y[0],f"{prefix}_h":y[1],f"{prefix}_Ricci":y[2],
        f"{prefix}_Ricci_prime":y[3],f"{prefix}_exact":exact.T,f"{prefix}_Radau":ry.T,
        f"{prefix}_rho":rho,f"{prefix}_pressure":p,f"{prefix}_constraint":constraint,
        f"{prefix}_curvature_density":curv,f"{prefix}_anomaly_density":y[1]**4,
        f"{prefix}_R2_density":curv-y[1]**4,
        f"{prefix}_radiation_density":rad,f"{prefix}_Kretschmann":kretschmann})
    result={"u":u,"beta":beta,"ell":ell,"D":D,"interval":[-6.,6.],
        "initial_state":[1.,0.,6*u,0.],"nfev":nfev,"Radau_nfev":nrad,
        "scaled_error_vs_exact":error,"Radau_scaled_error_vs_exact":radau_error,
        "constraint_max_absolute":float(np.max(abs(constraint))),"a_min":float(y[0].min()),
        "Ricci_range":[float(y[2].min()),float(y[2].max())],
        "Kretschmann_range":[float(kretschmann.min()),float(kretschmann.max())],
        "bounce_rho":float(center[0]),"bounce_pressure":float(center[1]),
        "bounce_curvature_density":float(center[4]),"bounce_radiation_density":float(center[5]),
        "finite_difference":fd,"tolerance_convergence":conv,
        "scope":"Exact bounce for the specified effective R^2 coupling. A compensated subtraction-convention change does not change the physics. No generic anomaly singularity-resolution claim or stability claim."}
    checks={"full_trace_ODE_vs_exact":error<BOUNDS["bounce_exact_scaled"],
        "independent_Radau_vs_exact":radau_error<BOUNDS["bounce_Radau_scaled"],
        "constraint":result["constraint_max_absolute"]<BOUNDS["bounce_constraint_absolute"],
        "positive_scale":result["a_min"]>0,
        "contraction_then_expansion":bool(np.all(y[1,BOUNCE_TIME<0]<0) and np.all(y[1,BOUNCE_TIME>0]>0)),
        "tolerance_convergence":conv[-1]["scaled_error_vs_exact"]<conv[0]["scaled_error_vs_exact"]/30,
        "source_continuity_FD":fd[-1]["continuity_residual"]<BOUNDS["bounce_FD_balance_absolute"],
        "trace_FD":fd[-1]["trace_residual"]<BOUNDS["bounce_FD_trace_absolute"],
        "source_cancellation":abs(center[0]+ell)+abs(center[4]+3*beta*u*u)+abs(center[5]-D)<BOUNDS["bounce_source_cancellation"]}
    result["checks"]={k:bool(v) for k,v in checks.items()}
    print(json.dumps(result,indent=2),flush=True)
    return result


def plot_results(root,stem,arrays,metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.08,top=.88,wspace=.29,hspace=.36)
    colors={0.:"#23658b",.02:"#b24d35"}
    for case in metrics["branches"]:
        label=case["label"];ell=case["ell"];branch=case["branch"]
        color=colors[ell];style="-" if branch==-1 else "--"
        legend=f"ell={ell:g}, {'unterer' if branch==-1 else 'oberer'} Ast"
        axes[0,0].plot(FUTURE_TIME,arrays[f"h_{label}"],color=color,ls=style,label=legend)
        axes[0,1].semilogy(FUTURE_TIME,arrays[f"a_{label}"],color=color,ls=style,label=legend)
        axes[1,0].semilogx(DELTAS,arrays[f"near_h_{label}"],color=color,ls=style,label=legend)
        axes[1,1].loglog(arrays[f"near_time_from_boundary_{label}"],arrays[f"near_Kretschmann_{label}"],color=color,ls=style,label=legend)
    axes[0,0].set(xlabel="Dimensionslose Zeit τ",ylabel="h = √α H",title="Beide expandierenden Zweige")
    axes[0,1].set(xlabel="Dimensionslose Zeit τ",ylabel="Skalenfaktor a",title="Selbstkonsistente Expansion")
    axes[1,0].axhline(1/np.sqrt(2),color="gray",ls=":",lw=1)
    axes[1,0].set(xlabel="δ = a/a* − 1",ylabel="h",title="Endliches h am vergangenen Rand")
    axes[1,1].set(xlabel="τ − τ*",ylabel="α² K",title="Krümmungspol: K ∝ (τ − τ*)⁻¹")
    for ax in axes.flat:ax.grid(alpha=.22);ax.legend(fontsize=7)
    fig.suptitle("Konforme renormierte Quelle: effektive R²-Kopplung β=0\nD=0,04; positives Strahlungsglied — kein regulärer Bounce an a*",fontsize=13)
    fig.savefig(root/"figures"/(stem+"_branches.png"),dpi=165)
    fig.savefig(root/"figures"/(stem+"_branches.svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.08,top=.88,wspace=.3,hspace=.36)
    tt=arrays["bounce_time"]
    axes[0,0].plot(tt,arrays["bounce_a"],color="#23658b",label="Volle Trace-ODE")
    axes[0,0].plot(tt,arrays["bounce_exact"][:,0],color="#b24d35",ls="--",label="Exakt: cosh τ")
    axes[0,0].set(xlabel="τ",ylabel="a",title="Global positive Geometrie")
    axes[0,1].plot(tt,arrays["bounce_h"],color="#23658b",label="h = tanh τ")
    axes[0,1].axhline(0,color="gray",lw=.8);axes[0,1].set(xlabel="τ",ylabel="h",title="Reguläre Kontraktion und Expansion")
    axes[1,0].plot(tt,arrays["bounce_radiation_density"],color="#b28725",label="Strahlung D/a⁴")
    axes[1,0].plot(tt,arrays["bounce_anomaly_density"],color="#3b8765",ls="--",label="Anomaliebeitrag h⁴")
    axes[1,0].plot(tt,arrays["bounce_R2_density"],color="#b24d35",label="Lokaler β-R²-Beitrag")
    axes[1,0].plot(tt,arrays["bounce_rho"],color="#23658b",label="Vollständige Quelle")
    axes[1,0].set(xlabel="τ",ylabel="α cG ρ",title="Quellenkompensation am Bounce")
    axes[1,1].plot(tt,arrays["bounce_Ricci"],color="#23658b",label="α R")
    axes[1,1].plot(tt,arrays["bounce_Kretschmann"],color="#b24d35",label="α² K")
    axes[1,1].set(xlabel="τ",ylabel="Krümmungsskalare",title="Beschränkte Krümmung")
    for ax in axes.flat:ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.suptitle("Exakter Bounce bei effektiver R²-Kopplung β=1/6\nell=0, D=1/2; minimale Skala a(0)=1",fontsize=13)
    fig.savefig(root/"figures"/(stem+"_bounce.png"),dpi=165)
    fig.savefig(root/"figures"/(stem+"_bounce.svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    fig.subplots_adjust(left=.09,right=.97,bottom=.17,top=.80,wspace=.30)
    axes[0].plot(tt,arrays["vacuum_bounce_h"],color="#23658b",label="h, volle Trace-ODE")
    axes[0].plot(tt,arrays["vacuum_bounce_exact"][:,1],color="#b24d35",ls="--",label="Exakt: √u tanh(√u τ)")
    axes[0].axhline(0,color="gray",lw=.8)
    axes[0].set(xlabel="τ",ylabel="h",title="Regulärer Bounce im konformen Vakuum")
    axes[1].plot(tt,arrays["vacuum_bounce_rho"],color="#b24d35",label="Renormierte Vakuumquelle")
    axes[1].axhline(3/16,color="#b28725",ls="--",label="Kosmologische Konstante ell")
    axes[1].plot(tt,arrays["vacuum_bounce_rho"]+3/16,color="#23658b",label="Summe = h²")
    axes[1].axhline(0,color="gray",lw=.8)
    axes[1].set(xlabel="τ",ylabel="Dimensionslose Quelldichte",title="Negative Vakuumenergie am Minimum")
    for ax in axes:ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.suptitle("Reines Vakuumbeispiel: u=3/4, β=1/9, ell=3/16, D=0\na=cosh(√u τ); keine Strahlungsanregung",fontsize=12)
    fig.savefig(root/"figures"/(stem+"_vacuum_bounce.png"),dpi=165)
    fig.savefig(root/"figures"/(stem+"_vacuum_bounce.svg"),metadata={"Date":None});plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    fig.subplots_adjust(left=.09,right=.97,bottom=.18,top=.84,wspace=.3)
    for case in metrics["branches"]:
        rows=case["tolerance_convergence"]
        axes[0].loglog([r["rtol"] for r in rows],[r["scaled_error_vs_independent_quadrature"] for r in rows],"o-",label=case["label"])
    for method,case in [("u=1","bounce"),("u=3/4","vacuum_bounce")]:
        rows=metrics[case]["tolerance_convergence"]
        axes[1].loglog([r["rtol"] for r in rows],[r["scaled_error_vs_exact"] for r in rows],"o-",label=method)
    axes[0].set(xlabel="ODE-Relativtoleranz",ylabel="Skalierter Fehler",title="β=0: Vergleich mit unabhängiger Quadratur")
    axes[1].set(xlabel="ODE-Relativtoleranz",ylabel="Skalierter Fehler",title="Vergleich mit beiden exakten Bounces")
    for ax in axes:ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.savefig(root/"figures"/(stem+"_validation.png"),dpi=165)
    fig.savefig(root/"figures"/(stem+"_validation.svg"),metadata={"Date":None});plt.close(fig)
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--tag",default="20261017")
    parser.add_argument("--pilot",action="store_true");args=parser.parse_args()
    if not args.tag.replace("_","").replace("-","").isalnum():parser.error("Invalid tag")
    root=Path(__file__).resolve().parents[1];stem="gravity_renormalized_"+args.tag
    paths=[root/"work"/(stem+e) for e in (".json",".npz",".md")]
    paths += [root/"figures"/(stem+s+e) for s in ("_branches","_bounce","_vacuum_bounce","_validation") for e in (".png",".svg")]
    if not args.pilot and any(p.exists() for p in paths):parser.error("Historical outputs exist; use a fresh tag")
    start=time.perf_counter();arrays={"future_time":FUTURE_TIME,"near_delta":DELTAS}
    metrics={"model":{"A":A_ANOMALY,"D":D_RADIATION,"ells":list(ELLS),"a0":1.,
        "future_interval":[0.,30.],"dimensionless_definitions":"alpha=c_G A; tau=t/sqrt(alpha), h=sqrt(alpha)H, ell=alpha Lambda/3, D=alpha c_G C, beta=b/A, rhohat=alpha c_G rho",
        "scope":"Massless conformally coupled real scalar, conformally flat FLRW, specified effective R^2 coupling; a compensated subtraction-convention change does not change the physics; no finite ultraviolet cutoff or experimental data fit"},
        "fixed_bounds":BOUNDS,"branches":[inspect_branch(ell,b,arrays) for ell in ELLS for b in (-1,1)],
        "deSitter":inspect_deSitter(),"bounce":inspect_bounce(arrays),"vacuum_bounce":inspect_bounce(arrays,.75)}
    groups=[x["checks"] for x in metrics["branches"]]+[metrics["deSitter"]["checks"],metrics["bounce"]["checks"],metrics["vacuum_bounce"]["checks"]]
    total=sum(len(x) for x in groups);passed=sum(sum(x.values()) for x in groups)
    metrics["checks_passed"]=passed;metrics["checks_total"]=total;metrics["all_checks_pass"]=passed==total
    if args.pilot:
        print(json.dumps({"passed":passed,"total":total,"all_checks_pass":passed==total,"elapsed":time.perf_counter()-start},indent=2));return
    mplversion=plot_results(root,stem,arrays,metrics);np.savez_compressed(paths[1],**arrays)
    metrics["provenance"]={"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":mplversion,"platform":platform.platform(),
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"arrays_sha256":hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        "command":f"python work/gravity_renormalized.py --tag {args.tag}","elapsed_seconds":time.perf_counter()-start}
    paths[0].write_text(json.dumps(metrics,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    report=f"""# Renormierte konforme FLRW-Quelle: numerische Gegenprüfungen

Lauf `{args.tag}`: **{passed}/{total}** feste Prüfungen bestanden.
Reproduktion: `python work/gravity_renormalized.py --tag NEUER_TAG`.

Das dimensionslose Modell verwendet die konforme Spuranomalie mit
`A=1/(960π²)`, `alpha=c_G A` und eine ausdrücklich festgelegte effektive
R²-Kopplung. Die Zeit ist `tau=t/sqrt(alpha)`, `h=sqrt(alpha) H`.
Für beta=0 werden D=.04 und ell=0,.02 auf beiden expandierenden Zweigen
betrachtet. Eine unabhängig durchgeführte Quadratur in log(a) prüft die
ODE. Die Vergangenheit wird nur bis delta=1e−6 numerisch integriert;
die Grenzasymptotik bis delta=1e−8 wird über den regulären Parameter
z=sqrt(discriminant) ausgewertet. Am Grenzpunkt bleiben a,h endlich,
während Ableitungen, Druck und Krümmung divergieren. Dies ist kein Bounce.

Die effektive Kopplung beta=1/6 erlaubt bei ell=0,D=.5 die exakte Lösung
a=cosh(tau), h=tanh(tau). Die vollständige Trace-ODE mit vier Variablen
wird vorwärts und rückwärts durch h=0 integriert, ohne durch h zu teilen.
Die lokale Krümmungsquelle kompensiert am Bounce die positive Strahlung.
Die Beispiele unterscheiden sich durch einen realen effektiven Modellparameter.
Ein bloßer Wechsel der Subtraktionskonvention mit entsprechender Kompensation
der Gravitationsparameter ändert die Physik nicht. Kein Ergebnis wird als
generischer Singularitätsauflösungssatz ausgegeben.

Die gleiche exakte Familie enthält für u=3/4 das reine Vakuumbeispiel
beta=1/9, ell=3/16, D=0, a=cosh(sqrt(u)*tau). Am Bounce ist die renormierte
Vakuumdichte -3/16; zusammen mit ell ergibt sich h²=0. Auch dieses Beispiel
wird mit DOP853, Radau und finiten Differenzen unabhängig kontrolliert.

Die C=0-de-Sitter-Kontrolle unterscheidet endliche affine Nullvergangenheit
des expandierenden flachen Patches von unendlicher komovierender Eigenzeit.
Der untere ell=0-Zweig ist Minkowski und wird ausdrücklich ausgenommen.

QuelleSHA256: `{metrics['provenance']['source_sha256']}`.
RoharraySHA256: `{metrics['provenance']['arrays_sha256']}`.
Alle Schwellen und Ergebnisse: `{paths[0].name}`.
"""
    paths[2].write_text(report,encoding="utf-8")
    print(json.dumps({"passed":passed,"total":total,"all_checks_pass":passed==total,"files":[str(p.relative_to(root)) for p in paths]},indent=2),flush=True)
    if passed!=total:raise SystemExit(1)


if __name__=="__main__":main()
