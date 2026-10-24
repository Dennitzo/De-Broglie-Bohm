"""Free minimally coupled massless scalar test modes on de Sitter FLRW.

python work/gravity_curved_field.py --tag 20261017

Natural units H=hbar=c=1, a=exp(t), spatially flat slicing, k>0 only.
The prescribed de Sitter metric is not solved with the actual Bohmian field as
an Einstein source. A nonzero classical stress-energy balance is measured.
Conformal time is a time-coordinate/lapse change of the same foliation.
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

WAVENUMBERS = np.array([.5, 1., 2., 4., 8.])
XI = np.array([.8, -1.1, .7, 1.25, -.9])
FINAL_TIME = 4.0
KG_START = -3.0
THRESHOLDS = {
    "analytic_ermakov_relative": 2e-12,
    "analytic_schrodinger_modal": 2e-11,
    "complex_kg_relative_error": 2e-8,
    "complex_kg_wronskian_error": 2e-8,
    "finite_difference_schrodinger_relative": 3e-5,
    "phase_integral_error": 2e-11,
    "rk4_fine_absolute": 2e-7,
    "lapse_covariance_absolute": 2e-8,
    "probability_log_jacobian_error": 2e-11,
    "deterministic_moment_error": 2e-11,
    "energy_identity_relative": 2e-12,
    "finite_difference_energy_absolute": 2e-6,
    "classical_energy_balance_absolute": 2e-10,
    "operator_energy_balance_absolute": 2e-10,
    "nonzero_bohm_source_minimum": 1e-2,
    "power_formula_relative": 2e-12,
}


def cosmic_parameters(t, k=WAVENUMBERS):
    a = np.exp(t)
    z = k*k/(a*a)
    mass = a**3
    sigma = np.sqrt((1+z)/(2*k**3))
    g = -z/(1+z)
    gdot = 2*z/(1+z)**2
    b = mass*g
    bdot = mass*(3*g+gdot)
    eta = -1/a
    theta = .5*(np.arctan(k*eta)-k*eta-np.arctan(-k)-k)
    return a, mass, z, sigma, g, gdot, b, bdot, theta


def conformal_parameters(eta, k=WAVENUMBERS):
    a = -1/eta
    effective_mass = a*a
    sigma = np.sqrt((1+k*k*eta*eta)/(2*k**3))
    geta = k*k*eta/(1+k*k*eta*eta)
    b = effective_mass*geta
    theta = .5*(np.arctan(k*eta)-k*eta-np.arctan(-k)-k)
    return a, effective_mass, sigma, geta, b, theta


def exact_mode_cosmic(t, k=WAVENUMBERS):
    eta = -np.exp(-t)
    coefficient = 1/np.sqrt(2*k**3)
    exponential = np.exp(-1j*k*eta)
    u = coefficient*(1+1j*k*eta)*exponential
    udot = -coefficient*k*k*eta*eta*exponential
    return u, udot


def exact_mode_conformal(eta, k=WAVENUMBERS):
    coefficient = 1/np.sqrt(2*k**3)
    exponential = np.exp(-1j*k*eta)
    return coefficient*(1+1j*k*eta)*exponential, coefficient*k*k*eta*exponential


def q_exact(t, xi=XI):
    return cosmic_parameters(t)[3]*xi


def velocity(t, q):
    return cosmic_parameters(t)[4]*q


def rk4(initial, dt):
    steps = int(round(FINAL_TIME/dt))
    assert abs(steps*dt-FINAL_TIME) < 1e-13
    result = np.empty((steps+1,len(WAVENUMBERS)))
    result[0] = initial
    for i in range(steps):
        t, q = i*dt, result[i]
        a = velocity(t,q)
        b = velocity(t+dt/2,q+dt*a/2)
        c = velocity(t+dt/2,q+dt*b/2)
        d = velocity(t+dt,q+dt*c)
        result[i+1] = q+dt*(a+2*b+2*c+d)/6
    return np.arange(steps+1)*dt,result


def wavefunction_log(coordinate, q, conformal=False):
    if conformal:
        _, _, sigma, _, b, theta = conformal_parameters(coordinate)
    else:
        *_, sigma, g, gd, b, bd, theta = cosmic_parameters(coordinate)
    return np.sum(-.25*np.log(2*np.pi*sigma*sigma)-q*q/(4*sigma*sigma)
                  +1j*(.5*b*q*q+theta))


def schrodinger_finite_difference(conformal=False):
    t = .37
    coordinate = -np.exp(-t) if conformal else t
    q = q_exact(t)
    psi = np.exp(wavefunction_log(coordinate,q,conformal))
    if conformal:
        a,effective_mass,*_ = conformal_parameters(coordinate)
        potential_coefficient = a*a*WAVENUMBERS**2
    else:
        a,effective_mass,*_ = cosmic_parameters(coordinate)
        potential_coefficient = a*WAVENUMBERS**2
    rows=[]
    for h in (.0004,.0002,.0001):
        dpsi = (np.exp(wavefunction_log(coordinate+h,q,conformal))
                -np.exp(wavefunction_log(coordinate-h,q,conformal)))/(2*h)
        laplacian=0j
        for j in range(len(q)):
            step=np.zeros_like(q);step[j]=h
            laplacian += (np.exp(wavefunction_log(coordinate,q+step,conformal))-2*psi
                          +np.exp(wavefunction_log(coordinate,q-step,conformal)))/h**2
        rhs=-laplacian/(2*effective_mass)+.5*np.sum(potential_coefficient*q*q)*psi
        rows.append({"h":h,"relative_residual":float(abs(1j*dpsi-rhs)/abs(psi))})
    return rows


def integrate_complex_modes(tgrid):
    count=len(WAVENUMBERS)
    u0,p0=exact_mode_cosmic(KG_START)
    def ode(t,y):
        u,p=y[:count],y[count:]
        return np.concatenate([p,-3*p-WAVENUMBERS**2*np.exp(-2*t)*u])
    rows=[];solutions=[]
    for tolerance in (1e-8,1e-10,1e-12):
        solution=solve_ivp(ode,(KG_START,FINAL_TIME),np.concatenate([u0,p0]),
                           method="DOP853",rtol=tolerance,atol=tolerance*1e-3,t_eval=tgrid)
        if not solution.success:raise RuntimeError(solution.message)
        numerical_u=solution.y[:count].T
        numerical_p=solution.y[count:].T
        exact_u=np.array([exact_mode_cosmic(t)[0] for t in tgrid])
        relative_error=float(np.max(abs(numerical_u-exact_u)/abs(exact_u)))
        wronskian=np.exp(3*tgrid[:,None])*(numerical_u*numerical_p.conj()-numerical_u.conj()*numerical_p)
        rows.append({"rtol":tolerance,"atol":tolerance*1e-3,"nfev":solution.nfev,
                     "relative_mode_error":relative_error,"wronskian_max_error":float(np.max(abs(wronskian-1j)))})
        solutions.append((numerical_u,numerical_p,wronskian))
    return rows,solutions


def independent_conformal_integration(tgrid,q0):
    count=len(WAVENUMBERS)
    u0,p0=exact_mode_conformal(-1.0)
    def ode(eta,y):
        u,p,q=y[:count],y[count:2*count],y[2*count:]
        return np.concatenate([p,2*p/eta-WAVENUMBERS**2*u,(p/u).real*q])
    etagrid=-np.exp(-tgrid)
    solution=solve_ivp(ode,(-1,float(etagrid[-1])),np.concatenate([u0,p0,q0.astype(complex)]),
                       method="DOP853",rtol=2e-13,atol=2e-15,t_eval=etagrid)
    if not solution.success:raise RuntimeError(solution.message)
    u,p,q=solution.y[:count].T,solution.y[count:2*count].T,solution.y[2*count:].T.real
    effective_mass=1/etagrid**2
    wronskian=effective_mass[:,None]*(u*p.conj()-u.conj()*p)
    return {"solver":"DOP853","rtol":2e-13,"atol":2e-15,"nfev":solution.nfev,
            "eta_interval":[-1,float(etagrid[-1])],
            "hamiltonian":"H_eta=sum[-1/(2a^2) d_q^2 + a^2 k^2 q^2/2] = a H_t",
            "wronskian_max_error":float(np.max(abs(wronskian-1j)))},(etagrid,u,p,q)


def classical_comparison(tgrid,q0):
    initial=np.concatenate([q0,velocity(0,q0)])
    def ode(t,y):
        q,p=y[:5],y[5:]
        return np.concatenate([p,-3*p-WAVENUMBERS**2*np.exp(-2*t)*q])
    solution=solve_ivp(ode,(0,FINAL_TIME),initial,method="DOP853",rtol=2e-13,atol=2e-15,t_eval=tgrid)
    if not solution.success:raise RuntimeError(solution.message)
    q,p=solution.y[:5].T,solution.y[5:].T
    acceleration=-3*p-WAVENUMBERS**2*np.exp(-2*tgrid[:,None])*q
    return q,p,acceleration


def energy_observables(tgrid,q,p,acceleration):
    z=WAVENUMBERS**2*np.exp(-2*tgrid[:,None])
    rho=.5*(p*p+z*q*q)
    pressure=.5*p*p-z*q*q/6
    rho_dot=p*acceleration-z*q*q+z*q*p
    balance=rho_dot+3*(rho+pressure)
    return rho,pressure,rho_dot,balance


def write_plots(root,stem,arrays,metrics):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,2,figsize=(12,8.4))
    fig.subplots_adjust(left=.08,right=.97,bottom=.08,top=.88,wspace=.29,hspace=.37)
    colors=plt.cm.viridis(np.linspace(.05,.92,len(WAVENUMBERS)))
    t=arrays["time"]
    for j,(k,color) in enumerate(zip(WAVENUMBERS,colors)):
        axes[0,0].plot(t,arrays["bohm_q"][:,j]/arrays["initial_q"][j],color=color,label=f"k={k:g}")
        axes[0,0].plot(t,arrays["classical_q"][:,j]/arrays["initial_q"][j],color=color,ls="--",lw=1)
        axes[0,1].semilogy(t,arrays["dimensionless_power"][:,j],color=color)
        axes[1,0].semilogy(t,-arrays["bohm_energy_balance"][:,j],color=color)
    axes[0,0].set(xlabel="Kosmische Zeit t",ylabel="q(t) / q(0)",title="Bohmsche Moden und klassische KG-Bahnen")
    legend1=axes[0,0].legend(ncol=3,fontsize=8,loc="upper right")
    axes[0,0].add_artist(legend1)
    axes[0,0].legend(handles=[Line2D([0],[0],color="black",label="Bohm"),
                              Line2D([0],[0],color="black",ls="--",label="klassisch")],fontsize=8,loc="lower right")
    limit=1/(4*np.pi**2)
    axes[0,1].axhline(limit,color="black",ls=":",label=r"$H^2/(4\pi^2)$")
    axes[0,1].set(xlabel="Kosmische Zeit t",ylabel=r"$\Delta_\phi^2=k^3\sigma_k^2/(2\pi^2)$",title="Modenleistung und Horizontgrenze")
    axes[0,1].legend(fontsize=9)
    axes[1,0].set(xlabel="Kosmische Zeit t",ylabel=r"$-[\dot\rho_B+3H(\rho_B+p_B)]$",title="Nichtverschwindender klassischer Quellenrest")
    axes[1,0].text(.04,.08,"Auf der Bohmbahn: −q̇ q/(4a⁶σ⁴)\nKlassische KG-Bahn: Rest = 0",transform=axes[1,0].transAxes,fontsize=9,
                    bbox={"facecolor":"white","alpha":.85,"edgecolor":"none"})
    j=2
    axes[1,1].plot(t,arrays["ensemble_bohm_balance"][:,j],label="Gaußmittel des Bohm-Quellenrests",color="#b44136")
    axes[1,1].plot(t,arrays["momentum_variance_balance"][:,j],label="Impulsvarianzbeitrag",color="#277a69",ls="--")
    axes[1,1].plot(t,arrays["operator_energy_balance"][:,j],label="Operator-Erwartungswert: Summe",color="black",lw=1.1)
    axes[1,1].set_yscale("symlog",linthresh=1e-9)
    axes[1,1].set_yticks([-1,-1e-3,-1e-6,0,1e-6,1e-3,1],
                        ["−1","−10⁻³","−10⁻⁶","0","10⁻⁶","10⁻³","1"])
    axes[1,1].set(xlabel="Kosmische Zeit t",ylabel="Energiequellenrest",title="Erhaltung des Operator-Mittelwerts (k=2)")
    axes[1,1].legend(fontsize=8,loc="upper right")
    for ax in axes.flat:ax.grid(alpha=.25)
    fig.suptitle("Bohmsches Skalar-Testfeld auf de-Sitter-FLRW\n"
                 "a(t)=eᵗ, H=1; gleiche Anfangswerte für Bohm und klassische Vergleichsbahn",fontsize=13)
    fig.savefig(root/"figures"/(stem+".png"),dpi=170)
    fig.savefig(root/"figures"/(stem+".svg"),metadata={"Date":None})
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.2))
    fig.subplots_adjust(left=.09,right=.98,bottom=.17,top=.86,wspace=.29)
    rows=metrics["rk4_convergence"]
    axes[0].loglog([r["dt"] for r in rows],[r["max_absolute_error"] for r in rows],"o-")
    axes[0].set(xlabel="Schrittweite Δt",ylabel="Maximaler Modenfehler",title="RK4 gegen exakten Gaußfluss")
    kg=metrics["complex_KG_convergence"]
    axes[1].loglog([r["rtol"] for r in kg],[r["relative_mode_error"] for r in kg],"o-",label="KG-Mode")
    axes[1].loglog([r["rtol"] for r in kg],[r["wronskian_max_error"] for r in kg],"s-",label="Kanonischer Wronskian")
    axes[1].set(xlabel="DOP853-Relativtoleranz",ylabel="Maximaler Fehler",title="Unabhängige Klein-Gordon-Integration")
    axes[1].legend(fontsize=9)
    for ax in axes:ax.grid(alpha=.25)
    fig.savefig(root/"figures"/(stem+"_validation.png"),dpi=170)
    fig.savefig(root/"figures"/(stem+"_validation.svg"),metadata={"Date":None})
    plt.close(fig)
    return matplotlib.__version__


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag",default="20261017")
    parser.add_argument("--pilot",action="store_true")
    args=parser.parse_args()
    if not args.tag.replace("_","").replace("-","").isalnum():parser.error("invalid tag")
    root=Path(__file__).resolve().parents[1];stem="gravity_curved_field_"+args.tag
    paths=[root/"work"/(stem+e) for e in (".json",".npz",".md")]+[root/"figures"/(stem+e) for e in (".png",".svg","_validation.png","_validation.svg")]
    if not args.pilot and any(p.exists() for p in paths):parser.error("Historical output exists; choose another tag.")
    start=time.perf_counter()
    tgrid=np.linspace(0,FINAL_TIME,401)
    parameters=[cosmic_parameters(t) for t in tgrid]
    sigma=np.array([p[3] for p in parameters]);g=np.array([p[4] for p in parameters]);gd=np.array([p[5] for p in parameters])
    q=sigma*XI;p=g*q;acceleration=(g*g+gd)*q
    q0=q[0]
    arrays={"time":tgrid,"k":WAVENUMBERS,"xi":XI,"sigma":sigma,"bohm_q":q,"bohm_qdot":p,
            "bohm_qddot":acceleration,"initial_q":q0}
    ermakov_error=schrodinger_error=0.0
    for i,t in enumerate(tgrid):
        a,m,z,s,gg,ggd,b,bd,theta=parameters[i]
        ermakov=(gg*gg+ggd)*s+3*gg*s+z*s-1/(4*m*m*s**3)
        ermakov_error=max(ermakov_error,float(np.max(abs(ermakov)/(1+z*s))))
        potential=.5*m*z*q[i]**2
        quantum=1/(4*m*s*s)-q[i]**2/(8*m*s**4)
        phase_dt=.5*bd*q[i]**2-1/(4*m*s*s)
        real=phase_dt+b*b*q[i]**2/(2*m)+potential+quantum
        imag=-.5*gg+q[i]**2*gg/(2*s*s)-b*q[i]**2/(2*m*s*s)+b/(2*m)
        schrodinger_error=max(schrodinger_error,float(np.max(abs(real+1j*imag))))
    metrics={"model":{"background":"spatially flat de Sitter FLRW, prescribed test background",
                       "H":1,"hbar":1,"c":1,"scalar_mass":0,"curvature_coupling":0,
                       "k":WAVENUMBERS.tolist(),"zero_mode":"excluded",
                       "initial_standardized_coordinates":XI.tolist(),"time_interval":[0,FINAL_TIME],
                       "KG_integration_start":KG_START,"scale_factor":"a=exp(t)",
                       "conformal_time":"eta=-exp(-t)","lapse":"dt/deta=a",
                       "energy_average":"spatially averaged real-mode contribution, comoving normalization volume=1; p is one third of the spatial stress trace, not a local perfect-fluid assertion"},
             "acceptance_thresholds":THRESHOLDS,
             "analytic":{"ermakov_scaled_max_error":ermakov_error,"schrodinger_modal_max_error":schrodinger_error},
             "schrodinger_finite_difference_cosmic":schrodinger_finite_difference(False),
             "schrodinger_finite_difference_conformal":schrodinger_finite_difference(True)}
    phase_error=0.0
    for j,k in enumerate(WAVENUMBERS):
        for tau in (.37,1.,4.):
            integrated=quad(lambda t:-1/(4*cosmic_parameters(t)[1]*cosmic_parameters(t)[3][j]**2),
                            0,tau,epsabs=2e-13,epsrel=2e-13)[0]
            phase_error=max(phase_error,abs(integrated-cosmic_parameters(tau)[-1][j]))
    metrics["analytic"]["global_phase_integral_error"]=float(phase_error)
    kgrows,kgsolutions=integrate_complex_modes(tgrid)
    metrics["complex_KG_convergence"]=kgrows
    arrays.update(KG_u=kgsolutions[-1][0],KG_udot=kgsolutions[-1][1],KG_wronskian=kgsolutions[-1][2])
    rows=[]
    for dt in (.08,.04,.02):
        ts,trajectory=rk4(q0,dt)
        exact=np.array([q_exact(t) for t in ts])
        rows.append({"dt":dt,"max_absolute_error":float(np.max(abs(trajectory-exact)))})
        arrays[f"rk4_times_dt{dt:g}"]=ts;arrays[f"rk4_q_dt{dt:g}"]=trajectory
    metrics["rk4_convergence"]=rows
    metrics["rk4_orders"]=[float(np.log2(rows[i]["max_absolute_error"]/rows[i+1]["max_absolute_error"])) for i in (0,1)]
    conformal_info,(etas,eta_u,eta_up,eta_q)=independent_conformal_integration(tgrid,q0)
    conformal_info["trajectory_covariance_absolute_error"]=float(np.max(abs(eta_q-q)))
    conformal_info["phase_coefficient_b_error"]=float(max(np.max(abs(conformal_parameters(eta)[4]-cosmic_parameters(t)[6])) for t,eta in zip(tgrid,etas)))
    arrays.update(conformal_eta=etas,conformal_u=eta_u,conformal_uprime=eta_up,conformal_q=eta_q)
    metrics["conformal_time_covariance"]=conformal_info
    classical_q,classical_p,classical_acceleration=classical_comparison(tgrid,q0)
    arrays.update(classical_q=classical_q,classical_qdot=classical_p,classical_qddot=classical_acceleration)
    rho,pressure,rhodot,balance=energy_observables(tgrid,q,p,acceleration)
    classical_rho,classical_pressure,classical_rhodot,classical_balance=energy_observables(tgrid,classical_q,classical_p,classical_acceleration)
    quantum_exchange=p*q/(4*np.exp(6*tgrid[:,None])*sigma**4)
    energy_error=float(np.max(abs(balance-quantum_exchange))/np.max(abs(quantum_exchange)))
    arrays.update(bohm_rho=rho,bohm_pressure=pressure,bohm_rhodot=rhodot,bohm_energy_balance=balance,
                  quantum_exchange=quantum_exchange,classical_rho=classical_rho,classical_pressure=classical_pressure,
                  classical_energy_balance=classical_balance)
    fd_energy=[]
    for h in (.002,.001,.0005):
        tau=.73
        def rho_at(s):
            _,_,z,ss,gg,*_=cosmic_parameters(s);qq=ss*XI
            return .5*((gg*qq)**2+z*qq*qq)
        rd=(rho_at(tau+h)-rho_at(tau-h))/(2*h)
        _,m,z,ss,gg,*_=cosmic_parameters(tau);qq=ss*XI;pp=gg*qq
        rho_here=.5*(pp*pp+z*qq*qq);pressure_here=.5*pp*pp-z*qq*qq/6
        exact_exchange=pp*qq/(4*m*m*ss**4)
        fd_energy.append({"h":h,"absolute_balance_error":float(np.max(abs(rd+3*(rho_here+pressure_here)-exact_exchange)))})
    metrics["energy_balance"]={"identity_relative_error":energy_error,
        "classical_balance_max_abs":float(np.max(abs(classical_balance))),
        "bohm_balance_max_abs":float(np.max(abs(balance))),
        "initial_bohm_balance":balance[0].tolist(),"final_bohm_balance":balance[-1].tolist(),
        "finite_difference_checks":fd_energy,
        "interpretation":"Nonzero classical stress-energy divergence on actual Bohm trajectories; no self-consistent Einstein solution is asserted."}
    # Quantum momentum variance supplies an independently defined finite-mode
    # operator expectation. It is not an invented counterterm for an actual path.
    mean_rho,mean_pressure,mean_rhodot,mean_balance=energy_observables(tgrid,sigma,g*sigma,(g*g+gd)*sigma)
    delta=1/(8*np.exp(6*tgrid[:,None])*sigma**2)
    delta_dot=delta*(-6-2*g)
    delta_balance=delta_dot+6*delta
    operator_balance=mean_balance+delta_balance
    operator_rho=mean_rho+delta
    operator_pressure=mean_pressure+delta
    arrays.update(ensemble_bohm_rho=mean_rho,ensemble_bohm_pressure=mean_pressure,ensemble_bohm_balance=mean_balance,
                  momentum_variance_density=delta,momentum_variance_balance=delta_balance,
                  operator_rho=operator_rho,operator_pressure=operator_pressure,operator_energy_balance=operator_balance)
    fd_operator=[]
    for h in (.002,.001,.0005):
        tau=.73
        def operator_at(s):
            _,m,z,ss,gg,*_=cosmic_parameters(s)
            correction=1/(8*m*m*ss*ss)
            return .5*((gg*ss)**2+z*ss*ss)+correction,.5*(gg*ss)**2-z*ss*ss/6+correction
        derivative=(operator_at(tau+h)[0]-operator_at(tau-h)[0])/(2*h)
        rho_here,pressure_here=operator_at(tau)
        fd_operator.append({"h":h,"absolute_balance_error":float(np.max(abs(derivative+3*(rho_here+pressure_here))))})
    metrics["operator_stress_expectation"]={
        "delta_rho_equals_delta_p":"1/(8a^6 sigma^2), from quantum momentum variance",
        "mean_bohm_source_max_abs":float(np.max(abs(mean_balance))),
        "operator_balance_max_abs":float(np.max(abs(operator_balance))),
        "correction_exchange_identity_max_abs":float(np.max(abs(delta_balance+g/(4*np.exp(6*tgrid[:,None])*sigma**2)))),
        "finite_difference_checks":fd_operator,
        "interpretation":"Finite-mode operator stress expectation is conserved; no individual-field counterterm or UV-renormalized infinite stress tensor is asserted."}
    logflow_error=0.0
    logp0=-.5*np.log(2*np.pi*sigma[0]**2).sum()-.5*np.sum((q0/sigma[0])**2)
    logdet=np.log(sigma/sigma[0]).sum(axis=1)
    for i in range(len(tgrid)):
        logp=-.5*np.log(2*np.pi*sigma[i]**2).sum()-.5*np.sum((q[i]/sigma[i])**2)
        logflow_error=max(logflow_error,abs(logp+logdet[i]-logp0))
    nodes,weights=np.polynomial.hermite.hermgauss(5);nodes*=np.sqrt(2);weights/=np.sqrt(np.pi)
    variance_measured=((sigma[:,:,None]*nodes)**2)@weights
    moment_error=float(np.max(abs(variance_measured-sigma*sigma)))
    arrays.update(flow_logdet=logdet,gauss_hermite_nodes=nodes,gauss_hermite_weights=weights)
    metrics["equivariance"]={"log_density_jacobian_error":float(logflow_error),
        "deterministic_variance_error":moment_error,"gauss_hermite_nodes":5,
        "min_width_on_sampled_interval":float(sigma.min()),"q_over_sigma_error":float(np.max(abs(q/sigma-XI)))}
    power=WAVENUMBERS**3*sigma**2/(2*np.pi**2)
    power_formula=(1+WAVENUMBERS**2*np.exp(-2*tgrid[:,None]))/(4*np.pi**2)
    arrays["dimensionless_power"]=power
    metrics["dimensionless_power"]={"formula_relative_error":float(np.max(abs(power-power_formula)/power_formula)),
        "late_limit":float(1/(4*np.pi**2)),"initial_values":power[0].tolist(),"final_values":power[-1].tolist(),
        "final_relative_excess_above_limit":(power[-1]*4*np.pi**2-1).tolist(),
        "horizon_crossing_times":np.log(WAVENUMBERS).tolist()}
    th=THRESHOLDS
    checks={
        "ermakov_equation":ermakov_error<th["analytic_ermakov_relative"],
        "analytic_modal_schrodinger":schrodinger_error<th["analytic_schrodinger_modal"],
        "five_mode_schrodinger_FD_cosmic":metrics["schrodinger_finite_difference_cosmic"][-1]["relative_residual"]<th["finite_difference_schrodinger_relative"],
        "five_mode_schrodinger_FD_conformal":metrics["schrodinger_finite_difference_conformal"][-1]["relative_residual"]<th["finite_difference_schrodinger_relative"],
        "global_phase":phase_error<th["phase_integral_error"],
        "complex_KG_independent_ODE":kgrows[-1]["relative_mode_error"]<th["complex_kg_relative_error"],
        "cosmic_canonical_wronskian":kgrows[-1]["wronskian_max_error"]<th["complex_kg_wronskian_error"],
        "conformal_canonical_wronskian":conformal_info["wronskian_max_error"]<th["complex_kg_wronskian_error"],
        "rk4_guidance":rows[-1]["max_absolute_error"]<th["rk4_fine_absolute"],
        "same_foliation_time_covariance":conformal_info["trajectory_covariance_absolute_error"]<th["lapse_covariance_absolute"],
        "density_jacobian_equivariance":logflow_error<th["probability_log_jacobian_error"],
        "gaussian_deterministic_moments":moment_error<th["deterministic_moment_error"],
        "energy_exchange_identity":energy_error<th["energy_identity_relative"],
        "energy_exchange_FD":fd_energy[-1]["absolute_balance_error"]<th["finite_difference_energy_absolute"],
        "classical_stress_balance":np.max(abs(classical_balance))<th["classical_energy_balance_absolute"],
        "finite_mode_operator_stress_conservation":np.max(abs(operator_balance))<th["operator_energy_balance_absolute"],
        "operator_stress_FD":fd_operator[-1]["absolute_balance_error"]<th["finite_difference_energy_absolute"],
        "bohm_classical_source_nonzero":np.max(abs(balance))>th["nonzero_bohm_source_minimum"],
        "dimensionless_power_formula":metrics["dimensionless_power"]["formula_relative_error"]<th["power_formula_relative"],
    }
    metrics["acceptance_checks"]={key:bool(value) for key,value in checks.items()}
    metrics["all_checks_pass"]=all(checks.values())
    if args.pilot:
        print(json.dumps(metrics,indent=2),flush=True);return
    version=write_plots(root,stem,arrays,metrics)
    np.savez_compressed(paths[1],**arrays)
    metrics["provenance"]={"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"scipy":scipy.__version__,
        "matplotlib":version,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arrays_sha256":hashlib.sha256(paths[1].read_bytes()).hexdigest(),
        "command":f"python work/gravity_curved_field.py --tag {args.tag}","elapsed_seconds":time.perf_counter()-start}
    paths[0].write_text(json.dumps(metrics,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    report=f"""# Bohmsches Testfeld auf de-Sitter-FLRW

Verbindlicher Lauf: **{args.tag}**, **{sum(checks.values())}/{len(checks)}** feste Prüfungen bestanden.
Reproduktion: `python work/gravity_curved_field.py --tag NEUER_TAG`.

Das räumlich flache FLRW-Testmodell setzt a(t)=exp(t), H=ħ=c=1 und verwendet ein
freies massloses minimal gekoppeltes reelles Skalarfeld. Die fünf positiven
Wellenzahlen sind .5,1,2,4,8; die Nullmode ist ausgeschlossen. Die Metrik bleibt
vorgegeben. Diese Rechnung ist keine gekoppelte Einstein-Materie-Lösung.

Für jede kanonische Modenkoordinate gilt M=a³, Ω²=k²/a² und
σ²=(1+k²/a²)/(2k³). Die normierte Gaußwellenfunktion hat Phase bq²/2+θ,
b=M σdot/σ. Mit η=−exp(−t) lautet die kontinuierliche globale Phase
θ=[atan(kη)−kη−atan(−k)−k]/2, sodass θ(0)=0.
Der exakte Modenfluss ist q(t)=σ(t)ξ. Er wird mit RK4, einer unabhängigen
komplexen Klein-Gordon-Integration ab t=−3 sowie einer separaten gemeinsamen
KG-/Führungsintegration in Konformzeit verglichen. Dort lautet der Hamiltonian
Hη=Σ[−∂q²/(2a²)+a²k²q²/2]=aHt. Dies prüft die Zeitkoordinate und Lapse
derselben Foliation, keine Unabhängigkeit von einer anderen Foliation.

Die klassische Vergleichsbahn löst die homogene KG-Gleichung mit denselben
q(0),qdot(0). Auf der Bohmbahn bleibt dagegen
qddot+3Hqdot+Ω²q=q/(4a⁶σ⁴). Der räumlich isotrop gemittelte klassische
Energiequellenrest beträgt rho_dot+3H(rho+p)=qdot*q/(4a⁶σ⁴).
Dabei sind rho und p räumlich gemittelte Einmodenbeiträge für komovierendes
Normierungsvolumen 1; p ist ein Drittel der räumlichen Spannungsspur. Eine
einzelne Mode wird nicht als lokales isotropes perfektes Fluid behandelt.
Ein klassischer Stress-Energie-Tensor dieser tatsächlichen Bohmbahn allein
erfüllt somit die für Einsteins Gleichung notwendige Erhaltung nicht.
Das Quellenglied wird ausdrücklich gemessen, nicht durch eine behauptete
unbekannte gravitative Zusatzquelle entfernt.

Der endliche Operator-Erwartungswert unterscheidet sich vom Gaußmittel des
klassischen Tensors auf den Bohmbahnen: Die quantenmechanische Impulsvarianz
liefert delta-rho=delta-p=1/(8a⁶σ²). Dieser fest definierte Operatorbeitrag
kompensiert den gemittelten Quellenrest exakt. Es handelt sich nicht um einen
frei hinzugefügten Term für die individuelle Feldrealisierung. Der volle
UV-renormierte unendlichdimensionale Stress-Energie-Tensor wird nicht berechnet.

- Maximales skaliertes Ermakov-Residuum: {ermakov_error:.6e}.
- Maximales analytisches modales Schrödingerresiduum: {schrodinger_error:.6e}.
- Vollständiges Fünfmoden-FD-Residuum kosmisch/konform:
  {metrics['schrodinger_finite_difference_cosmic'][-1]['relative_residual']:.6e} /
  {metrics['schrodinger_finite_difference_conformal'][-1]['relative_residual']:.6e}.
- Komplexe KG-ODE: relativer Fehler {kgrows[-1]['relative_mode_error']:.6e}.
- Zeitkoordinatenvergleich: {conformal_info['trajectory_covariance_absolute_error']:.6e}.
- RK4-Fehler: {rows[-1]['max_absolute_error']:.6e}.
- Quellenidentität, relativer Fehler: {energy_error:.6e}.

Das Modenspektrum Δφ²=k³σ²/(2π²)=(1+k²/a²)/(4π²) nähert sich
H²/(4π²). Der Verlauf zeigt den UV-/Horizontterm ausdrücklich. Eine
IR-vollständige de-Sitter-invariante Vakuumbehauptung wird durch diese fünf
positiven Moden nicht gemacht. Es werden keine experimentellen Daten erzeugt.

QuelleSHA256: `{metrics['provenance']['source_sha256']}`.
RoharraySHA256: `{metrics['provenance']['arrays_sha256']}`.
Alle Werte, Konvergenzreihen und Schwellen stehen in `{paths[0].name}`.
"""
    paths[2].write_text(report,encoding="utf-8")
    print(json.dumps({"all_checks_pass":all(checks.values()),"checks":metrics["acceptance_checks"],
                      "files":[str(p.relative_to(root)) for p in paths]},indent=2),flush=True)
    if not all(checks.values()):raise SystemExit(1)


if __name__=="__main__":main()
