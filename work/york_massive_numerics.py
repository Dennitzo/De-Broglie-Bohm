# -*- coding: utf-8 -*-
"""
Phase 0.5 / Phase 1.3: Massive York-Dynamik mit Ereigniserkennung.

Modell:
  FLRW k=0, Lambda=0, G = 1, kappa = 8*pi*G,
  V(Phi) = 1/2 m^2 Phi^2.

Kosmische Zeit (N=1):
  adot/a = H
  Phi_dot = v
  v_dot = -3 H v - m^2 Phi
  H_dot = -4*pi*G*v^2
  Friedmann: H^2 = (8*pi*G/3)*(v^2/2 + m^2 Phi^2/2)

York-Zeit (Expansionzweig, H>0):
  T = -H/(4*pi*G),  T < 0
  dT/dt = v^2
  A(T,Phi) = 12*pi*G*T^2 - 2V = 12*pi*G*T^2 - m^2 Phi^2 = v^2

Reduzierter York-Hamiltonoperator auf einem Zweig sigma = sign(p):
  H_T = |p|/sqrt(A), p = a^3 v
  Phi' = sigma/sqrt(A)
  p' = -sigma*p*V'(Phi)/A^(3/2)

Der Code integriert die kosmischen Gleichungen durch die Umkehrpunkte hindurch
und die York-Gleichungen nur auf Intervallen mit A>0. Eine York-Integration
ueber eine Singularitaet A=0 wird durch Ereignisterminierung verhindert.
"""
import json
import math
import os
import datetime

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "york_massive_numerics.json")

G = 1.0
KAPPA = 8.0 * math.pi * G
m = 1.0

phi0 = 1.0
v0 = -0.7
a0 = 1.0
rho0 = 0.5 * v0**2 + 0.5 * m**2 * phi0**2
H0 = math.sqrt(KAPPA * rho0 / 3.0)
t_end = 28.0

results = []
intervals = []
events = []


def record(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("PASS" if ok else "FAIL"), "|", name, "|", detail)


def cosmic_rhs(t, y):
    a, phi, v, H = y
    adot = a * H
    phidot = v
    vdot = -3.0 * H * v - m**2 * phi
    Hdot = -4.0 * math.pi * G * v**2
    return [adot, phidot, vdot, Hdot]


def turning_event(t, y):
    return y[2]


turning_event.terminal = False
turning_event.direction = 0.0


y0 = [a0, phi0, v0, H0]
cosmic_sol = solve_ivp(
    cosmic_rhs,
    (0.0, t_end),
    y0,
    dense_output=True,
    rtol=1e-11,
    atol=1e-12,
    max_step=0.02,
    events=[turning_event],
)

if not cosmic_sol.success:
    record("cosmische Integration erfolgreich", False, cosmic_sol.message)
    raise SystemExit(1)

record("cosmische Integration erfolgreich", True, cosmic_sol.message)

raw_event_times = list(cosmic_sol.t_events[0])
filtered_events = []
for te in raw_event_times:
    if not filtered_events or abs(te - filtered_events[-1]) > 1e-7:
        filtered_events.append(float(te))
events = filtered_events
record(
    "massives Potential: numerische Umkehrpunkte (A -> 0 / Phi_dot = 0) gefunden",
    len(events) >= 2,
    "Anzahl Ereignisse = %d; Zeiten = %s" % (len(events), [round(x, 6) for x in events[:8]]),
)

# ----------------- globale Kontrollen auf dichtem kosmischem Gitter --------
grid_t = np.linspace(0.0, t_end, 14000)
grid = cosmic_sol.sol(grid_t)
a_cos, phi_cos, v_cos, H_cos = grid
T_cos = -H_cos / (4.0 * math.pi * G)
rho_cos = 0.5 * v_cos**2 + 0.5 * m**2 * phi_cos**2
friedmann_res = H_cos**2 - KAPPA * rho_cos / 3.0
A_cos = 12.0 * math.pi * G * T_cos**2 - m**2 * phi_cos**2
A_identity_res = A_cos - v_cos**2

friedmann_max = float(np.max(np.abs(friedmann_res)))
A_identity_max = float(np.max(np.abs(A_identity_res)))
dT_min = float(np.min(np.diff(T_cos)))
record(
    "Friedmann-Zwang entlang der kosmischen Zeit",
    friedmann_max < 1e-6,
    "max |H^2 - kappa rho/3| = %.3e" % friedmann_max,
)
record(
    "York-Identitaet A(T,Phi) = Phi_dot^2",
    A_identity_max < 1e-6,
    "max |A - v^2| = %.3e" % A_identity_max,
)
record(
    "York-Zeit ist monoton nichtfallend (dT/dt = v^2)",
    dT_min > -1e-9,
    "min Delta T = %.3e; T-Bereich = [%.6f, %.6f]" % (dT_min, float(T_cos[0]), float(T_cos[-1])),
)

event_A_values = []
event_v_values = []
for te in events:
    state = cosmic_sol.sol(te)
    _, phi_e, v_e, H_e = state
    T_e = -H_e / (4.0 * math.pi * G)
    A_e = 12.0 * math.pi * G * T_e**2 - m**2 * phi_e**2
    event_A_values.append(float(A_e))
    event_v_values.append(float(v_e))

event_A_max = max(abs(x) for x in event_A_values) if event_A_values else float("nan")
event_v_max = max(abs(x) for x in event_v_values) if event_v_values else float("nan")
record(
    "An den Ereignissen gilt Phi_dot = 0 und A = 0 (numerisch)",
    event_v_max < 1e-6 and event_A_max < 1e-4,
    "max |v(event)| = %.3e; max |A(event)| = %.3e" % (event_v_max, event_A_max),
)


def A_of_T_phi(T, phi):
    return 12.0 * math.pi * G * T**2 - m**2 * phi**2


def york_rhs_factory(sigma):
    def york_rhs(t, y):
        phi, p = y
        A = A_of_T_phi(t, phi)
        A_safe = max(A, 1e-14)
        sqrtA = math.sqrt(A_safe)
        dphi = sigma / sqrtA
        dp = -sigma * p * (m**2 * phi) / (A_safe * sqrtA)
        return [dphi, dp]

    def event_A(t, y):
        return A_of_T_phi(t, y[1]) - 1e-9

    event_A.terminal = True
    event_A.direction = -1.0
    return york_rhs, event_A


# ----------------- York-Integration auf Regularitaetsintervallen ------------
boundaries = [0.0] + events + [float(cosmic_sol.t[-1])]
max_phi_err = 0.0
max_p_err = 0.0
max_p_rel_err = 0.0
max_a_err = 0.0
min_abs_v_regular = float("inf")
comparison_count = 0

for interval_index in range(len(boundaries) - 1):
    left = boundaries[interval_index]
    right = boundaries[interval_index + 1]
    span = right - left
    if span <= 1e-4:
        continue
    eps = max(0.18 * span, 1e-4)
    tl = left + eps
    tr = right - eps
    if tr - tl <= 1e-4:
        continue

    st_l = cosmic_sol.sol(tl)
    st_r = cosmic_sol.sol(tr)
    T_l = -st_l[3] / (4.0 * math.pi * G)
    T_r = -st_r[3] / (4.0 * math.pi * G)
    if not (math.isfinite(T_l) and math.isfinite(T_r)) or T_r <= T_l + 1e-8:
        continue

    a_l, phi_l, v_l, _ = st_l
    p_l = a_l**3 * v_l
    sigma = 1.0 if p_l > 0.0 else -1.0

    interval_t = np.linspace(tl, tr, 601)
    interval_states = cosmic_sol.sol(interval_t)
    ia, iphi, iv, iH = interval_states
    min_abs_v = float(np.min(np.abs(iv)))
    if min_abs_v < min_abs_v_regular:
        min_abs_v_regular = min_abs_v

    A_l = A_of_T_phi(T_l, phi_l)
    if A_l <= 1e-10:
        continue

    rhs, evA = york_rhs_factory(sigma)
    try:
        yk = solve_ivp(
            rhs,
            (float(T_l), float(T_r)),
            [float(phi_l), float(p_l)],
            dense_output=True,
            rtol=1e-10,
            atol=1e-12,
            max_step=(float(T_r) - float(T_l)) / 800.0,
            events=[evA],
        )
    except Exception as exc:
        record(
            "York-Integration Intervall %d" % interval_index,
            False,
            "Exception: %s" % exc,
        )
        continue

    yk_t_stop = float(yk.t[-1])
    usable_T_max = min(float(T_r), yk_t_stop - 1e-6)
    usable_T_min = float(T_l) + 1e-8
    if usable_T_max <= usable_T_min:
        continue

    # Vergleichspunkte aus der kosmischen Loesung, dort York-Dichte auswerten.
    comp_times = interval_t[
        (interval_t >= tl + 1e-6) & (interval_t <= tr - 1e-6)
    ]
    if comp_times.size < 20:
        comp_times = interval_t
    comp_states = cosmic_sol.sol(comp_times)
    ca, cphi, cv, cH = comp_states
    cT = -cH / (4.0 * math.pi * G)
    mask = (cT >= usable_T_min) & (cT <= usable_T_max)
    if np.count_nonzero(mask) < 20:
        continue

    T_samples = cT[mask]
    # T_samples ist monoton, aber minimal bereinigen.
    order = np.argsort(T_samples)
    T_samples = np.asarray(T_samples, dtype=float)[order]
    cphi_s = cphi[mask][order]
    ca_s = ca[mask][order]
    cv_s = cv[mask][order]

    yk_values = yk.sol(T_samples)
    yphi = np.asarray(yk_values[0], dtype=float)
    yp = np.asarray(yk_values[1], dtype=float)
    yA = A_of_T_phi(T_samples, yphi)
    yA_safe = np.maximum(yA, 1e-14)
    ya = (np.abs(yp) / np.sqrt(yA_safe)) ** (1.0 / 3.0)

    pcos = (ca_s**3) * cv_s
    phi_err = float(np.max(np.abs(yphi - cphi_s)))
    p_err = float(np.max(np.abs(yp - pcos)))
    p_rel = float(np.max(np.abs(yp - pcos) / np.maximum(np.abs(pcos), 1e-12)))
    a_err = float(np.max(np.abs(ya - ca_s)))
    A_pos = float(np.min(yA))
    p_sign_match = int(np.sum(np.sign(yp) == np.sign(pcos)))

    max_phi_err = max(max_phi_err, phi_err)
    max_p_err = max(max_p_err, p_err)
    max_p_rel_err = max(max_p_rel_err, p_rel)
    max_a_err = max(max_a_err, a_err)
    comparison_count += 1

    intervals.append(
        {
            "index": interval_index,
            "t_left": float(tl),
            "t_right": float(tr),
            "T_left": float(T_l),
            "T_right": float(T_r),
            "sigma": float(sigma),
            "min_abs_v": float(min_abs_v),
            "min_A_york_samples": float(A_pos),
            "phi_max_err": phi_err,
            "p_max_err": p_err,
            "p_relative_err": p_rel,
            "p_york_range": [float(np.min(yp)), float(np.max(yp))],
            "p_cosmic_range": [float(np.min(pcos)), float(np.max(pcos))],
            "p_sign_match": p_sign_match,
            "samples": int(np.count_nonzero(mask)),
            "status": int(yk.status),
            "a_max_err": a_err,
        }
    )

record(
    "York- und kosmische Integration stimmen auf regularen Intervallen ueberein",
    comparison_count >= 2
    and max_phi_err < 1e-6
    and max_p_rel_err < 1e-5
    and max_a_err < 5e-4,
    "Intervalle=%d; max Phi-Fehler=%.3e; max relativer p-Fehler=%.3e; max absoluter p-Fehler=%.3e; max a-Fehler=%.3e"
    % (comparison_count, max_phi_err, max_p_rel_err, max_p_err, max_a_err),
)
record(
    "Zeittransformation ist nur auf Intervallen mit Phi_dot != 0 invertierbar",
    np.isfinite(min_abs_v_regular) and min_abs_v_regular > 1e-4,
    "min |Phi_dot| auf verglichenen Regularitaetsintervallen = %.3e" % min_abs_v_regular,
)
record(
    "Keine York-Integration ueber A=0: Intervalle enden vor Ereignissen",
    comparison_count >= 2 and all(iv["min_A_york_samples"] > 0.0 for iv in intervals),
    "min A in verglichenen York-Intervallen = %.3e"
    % (min((iv["min_A_york_samples"] for iv in intervals), default=float("nan"))),
)

all_ok = all(r["pass"] for r in results)
report = {
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "model": {
        "k": 0,
        "Lambda": 0,
        "G": G,
        "kappa": KAPPA,
        "m": m,
        "potential": "0.5*m^2*Phi^2",
        "initial": {"a": a0, "Phi": phi0, "Phi_dot": v0, "H": H0},
        "t_end": t_end,
    },
    "events": {
        "count": len(events),
        "times": events,
        "A_values": event_A_values,
        "v_values": event_v_values,
    },
    "cosmic_checks": {
        "friedmann_max": friedmann_max,
        "A_identity_max": A_identity_max,
        "dT_min": dT_min,
    },
    "comparison_intervals": intervals,
    "max_errors": {
        "phi": max_phi_err,
        "p_absolute": max_p_err,
        "p_relative": max_p_rel_err,
        "a": max_a_err,
    },
    "all_pass": all_ok,
    "checks": results,
}

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)

print("\n==== SUMMARY ====")
print("all_pass =", all_ok)
print("checks   =", len(results))
print("events   =", len(events))
print("json     =", OUT)
raise SystemExit(0 if all_ok else 1)
