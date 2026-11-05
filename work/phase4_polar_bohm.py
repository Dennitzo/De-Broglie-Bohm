"""Phase 4/P2: nachvollziehbare Polarzerlegung des Modells aus S19.7.

G ist die Impuls-Hesse-Matrix von Hbar. Daher lautet der kinetische
Quantenoperator -hbar**2/2 * Delta und Q = -hbar**2/(2*R) * Delta R.
Die Polaridentitaet fuer -hbar**2 * Delta wird davon getrennt geprueft.
Die Familien sind Testfunktionen, keine schon nachgewiesenen Loesungen.
Ausgabe liegt neben diesem Skript; jeder abgeschlossene Teil wird gesichert.
"""

from functools import lru_cache
import json
from pathlib import Path
import time
import traceback

import mpmath as mp
import sympy as sp


OUT = Path(__file__).resolve().with_suffix(".json")
START = time.monotonic()
res = {}


def progress(message, expr=None):
    suffix = "" if expr is None else f"; {sp.count_ops(expr)} symbolische Operationen"
    print(f"[{time.monotonic() - START:9.2f}s] {message}{suffix}", flush=True)


def dump():
    temporary = OUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    temporary.replace(OUT)


def zero_report(expr, label, substitutions=None):
    """Keine globale simplify/expand-Kaskade und kein numerischer Scheinbeweis."""
    result = {"label": label, "symbolisch_null": expr == 0}
    if expr == 0:
        result["numeric_max_rest"] = "0"
    elif substitutions is not None:
        remaining = expr.subs(substitutions)
        if remaining.free_symbols:
            result["offene_symbole"] = sorted(map(str, remaining.free_symbols))
        else:
            result["numeric_rest"] = str(abs(sp.N(remaining, 30)))
            result["hinweis"] = "Eine Zahlenprobe ist kein allgemeiner Nullnachweis."
    return result


def generic_operator_checks():
    """Produktregel unabhaengig an freien Funktionen statt an zwei Kopien
    derselben geschlossenen Formel pruefen. Einzelterme skalieren auf jede
    Dimension; formale Symmetrie setzt geeignete Randbedingungen voraus.
    """
    x, y, h = sp.symbols("x y h", real=True, nonzero=True)
    r = sp.Function("R")(x, y)
    s = sp.Function("S")(x, y)
    density = sp.Function("mu")(x, y)
    coefficient = sp.Function("g")(x, y)
    psi = r * sp.exp(sp.I * s / h)
    checks = []
    for second in (x, y):
        direct = sp.diff(density * coefficient * sp.diff(psi, second), x) / (density * psi)
        drx, dry = sp.diff(r, x) / r, sp.diff(r, second) / r
        dsx, dsy = sp.diff(s, x), sp.diff(s, second)
        drift = sp.diff(coefficient, x) + coefficient * sp.diff(density, x) / density
        product = coefficient * (
            sp.diff(r, x, second) / r
            + sp.I / h * (drx * dsy + dry * dsx + sp.diff(s, x, second))
            - dsx * dsy / h**2
        ) + drift * (dry + sp.I * dsy / h)
        checks.append(zero_report(sp.expand(direct - product), f"Produktregel d_x(mu g d_{second} Psi)/(mu Psi)"))

    f, g = sp.Function("f")(x, y), sp.Function("v")(x, y)
    a, b, d = (sp.Function(name)(x, y) for name in ("a", "b", "d"))
    matrix = ((a, b), (b, d))
    coords = (x, y)

    def delta(value):
        return sum(sp.diff(density * matrix[i][j] * sp.diff(value, coords[j]), coords[i])
                   for i in range(2) for j in range(2)) / density

    boundary = sum(sp.diff(density * matrix[i][j] * (g * sp.diff(f, coords[j]) - f * sp.diff(g, coords[j])), coords[i])
                   for i in range(2) for j in range(2))
    formal = zero_report(sp.expand((delta(f) * g - f * delta(g)) * density - boundary),
                         "Formale Greensche Identitaet fuer symmetrisches G")
    formal["geltungsbereich"] = "Formale Symmetrie nach Randtermkontrolle; kein Nachweis der Selbstadjungiertheit einer Operatordomaene."
    res["selbstadjungiertheit"] = formal
    res["produktregel_allgemein"] = checks
    dump()


def build_model():
    progress("Phase 1: Metrik aus Impuls-Hessian aufbauen")
    a0, a1, ac, c, q = sp.symbols("alpha_0 alpha_1 alpha_c c q", real=True)
    coordinates = (a0, a1, ac, c, q)
    momenta = sp.symbols("P_0 P_1 P_c p_0 p_q", real=True)
    p0, p1, pc, u0, uq = momenta
    beta, mass, hbar = sp.symbols("beta m hbar", positive=True)
    i0, i1, i2 = (sp.besseli(n, 3 * a1) for n in range(3))
    cz = (sp.besseli(0, a1) - sp.besseli(2, a1)) / 2
    geo = -(p0**2 * i0 - 4 * p0 * p1 * i1 + 2 * p1**2 * (i0 + i2) + 2 * pc**2 * (i0 - i2))
    matter = u0**2 * i0 - 4 * u0 * uq * i1 + 2 * uq**2 * (i0 + i2)
    potential = (-sp.exp(a0) * a1**2 * cz / (6 * beta) + sp.exp(a0) * q**2 * cz / 2
                 + mass**2 * sp.exp(3 * a0) * (c**2 * i0 + 2 * c * q * i1 + q**2 * (i0 + i2) / 2) / 2)
    hamiltonian = beta * sp.exp(-3 * a0) * (geo + matter) / 2 + potential
    metric = sp.hessian(hamiltonian, momenta)
    support = [(i, j) for i in range(5) for j in range(5) if metric[i, j] != 0]
    determinant_piece = i0 * (i0 + i2) - 2 * i1**2
    density = beta**(-sp.Rational(5, 2)) * sp.exp(sp.Rational(15, 2) * a0) / (sp.sqrt(8) * sp.sqrt(i0 - i2) * determinant_piece)
    # The positive density branch is used where I0-I2>0 and determinant_piece>0.
    log_density = (sp.Rational(15, 2),
                   -sp.diff(i0 - i2, a1) / (2 * (i0 - i2)) - sp.diff(determinant_piece, a1) / determinant_piece,
                   sp.S.Zero, sp.S.Zero, sp.S.Zero)
    drift = tuple(sum(sp.diff(metric[i, j], coordinates[i]) + log_density[i] * metric[i, j]
                      for i in range(5) if metric[i, j] != 0) for j in range(5))
    res["operator"] = {
        "kinetik": "H_kin = (1/2) G^{AB} p_A p_B; H_hat_kin = -(hbar^2/2) Delta",
        "dichte_zweig": "I0(3 alpha_1)-I2(3 alpha_1)>0 und Dp>0",
        "mu": str(density), "nichtnull_metrikeintraege": len(support),
        "faktor_korrektur": "Q = -(hbar^2/2) Delta R/R; die vorherige Fassung enthielt den doppelten Wert.",
    }
    progress("Phase 1: Metrik fertig", metric)
    dump()
    return dict(coords=coordinates, beta=beta, mass=mass, hbar=hbar, metric=metric,
                support=support, density=density, log_density=log_density, drift=drift, potential=potential)


def symbolic_families(model):
    progress("Phase 2: Gemeinsame logarithmische Amplitudenableitungen vorbereiten")
    coordinates, hbar = model["coords"], model["hbar"]
    a0, a1, ac, c, q = coordinates
    metric, support, drift = model["metric"], model["support"], model["drift"]
    s0, s1, sc, sq, scc = sp.symbols("sigma_0 sigma_1 sigma_c sigma_q sigmac_", positive=True)
    centers = (s0, s1, scc, sc, sq)
    pi0, piq, p0c = sp.symbols("pi_0 pi_q p_0c", real=True)
    log_gauss = -sum((coord - width)**2 / (4 * width**2) for coord, width in zip(coordinates, centers))
    amplitude = model["density"]**(-sp.Rational(1, 2)) * sp.exp(log_gauss)
    r = tuple(-model["log_density"][i] / 2 + sp.diff(log_gauss, coord) for i, coord in enumerate(coordinates))
    r_second = {(i, j): sp.diff(r[j], coordinates[i]) for i, j in support}
    laplace_ratio = sum(metric[i, j] * (r_second[i, j] + r[i] * r[j]) for i, j in support) + sum(drift[j] * r[j] for j in range(5))
    quantum = -hbar**2 * laplace_ratio / 2
    progress("Phase 2: Delta R/R einmalig berechnet", laplace_ratio)
    res["gemeinsames_quantenpotential"] = {"Delta_R_durch_R": str(laplace_ratio), "Q": str(quantum),
                                          "koeffizient_hbar2_unabhaengig_von_hbar": not laplace_ratio.has(hbar)}
    dump()
    laplace_psi_ratios = {}
    for tag, phase in (("A", piq * q), ("C", pi0 * c + piq * q + sp.Rational(3, 2) * p0c * q * a1)):
        progress(f"Phase 2/{tag}: Phase, Strom und Polarzerlegung")
        ds = tuple(sp.diff(phase, coord) for coord in coordinates)
        velocity = tuple(sum(metric[i, j] * ds[j] for j in range(5) if metric[i, j] != 0) for i in range(5))
        phase_laplace = sum(metric[i, j] * sp.diff(ds[j], coordinates[i]) for i, j in support) + sum(drift[j] * ds[j] for j in range(5))
        gradient_square = sp.expand(sum(metric[i, j] * ds[i] * ds[j] for i, j in support))
        continuity_ratio = phase_laplace + 2 * sum(r[i] * velocity[i] for i in range(5))
        logarithmic_psi = tuple(r[i] + sp.I * ds[i] / hbar for i in range(5))
        psi_ratio = sum(metric[i, j] * (sp.diff(logarithmic_psi[j], coordinates[i]) + logarithmic_psi[i] * logarithmic_psi[j])
                        for i, j in support) + sum(drift[j] * logarithmic_psi[j] for j in range(5))
        laplace_psi_ratios[tag] = psi_ratio
        # This checks assembled coefficients; the independent direct differential
        # check is generic_operator_checks(), and the finite-difference route below.
        polar_residual = sp.expand(-hbar**2 * psi_ratio - (gradient_square - hbar**2 * laplace_ratio - sp.I * hbar * continuity_ratio))
        if polar_residual != 0:
            # Cancel only this already small residual (637 ops in family C),
            # never the full Bessel/exponential wavefunction or operator.
            polar_residual = sp.cancel(polar_residual)
        progress(f"Phase 2/{tag}: Polaridentitaet ausgewertet")
        res[f"polaridentitaet_{tag}"] = zero_report(polar_residual, f"Polaridentitaet Familie {tag}, dividiert durch Psi")
        res[f"polaridentitaet_{tag}"]["pruefroute"] = "Logarithmische Ableitungen; allgemeine Produktregel separat differenziert."
        res[f"kontinuitaet_{tag}"] = zero_report(continuity_ratio, f"div(mu J)/(mu R^2), Familie {tag}")
        res[f"kontinuitaet_{tag}"]["status"] = "ungeklaert" if continuity_ratio != 0 else "identisch_null"
        res[f"kontinuitaet_{tag}"]["ausdruck"] = str(continuity_ratio)
        res[f"kontinuitaet_{tag}"]["hinweis"] = "Ein nicht verschwindender Rest bedeutet: Diese zeitunabhaengige Testfunktion ist keine bestaetigte stationaere Loesung."
        cs = -sp.diff(phase, ac) / 3
        cc = sp.expand(sp.diff(phase, a1) / 3 - sp.diff(phase, c) * q / 2)
        res[f"constraints_{tag}"] = {"C_s_auf_Phase": str(cs), "C_c_auf_Phase": str(cc), "C_s_null": cs == 0, "C_c_null": cc == 0}
        if tag == "C":
            res[f"constraints_{tag}"]["erforderliche_relation"] = "p_0c = pi_0 (oder spezieller Punkt q=0); fuer unabhaengige Parameter ist C_c nicht null."
        res[f"fuehrung_{tag}"] = dict(zip(("V_alpha_0", "V_alpha_1", "V_alpha_c", "V_c", "V_q"), map(str, velocity)))
        res[f"fuehrung_{tag}"]["alpha_c_eingefroren"] = velocity[2] == 0
        res[f"gefuehrte_momente_{tag}"] = dict(zip(("P_0", "P_1", "P_c", "p_0", "p_q"), map(str, ds)))
        classical = gradient_square / 2 + model["potential"]
        total = classical + quantum
        progress(f"Phase 2/{tag}: c/q-Polynomkoeffizienten sammeln")
        # Coefficients contain Bessel functions and rational log derivatives.
        # Inferring an algebraic coefficient domain tries to normalise that
        # unrelated special-function algebra. EX keeps coefficients as exact
        # expressions and only collects the requested generators c and q.
        polynomial = sp.Poly(total, c, q, domain="EX")
        res[f"hj_und_qpotential_{tag}"] = {
            "klassischer_hj_rest": str(classical), "quantenpotential_geschlossen": str(quantum),
            "hj_rest_plus_q": str(total), "quantenpotential_ohbar2_koeffizient": str(-laplace_ratio / 2),
            "koeffizient_hbar2_unabhaengig_von_hbar": not laplace_ratio.has(hbar),
            "polynomkoeffizienten_c_q": {f"c^{powers[0]} q^{powers[1]}": str(coefficient) for powers, coefficient in polynomial.terms()},
            "loesungsstatus": "Testfamilie; HJ- und Kontinuitaetsreste muessen verschwinden, bevor sie als Loesung gelten kann.",
        }
        dump()
        progress(f"Phase 2/{tag}: Teilresultate gesichert")
    res["gemischteranteil"] = {"G_cq": str(metric[3, 4]),
        "bemerkung": "G^{cq} erzeugt gemischte c/q-Ableitungen. alpha_1-abhaengige Koeffizienten allein sind keine gemischte alpha_1/q-Zweitableitung; G^{alpha_1,q}=0."}
    dump()
    return dict(amplitude=amplitude, laplace_ratio=laplace_ratio, psi_ratios=laplace_psi_ratios,
                widths=(s0, s1, sc, sq, scc), piq=piq, pi0=pi0, p0c=p0c)


@lru_cache(maxsize=256)
def mp_model(a0v, a1v, betav):
    a0v, a1v, betav = map(mp.mpf, (a0v, a1v, betav))
    i0, i1, i2 = (mp.besseli(n, 3 * a1v) for n in range(3))
    d = i0 * (i0 + i2) - 2 * i1**2
    pre = betav * mp.exp(-3 * a0v)
    density = betav**(-mp.mpf("2.5")) * mp.exp(mp.mpf("7.5") * a0v) / (mp.sqrt(8) * mp.sqrt(i0 - i2) * d)
    matrix = mp.zeros(5, 5)
    matrix[0, 0], matrix[0, 1], matrix[1, 0] = -pre * i0, 2 * pre * i1, 2 * pre * i1
    matrix[1, 1], matrix[2, 2] = -2 * pre * (i0 + i2), -2 * pre * (i0 - i2)
    matrix[3, 3], matrix[3, 4], matrix[4, 3] = pre * i0, -2 * pre * i1, -2 * pre * i1
    matrix[4, 4] = 2 * pre * (i0 + i2)
    return density, matrix


def lb_numeric(function, point, step, beta_value):
    """Originaler Divergenzoperator aus zentralen Differenzen, ohne log-R-Formel."""
    density, _ = mp_model(point[0], point[1], beta_value)

    def flux(index, at):
        local_density, matrix = mp_model(at[0], at[1], beta_value)
        result = mp.mpf(0)
        for j in range(5):
            if matrix[index, j] == 0:
                continue
            plus, minus = list(at), list(at)
            plus[j] += step
            minus[j] -= step
            # Functions take one coordinate vector, not five positional arguments.
            result += local_density * matrix[index, j] * (function(plus) - function(minus)) / (2 * step)
        return result

    result = mp.mpf(0)
    for i in range(5):
        plus, minus = list(point), list(point)
        plus[i] += step
        minus[i] -= step
        result += (flux(i, plus) - flux(i, minus)) / (2 * step)
    return result / density


def numeric_checks(model, families):
    progress("Phase 3: Unabhaengige numerische Divergenzroute und Schrittweitenvergleich")
    mp.mp.dps = 40
    beta_value, hbar_value, phase_value = map(mp.mpf, ("1.7", "0.13", "0.7"))
    width_values = tuple(map(mp.mpf, ("0.3", "0.2", "0.25", "0.35", "0.3")))
    points = [tuple(map(mp.mpf, values)) for values in (("0.3", "0.25", "-0.15", "0.6", "0.4"), ("-0.2", "0.1", "0.05", "-0.3", "0.55"))]

    def amplitude(point):
        density, _ = mp_model(point[0], point[1], beta_value)
        centers = (width_values[0], width_values[1], width_values[4], width_values[2], width_values[3])
        return mp.exp(-sum((coordinate - width)**2 / (4 * width**2) for coordinate, width in zip(point, centers))) / mp.sqrt(density)

    def psi(point, family="A"):
        phase = phase_value * point[4]
        if family == "C":
            phase += mp.mpf("0.4") * point[3] + mp.mpf("1.5") * mp.mpf("0.2") * point[4] * point[1]
        return amplitude(point) * mp.exp(mp.j * phase / hbar_value)

    rows, comparison, metric_checks = [], [], []
    for index, point in enumerate(points):
        progress(f"Phase 3/{index + 1}: Symbolische Koeffizienten numerisch auswerten")
        substitutions = dict(zip(model["coords"], map(lambda value: sp.Float(str(value), 40), point)))
        substitutions.update(dict(zip(families["widths"], map(lambda value: sp.Float(str(value), 40), width_values))))
        substitutions.update({model["beta"]: sp.Float(str(beta_value), 40), model["mass"]: sp.Float("0.9", 40),
                              model["hbar"]: sp.Float(str(hbar_value), 40), families["piq"]: sp.Float(str(phase_value), 40),
                              families["pi0"]: sp.Float("0.4", 40), families["p0c"]: sp.Float("0.2", 40)})

        def value(expression):
            number = sp.N(expression.subs(substitutions), 35)
            return mp.mpc(str(sp.re(number)), str(sp.im(number)))

        symbolic_r = amplitude(point) * value(families["laplace_ratio"])
        symbolic_psi = psi(point) * value(families["psi_ratios"]["A"])
        symbolic_psi_c = psi(point, "C") * value(families["psi_ratios"]["C"])
        density, independent_metric = mp_model(point[0], point[1], beta_value)
        metric_error = max(abs(value(model["metric"][i, j]) - independent_metric[i, j]) for i in range(5) for j in range(5))
        density_error = abs(density**-2 - abs(mp.det(independent_metric)))
        metric_checks.append({"punkt": index, "metric_max_rest": str(metric_error), "density_det_rest": str(density_error)})
        previous_errors = None
        for step in map(mp.mpf, ("0.05", "0.025", "0.0125")):
            progress(f"Phase 3/{index + 1}: Direkter Divergenzoperator, Schrittweite {step}")
            direct_r = lb_numeric(amplitude, list(point), step, beta_value)
            direct_psi = lb_numeric(psi, list(point), step, beta_value)
            direct_psi_c = lb_numeric(lambda at: psi(at, "C"), list(point), step, beta_value)
            errors = (abs(symbolic_r - direct_r), abs(symbolic_psi - direct_psi), abs(symbolic_psi_c - direct_psi_c))
            rows.append({"punkt": index, "h": str(step), "Delta_R_direkt": str(direct_r),
                         "Delta_Psi_direkt": str(direct_psi), "Delta_Psi_C_direkt": str(direct_psi_c)})
            row = {"punkt": index, "h": str(step), "Delta_R_symbolisch": str(symbolic_r), "Delta_R_direkt": str(direct_r),
                   "Delta_Psi_symbolisch": str(symbolic_psi), "Delta_Psi_direkt": str(direct_psi),
                   "Delta_Psi_C_symbolisch": str(symbolic_psi_c), "Delta_Psi_C_direkt": str(direct_psi_c),
                   "abs_fehler_R": str(errors[0]), "abs_fehler_Psi": str(errors[1]), "abs_fehler_Psi_C": str(errors[2])}
            if previous_errors is not None:
                row["fehlerverhaeltnis_R"] = str(previous_errors[0] / errors[0]) if errors[0] else "null"
                row["fehlerverhaeltnis_Psi"] = str(previous_errors[1] / errors[1]) if errors[1] else "null"
                row["fehlerverhaeltnis_Psi_C"] = str(previous_errors[2] / errors[2]) if errors[2] else "null"
                row["fehler_nimmt_ab"] = all(current < previous for current, previous in zip(errors, previous_errors))
            previous_errors = errors
            comparison.append(row)
            res["numerik_operator"], res["vergleich_symbolisch_direkt"] = rows, comparison
            res["metrische_gegenkontrolle"] = metric_checks
            dump()
    res["numerische_pruefgrenzen"] = "Zwei lokale Punkte und drei Schrittweiten pruefen Implementierung/Konvergenz, nicht globale physikalische Gueltigkeit oder Loesungseigenschaft."
    dump()


def main():
    res["laufstatus"] = "in_progress"
    progress("Phase 0: Allgemeine Produktregel und Greensche Identitaet pruefen")
    generic_operator_checks()
    model = build_model()
    families = symbolic_families(model)
    numeric_checks(model, families)
    res["laufstatus"] = "completed_with_open_scientific_questions"
    res["offene_grenzen"] = "Formale Identitaeten und lokale Implementierungschecks sind keine bestaetigte Loesung der stationaeren HJ-/Kontinuitaetsgleichungen. Nicht verschwindende Reste bleiben offen."
    dump()
    progress(f"Fertig: Teilresultate und offene Grenzen in {OUT}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        res["laufstatus"] = "failed_or_interrupted"
        res["fehler"] = f"{type(error).__name__}: {error}"
        dump()
        traceback.print_exc()
        raise
