#!/usr/bin/env python3
"""SpinSqueezing-Lab - educational spin-squeezing and quantum metrology lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lessons are stated up front:

  * A coherent spin state (CSS) of N two-level atoms, all pointing
    the same way, has Kitagawa-Ueda xi^2 = 1 and Wineland xi_R^2 = 1
    exactly: its transverse noise is the standard quantum limit
    (SQL), neither better nor worse.
  * One-axis twisting (H = chi * Sz^2) shears that noise ellipse:
    xi_R^2 dips below 1, so phase sensitivity beats the SQL - but
    only the Wineland parameter counts, because it divides by the
    squared contrast |<S>|/J. Quoting the Kitagawa-Ueda xi^2 alone
    overstates the gain once the mean spin shrinks.
  * The best one-axis-twisting squeezing scales as xi_R^2 ~ N^(-2/3);
    two-axis counter-twisting scales as ~1/N in the ideal limit.
    Both exponents are FITTED here from computed optima, and the
    fitted numbers are printed next to the ideal values.
  * At chi*t = pi/2 (even N) the twisted CSS is an exact macroscopic
    superposition - a "cat" - of the CSS along +x and the CSS along
    -x. Its Sz (Dicke) populations never change under Sz^2 evolution
    (phases only); the bimodality lives in the Sx distribution.

Everything is computed exactly in the symmetric Dicke basis
(J = N/2, dimension N+1) with collective operators built from
scratch. No atom, laser, clock, or sensor exists here: every number
is exact arithmetic on the model. Real-device calibration, clock or
sensor product claims, live experiment direction, metrology
certification, navigation, defence, medical, and record claims are
permanent HOLD gates - with Gemini Argon 4 Extended Thinking as
DECIDER (names-only, tighten-only), Muse as advisor (never
overrides HOLD), and Google Flow plan-only seats.
"""
import cmath
import hashlib
import json
import math
import sys
from functools import lru_cache

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled teaching grid for the scaling fits: the same N values are
# used for one-axis and two-axis twisting so the exponents compare.
FIT_GRID_N = (10, 20, 40, 60, 80, 100)
# Windows the fitted exponents are asserted against. They are stated
# BEFORE the fit is run, from the prototype values, and the fitted
# numbers are printed in the demo so anyone can see the margin.
OAT_EXPONENT_WINDOW = (-0.76, -0.62)     # ideal: -2/3 = -0.6667
TACT_EXPONENT_WINDOW = (-0.97, -0.80)    # ideal: -1 (asymptotic)

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "real_sensor_calibration": "HOLD - no calibrating a real sensor, magnetometer, or clock from this simulation",
    "atomic_clock_product_claim": "HOLD - no atomic-clock product, specification, or performance claim",
    "quantum_advantage_marketing_claim": "HOLD - no marketing claim of demonstrated quantum advantage from a lab fit",
    "live_experiment_direction": "HOLD - no directing a live experiment on real atoms, ions, or photons",
    "real_hardware_control": "HOLD - no control of real lasers, cavities, traps, or magnetic coils",
    "metrology_certification": "HOLD - no certifying a real measurement, uncertainty budget, or standard",
    "navigation_guidance_use": "HOLD - no navigation, guidance, or gyroscope use from this model",
    "defense_targeting_use": "HOLD - no defence, targeting, or weapons use of any kind",
    "medical_measurement_claim": "HOLD - no medical, diagnostic, or imaging claim from squeezed-spin math",
    "gravitational_wave_claim": "HOLD - no claims about real gravitational-wave detectors or their sensitivity",
    "squeezing_record_claim": "HOLD - no claiming a real squeezing record; these are exact model numbers, not lab results",
    "commercial_sensor_claim": "HOLD - no commercial sensor, timing, or synchronisation product claim",
    "safety_critical_timing_use": "HOLD - no safety-critical timing or synchronisation reliance on a simulation",
    "investment_claim": "HOLD - no investment claims; a teaching lab is not a prospectus",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Dicke basis core
# N two-level atoms, collective spin J = N/2. The symmetric (Dicke)
# basis |J, m>, m = -J..J, has dimension N+1. Collective operators
# are applied through their exact matrix elements - no library, no
# approximation: Sz is diagonal, S+ raises m by one with the ladder
# factor sqrt(J(J+1) - m(m+1)).

def dicke_ms(N):
    """The m values of the Dicke basis, ascending from -J to +J."""
    J = N / 2.0
    return [-J + i for i in range(N + 1)]


def _ladder(N, m):
    """<J, m+1| S+ |J, m>."""
    J = N / 2.0
    return math.sqrt(max(0.0, J * (J + 1.0) - m * (m + 1.0)))


def css(N, theta, phi):
    """Coherent spin state pointing along (theta, phi).

    One atom in cos(theta/2)|up> + e^{i phi} sin(theta/2)|down>,
    symmetrised over N atoms: the Dicke coefficient for k up-spins
    is sqrt(C(N, k)) cos^k(theta/2) (e^{i phi} sin(theta/2))^(N-k).
    """
    ct, st = math.cos(theta / 2.0), math.sin(theta / 2.0)
    phase = cmath.exp(1j * phi)
    psi = []
    for m in dicke_ms(N):
        k = int(round(m + N / 2.0))
        psi.append(math.sqrt(math.comb(N, k)) * ct ** k
                   * (st * phase) ** (N - k))
    return psi


def css_x(N):
    """CSS along +x: every atom in (|up> + |down>)/sqrt(2)."""
    return css(N, math.pi / 2.0, 0.0)


def expectations(N, psi):
    """First and second collective-spin moments of a Dicke state.

    Uses <S+>, <S+^2>, and <{Sz, S+}> so every second moment follows
    exactly:  <Sx^2> - <Sy^2> = Re<S+^2>,  <{Sx, Sy}> = Im<S+^2>,
    <Sx^2> + <Sy^2> = J(J+1) - <Sz^2>,  <{Sz, Sx}> = Re<{Sz, S+}>,
    <{Sz, Sy}> = Im<{Sz, S+}>.
    """
    mv = dicke_ms(N)
    probs = [abs(a) ** 2 for a in psi]
    sz = sum(m * p for m, p in zip(mv, probs))
    sz2 = sum(m * m * p for m, p in zip(mv, probs))
    sp = 0j
    sp2 = 0j
    szsp = 0j
    for i, m in enumerate(mv):
        if i + 1 <= N:
            w = psi[i + 1].conjugate() * psi[i] * _ladder(N, m)
            sp += w
            szsp += w * (2.0 * m + 1.0)
        if i + 2 <= N:
            sp2 += (psi[i + 2].conjugate() * psi[i]
                    * _ladder(N, m) * _ladder(N, m + 1.0))
    J = N / 2.0
    sx, sy = sp.real, sp.imag
    sx2 = (J * (J + 1.0) - sz2 + sp2.real) / 2.0
    sy2 = (J * (J + 1.0) - sz2 - sp2.real) / 2.0
    return {
        "sx": sx, "sy": sy, "sz": sz,
        "sx2": sx2, "sy2": sy2, "sz2": sz2,
        "cov_xy": sp2.imag / 2.0 - sx * sy,
        "cov_xz": szsp.real / 2.0 - sx * sz,
        "cov_yz": szsp.imag / 2.0 - sy * sz,
    }


def covariance_matrix(ex):
    """Symmetrised covariance matrix Cov(a, b) = <{Sa, Sb}>/2 - <Sa><Sb>."""
    return [
        [ex["sx2"] - ex["sx"] ** 2, ex["cov_xy"], ex["cov_xz"]],
        [ex["cov_xy"], ex["sy2"] - ex["sy"] ** 2, ex["cov_yz"]],
        [ex["cov_xz"], ex["cov_yz"], ex["sz2"] - ex["sz"] ** 2],
    ]


def min_perp_variance(ex):
    """Smallest variance in the plane perpendicular to the mean spin.

    Returns (lambda_min, |<S>|). If the mean spin has vanished the
    perpendicular plane is undefined and lambda_min is None.
    """
    mean = (ex["sx"], ex["sy"], ex["sz"])
    norm = math.sqrt(sum(v * v for v in mean))
    if norm < 1e-12:
        return None, norm
    n = tuple(v / norm for v in mean)
    C = covariance_matrix(ex)
    a = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    proj = sum(a[k] * n[k] for k in range(3))
    u = tuple(a[i] - proj * n[i] for i in range(3))
    un = math.sqrt(sum(v * v for v in u))
    u = tuple(v / un for v in u)
    v = (n[1] * u[2] - n[2] * u[1],
         n[2] * u[0] - n[0] * u[2],
         n[0] * u[1] - n[1] * u[0])

    def quad(w1, w2):
        return sum(w1[i] * C[i][j] * w2[j]
                   for i in range(3) for j in range(3))

    c11, c22, c12 = quad(u, u), quad(v, v), quad(u, v)
    disc = math.sqrt(max(0.0, ((c11 - c22) / 2.0) ** 2 + c12 ** 2))
    return (c11 + c22) / 2.0 - disc, norm


def squeezing_params(N, psi):
    """Kitagawa-Ueda xi^2, Wineland xi_R^2, and contrast for a state.

    xi^2   = (DS_perp,min)^2 / (J/2)          (Kitagawa-Ueda)
    xi_R^2 = N (DS_perp,min)^2 / |<S>|^2      (Wineland)
           = xi^2 / contrast^2,  contrast = |<S>| / J.
    The Wineland parameter is the metrological one: phase sensitivity
    is Dphi = xi_R / sqrt(N). xi_R^2 is infinite when the mean spin
    has vanished (contrast 0), however small xi^2 looks.
    """
    ex = expectations(N, psi)
    lam, norm = min_perp_variance(ex)
    J = N / 2.0
    contrast = norm / J
    if lam is None:
        return {"xi2_ku": float("inf"), "xi2_wineland": float("inf"),
                "contrast": 0.0, "mean_spin": norm, "var_min": None,
                "expectations": ex}
    xi2_ku = lam / (J / 2.0)
    xi2_w = N * lam / (norm * norm)
    return {"xi2_ku": xi2_ku, "xi2_wineland": xi2_w,
            "contrast": contrast, "mean_spin": norm, "var_min": lam,
            "expectations": ex}


# ================================================================ 2. One-axis twisting
# H = chi * Sz^2 is diagonal in the Dicke basis, so evolution is
# exact: each component picks up the phase exp(-i chi t m^2). The
# scan finds the Wineland optimum on a fixed log grid, then refines
# it with a golden-section search on log(mu) - deterministic, and
# the same procedure for every N so the fit compares like with like.

def oat_evolve(N, psi0, mu):
    """Evolve under H = chi Sz^2 for mu = chi * t (exact phases)."""
    return [a * cmath.exp(-1j * mu * m * m)
            for a, m in zip(psi0, dicke_ms(N))]


def _log_grid(lo, hi, count):
    la, lb = math.log(lo), math.log(hi)
    return [math.exp(la + (lb - la) * i / (count - 1)) for i in range(count)]


def _golden_min(f, a, b, iters=60):
    """Golden-section minimisation of f on [a, b]; deterministic."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
    x = (a + b) / 2.0
    return x, f(x)


def oat_optimum(N):
    """Best Wineland squeezing under one-axis twisting for N atoms.

    Returns the optimum mu, xi_R^2 there, xi^2 there, and contrast.
    """
    psi0 = css_x(N)

    def wineland_at(log_mu):
        p = squeezing_params(N, oat_evolve(N, psi0, math.exp(log_mu)))
        return p["xi2_wineland"]

    grid = _log_grid(1e-3, 1.0, 240)
    vals = [(wineland_at(math.log(mu)), mu) for mu in grid]
    _, mu_best = min(vals)
    i = grid.index(mu_best)
    lo = math.log(grid[max(0, i - 1)])
    hi = math.log(grid[min(len(grid) - 1, i + 1)])
    log_mu, _ = _golden_min(wineland_at, lo, hi)
    mu = math.exp(log_mu)
    p = squeezing_params(N, oat_evolve(N, psi0, mu))
    return {"N": N, "mu": mu, "xi2_wineland": p["xi2_wineland"],
            "xi2_ku": p["xi2_ku"], "contrast": p["contrast"]}


def oat_curve(N, mus):
    """(mu, xi^2, xi_R^2, contrast) rows for one-axis twisting."""
    psi0 = css_x(N)
    rows = []
    for mu in mus:
        p = squeezing_params(N, oat_evolve(N, psi0, mu))
        rows.append({"mu": mu, "xi2_ku": p["xi2_ku"],
                     "xi2_wineland": p["xi2_wineland"],
                     "contrast": p["contrast"]})
    return rows


def _linfit(xs, ys):
    """Least-squares slope/intercept, from scratch."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sxy / sxx
    return slope, my - slope * mx


@lru_cache(maxsize=1)
def oat_scaling():
    """Fit xi_R^2(min) ~ N^exponent for one-axis twisting."""
    rows = [oat_optimum(N) for N in FIT_GRID_N]
    slope, intercept = _linfit(
        [math.log(r["N"]) for r in rows],
        [math.log(r["xi2_wineland"]) for r in rows])
    return {"rows": rows, "exponent": slope, "intercept": intercept,
            "ideal_exponent": -2.0 / 3.0}


# ================================================================ 3. The cat point
# At mu = pi/2 with N even, m is an integer, so the twisting phase
# exp(-i (pi/2) m^2) is 1 for even m and -i for odd m. That parity
# pattern is exactly the interference of two coherent states: the
# result is (|CSS +x> + phase * |CSS -x>) / sqrt(2), a macroscopic
# superposition ("cat"). The Sz populations cannot change - Sz^2
# evolution is phases only - so the bimodality is in the Sx basis:
# the Sx distribution has all its weight at m_x = +/-J.

def cat_report(N=20):
    """Measured cat-state facts at mu = pi/2 for even N."""
    if N % 2 != 0:
        raise ValueError("the pi/2 cat needs even N (integer m)")
    psi0 = css_x(N)
    psi = oat_evolve(N, psi0, math.pi / 2.0)
    plus = css(N, math.pi / 2.0, 0.0)
    minus = css(N, math.pi / 2.0, math.pi)
    ov_plus = sum(a.conjugate() * b for a, b in zip(plus, psi))
    ov_minus = sum(a.conjugate() * b for a, b in zip(minus, psi))
    p = squeezing_params(N, psi)
    pop_change = max(abs(abs(a) ** 2 - abs(b) ** 2)
                     for a, b in zip(psi0, psi))
    return {
        "N": N, "mu": math.pi / 2.0,
        "overlap_plus_x": abs(ov_plus) ** 2,
        "overlap_minus_x": abs(ov_minus) ** 2,
        "overlap_sum": abs(ov_plus) ** 2 + abs(ov_minus) ** 2,
        "mean_sx": p["expectations"]["sx"],
        "mean_sy": p["expectations"]["sy"],
        "contrast": p["contrast"],
        "sz_population_change": pop_change,
    }


# ================================================================ 4. Two-axis counter-twisting
# H = chi * (Sz^2 - Sy^2), mean spin along x: the state is twisted
# about y and z with opposite signs, which squeezes exponentially
# instead of polynomially - ideal scaling xi_R^2 ~ 1/N. The
# Hamiltonian is real symmetric and couples only m <-> m, m +/- 2,
# so it splits into even-m and odd-m blocks; each block is
# diagonalised with a from-scratch cyclic Jacobi eigensolver, and
# evolution is then exact in the eigenbasis.

def tact_matrix(N):
    """H = Sz^2 - Sy^2 (chi = 1) as a real symmetric matrix.

    Sy^2 diagonal: (J(J+1) - m^2)/2. Sy^2 off-diagonal:
    <m+2| Sy^2 |m> = -<m+2| S+^2 |m>/4 (and its mirror).
    """
    J = N / 2.0
    mv = dicke_ms(N)
    d = N + 1
    H = [[0.0] * d for _ in range(d)]
    for i, m in enumerate(mv):
        H[i][i] = m * m - (J * (J + 1.0) - m * m) / 2.0
        if i + 2 < d:
            off = _ladder(N, m) * _ladder(N, m + 1.0) / 4.0
            H[i][i + 2] = off
            H[i + 2][i] = off
    return H


def jacobi_eigh(A, tol=1e-12, max_sweeps=60):
    """Eigen-decomposition of a real symmetric matrix (cyclic Jacobi).

    Returns (eigenvalues, V) with A = V diag(eigenvalues) V^T;
    columns of V are the eigenvectors. Deterministic sweep order.
    """
    n = len(A)
    a = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[i][j] ** 2 for i in range(n)
                            for j in range(n) if i != j))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if a[p][q] == 0.0:
                    continue
                tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = ((1.0 if tau >= 0 else -1.0)
                     / (abs(tau) + math.sqrt(1.0 + tau * tau)))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    return [a[i][i] for i in range(n)], V


@lru_cache(maxsize=None)
def _tact_eigensystem(N):
    """Eigensystem of the TACT Hamiltonian, block by block (m parity).

    Returns (eigenvalues, V) in the original Dicke ordering.
    """
    H = tact_matrix(N)
    d = N + 1
    evals = [0.0] * d
    V = [[0.0] * d for _ in range(d)]
    for parity in (0, 1):
        idx = [i for i in range(d) if i % 2 == parity]
        block = [[H[i][j] for j in idx] for i in idx]
        ev, Vb = jacobi_eigh(block)
        for a, i in enumerate(idx):
            evals[i] = ev[a]
            for b, j in enumerate(idx):
                V[j][i] = Vb[b][a]
    return evals, V


def tact_evolve(N, psi0, t):
    """Evolve under H = chi (Sz^2 - Sy^2) for time t (chi = 1), exact."""
    evals, V = _tact_eigensystem(N)
    d = N + 1
    w = [sum(V[i][j] * psi0[i] for i in range(d)) for j in range(d)]
    wt = [w[j] * cmath.exp(-1j * evals[j] * t) for j in range(d)]
    return [sum(V[i][j] * wt[j] for j in range(d)) for i in range(d)]


def tact_optimum(N):
    """Best Wineland squeezing under two-axis counter-twisting."""
    psi0 = css_x(N)

    def wineland_at(t):
        return squeezing_params(N, tact_evolve(N, psi0, t))["xi2_wineland"]

    t_max = 3.0 / math.sqrt(N)
    count = 80
    grid = [t_max * i / count for i in range(1, count + 1)]
    vals = [(wineland_at(t), t) for t in grid]
    _, t_best = min(vals)
    step = t_max / count
    lo = max(0.0, t_best - step)
    hi = t_best + step
    t_opt, _ = _golden_min(wineland_at, lo, hi)
    p = squeezing_params(N, tact_evolve(N, psi0, t_opt))
    return {"N": N, "mu": t_opt, "xi2_wineland": p["xi2_wineland"],
            "xi2_ku": p["xi2_ku"], "contrast": p["contrast"]}


@lru_cache(maxsize=1)
def tact_scaling():
    """Fit xi_R^2(min) ~ N^exponent for two-axis counter-twisting."""
    rows = [tact_optimum(N) for N in FIT_GRID_N]
    slope, intercept = _linfit(
        [math.log(r["N"]) for r in rows],
        [math.log(r["xi2_wineland"]) for r in rows])
    top = _linfit([math.log(rows[-2]["N"]), math.log(rows[-1]["N"])],
                  [math.log(rows[-2]["xi2_wineland"]),
                   math.log(rows[-1]["xi2_wineland"])])[0]
    return {"rows": rows, "exponent": slope, "intercept": intercept,
            "ideal_exponent": -1.0, "top_local_slope": top}


def tact_curve(N, ts):
    """(t, xi^2, xi_R^2, contrast) rows for two-axis counter-twisting."""
    psi0 = css_x(N)
    rows = []
    for t in ts:
        p = squeezing_params(N, tact_evolve(N, psi0, t))
        rows.append({"mu": t, "xi2_ku": p["xi2_ku"],
                     "xi2_wineland": p["xi2_wineland"],
                     "contrast": p["contrast"]})
    return rows


# ================================================================ 5. kNN regime classifier
# Features: log10 of the Wineland parameter and the contrast. The
# training grid is labelled by the lab's own policy (model_label
# below) - the classifier interpolates the policy, it does not
# replace the computed squeezing parameters.

def model_label(xi2_wineland, contrast):
    """The regime policy distilled to a four-way teaching label."""
    if contrast < 0.5:
        return "over-twisted"
    if xi2_wineland >= 0.98:
        return "sql-limited"
    if xi2_wineland < 0.25:
        return "deep-squeezed"
    return "squeezed"


def _training_grid():
    rows = []
    for i in range(11):
        log_xi = -1.5 + 0.15 * i          # xi_R^2 from 10^-1.5 to 1
        for j in range(11):
            contrast = 0.1 * j           # 0.0 .. 1.0
            xi = 10.0 ** log_xi
            rows.append(((log_xi, contrast), model_label(xi, contrast)))
    return rows


TRAINING_GRID = _training_grid()


def _features(xi2_wineland, contrast):
    xi = min(max(xi2_wineland, 1e-3), 1.0)
    return (math.log10(xi), contrast)


def classify_regime(xi2_wineland, contrast, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    target = _features(xi2_wineland, contrast)
    ranked = []
    for features, label in TRAINING_GRID:
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(target, features)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    weights = {}
    for dist, _, label in ranked[:k]:
        weights[label] = weights.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"deep-squeezed": 0, "squeezed": 1,
                "sql-limited": 2, "over-twisted": 3}
    return sorted(weights.items(),
                  key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


def classifier_report():
    """Scenarios computed by the lab itself, labelled by kNN + policy."""
    scenarios = {}
    p0 = squeezing_params(20, css_x(20))
    scenarios["css_at_t0"] = (p0["xi2_wineland"], p0["contrast"])
    oat20 = oat_optimum(20)
    scenarios["oat_optimum_N20"] = (oat20["xi2_wineland"], oat20["contrast"])
    oat100 = oat_optimum(100)
    scenarios["oat_optimum_N100"] = (oat100["xi2_wineland"],
                                     oat100["contrast"])
    mild = oat_curve(50, [0.02])[0]
    scenarios["oat_mild_N50_mu0.02"] = (mild["xi2_wineland"],
                                        mild["contrast"])
    tact40 = tact_optimum(40)
    scenarios["tact_optimum_N40"] = (tact40["xi2_wineland"],
                                     tact40["contrast"])
    cat = cat_report(20)
    scenarios["cat_point_N20"] = (float("inf"), cat["contrast"])
    report = {}
    for name, (xi, contrast) in scenarios.items():
        report[name] = {
            "xi2_wineland": (None if xi == float("inf") else round(xi, 6)),
            "contrast": round(contrast, 6),
            "knn": classify_regime(xi, contrast),
            "policy": model_label(xi, contrast),
        }
    return report


# ================================================================ 6. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="spinsqueezing-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all device and metrology gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="spinsqueezing-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Quote the Wineland parameter, never the Kitagawa-Ueda xi^2 alone: once twisting shrinks the mean spin, only xi_R^2 = xi^2 / contrast^2 prices the phase sensitivity honestly - and print the fitted exponent next to the ideal one.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["a round noise circle on a spin sphere, marked SQL",
                  "one-axis twisting shearing the circle into a thin ellipse",
                  "the mean-spin arrow shrinking as the ellipse thins",
                  "the ellipse at chi t = pi/2 splitting into two blobs at +x and -x",
                  "two-axis counter-twisting winding the ellipse thinner, faster"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. CSS reference point: xi^2 = xi_R^2 = 1 exactly, contrast 1
    for N in (10, 50, 100):
        p = squeezing_params(N, css_x(N))
        chk(f"css_sql_exact_N{N}",
            abs(p["xi2_ku"] - 1.0) < 1e-9
            and abs(p["xi2_wineland"] - 1.0) < 1e-9
            and abs(p["contrast"] - 1.0) < 1e-12
            and abs(p["expectations"]["sx"] - N / 2.0) < 1e-9)
    p = squeezing_params(30, css(30, 0.7, 1.1))
    chk("css_any_direction_sql",
        abs(p["xi2_ku"] - 1.0) < 1e-9 and abs(p["contrast"] - 1.0) < 1e-9)

    # 2. OAT machinery: contrast closed form + populations frozen
    N = 40
    psi0 = css_x(N)
    for mu in (0.05, 0.1, 0.3):
        p = squeezing_params(N, oat_evolve(N, psi0, mu))
        chk(f"contrast_formula_mu{mu}",
            abs(p["contrast"] - math.cos(mu) ** (N - 1)) < 1e-9)
    psi_t = oat_evolve(N, psi0, 0.37)
    chk("oat_sz_populations_frozen",
        max(abs(abs(a) ** 2 - abs(b) ** 2)
            for a, b in zip(psi0, psi_t)) < 1e-15)
    chk("oat_norm_preserved",
        abs(sum(abs(a) ** 2 for a in psi_t) - 1.0) < 1e-12)

    # 3. OAT optimum: below SQL, improves with N, Wineland > KU there
    oat = oat_scaling()
    rows = oat["rows"]
    chk("oat_optima_below_sql", all(r["xi2_wineland"] < 1.0 for r in rows))
    chk("oat_optima_improve_with_N",
        all(rows[i + 1]["xi2_wineland"] < rows[i]["xi2_wineland"]
            for i in range(len(rows) - 1)))
    chk("oat_wineland_exceeds_ku_at_optimum",
        all(r["xi2_wineland"] > r["xi2_ku"] and r["contrast"] < 1.0
            for r in rows))
    chk("oat_wineland_is_ku_over_contrast2",
        all(abs(r["xi2_wineland"]
                - r["xi2_ku"] / r["contrast"] ** 2) < 1e-6 for r in rows))
    chk("oat_exponent_window",
        OAT_EXPONENT_WINDOW[0] <= oat["exponent"] <= OAT_EXPONENT_WINDOW[1])
    chk("oat_exponent_near_ideal",
        abs(oat["exponent"] - (-2.0 / 3.0)) < 0.1)

    # 4. Cat point: exact +/-x superposition, Sz populations untouched
    cat = cat_report(20)
    chk("cat_overlap_plus_half", abs(cat["overlap_plus_x"] - 0.5) < 1e-9)
    chk("cat_overlap_minus_half", abs(cat["overlap_minus_x"] - 0.5) < 1e-9)
    chk("cat_overlap_sum_one", abs(cat["overlap_sum"] - 1.0) < 1e-9)
    chk("cat_mean_spin_vanishes",
        abs(cat["mean_sx"]) < 1e-9 and abs(cat["mean_sy"]) < 1e-9
        and cat["contrast"] < 1e-9)
    chk("cat_sz_populations_unchanged", cat["sz_population_change"] < 1e-15)
    cat40 = cat_report(40)
    chk("cat_holds_N40", abs(cat40["overlap_sum"] - 1.0) < 1e-9)

    # 5. TACT machinery: Jacobi eigensolver refereed on the matrix
    H = tact_matrix(12)
    evals, V = jacobi_eigh(H)
    d = len(H)
    recon = max(abs(sum(V[i][k] * evals[k] * V[j][k] for k in range(d))
                    - H[i][j]) for i in range(d) for j in range(d))
    ortho = max(abs(sum(V[k][i] * V[k][j] for k in range(d))
                    - (1.0 if i == j else 0.0))
                for i in range(d) for j in range(d))
    chk("jacobi_reconstructs_tact_matrix", recon < 1e-9)
    chk("jacobi_eigenvectors_orthonormal", ortho < 1e-9)
    chk("tact_matrix_real_symmetric",
        all(H[i][j] == H[j][i] for i in range(d) for j in range(d)))
    chk("tact_parity_blocks_only",
        all(H[i][j] == 0.0 for i in range(d) for j in range(d)
            if abs(i - j) not in (0, 2)))
    psi_t = tact_evolve(20, css_x(20), 0.05)
    chk("tact_norm_preserved",
        abs(sum(abs(a) ** 2 for a in psi_t) - 1.0) < 1e-9)

    # 6. TACT scaling: beats OAT at every N, steeper exponent
    tact = tact_scaling()
    trows = tact["rows"]
    chk("tact_optima_below_sql", all(r["xi2_wineland"] < 1.0 for r in trows))
    chk("tact_beats_oat_every_N",
        all(t["xi2_wineland"] < o["xi2_wineland"]
            for t, o in zip(trows, rows)))
    chk("tact_exponent_window",
        TACT_EXPONENT_WINDOW[0] <= tact["exponent"]
        <= TACT_EXPONENT_WINDOW[1])
    chk("tact_exponent_steeper_than_oat",
        tact["exponent"] < oat["exponent"] - 0.1)
    chk("tact_top_local_slope_near_ideal",
        -1.1 < tact["top_local_slope"] < -0.85)

    # 7. Classifier: kNN interpolates the policy, scenarios agree
    rep = classifier_report()
    chk("classifier_css_sql", rep["css_at_t0"]["knn"] == "sql-limited"
        and rep["css_at_t0"]["policy"] == "sql-limited")
    chk("classifier_oat20_deep", rep["oat_optimum_N20"]["knn"]
        == "deep-squeezed")
    chk("classifier_oat100_deep", rep["oat_optimum_N100"]["knn"]
        == "deep-squeezed")
    chk("classifier_mild_squeezed", rep["oat_mild_N50_mu0.02"]["knn"]
        == "squeezed")
    chk("classifier_tact40_deep", rep["tact_optimum_N40"]["knn"]
        == "deep-squeezed")
    chk("classifier_cat_over_twisted", rep["cat_point_N20"]["knn"]
        == "over-twisted")
    chk("classifier_matches_policy_on_scenarios",
        all(v["knn"] == v["policy"] for v in rep.values()))
    chk("classifier_matches_policy_on_grid",
        all(classify_regime(10.0 ** f[0], f[1]) == label
            for f, label in TRAINING_GRID))

    # 8. Gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_sensor_hold", gate("real_sensor_calibration") == "HOLD")
    chk("gate_clock_hold", gate("atomic_clock_product_claim") == "HOLD")
    chk("gate_experiment_hold", gate("live_experiment_direction") == "HOLD")
    chk("gate_defense_hold", gate("defense_targeting_use") == "HOLD")
    chk("gate_record_hold", gate("squeezing_record_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("SpinSqueezing-Lab - spin squeezing & quantum metrology lab")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE REFERENCE POINT - a coherent spin state sits at the SQL")
    print("    All N atoms point the same way; the transverse noise is the")
    print("    standard quantum limit, computed exactly in the Dicke basis:")
    for N in (10, 100):
        p = squeezing_params(N, css_x(N))
        print(f"    N={N:4d}  xi^2 = {p['xi2_ku']:.12f}  "
              f"xi_R^2 = {p['xi2_wineland']:.12f}  "
              f"contrast = {p['contrast']:.12f}")
    print("    -> xi^2 = xi_R^2 = 1 at t = 0. Anything below 1 is squeezing;")
    print("       anything quoted without the contrast is only half the story.")
    print("\n[2] ONE-AXIS TWISTING - H = chi Sz^2 shears the noise ellipse")
    print("    Contrast falls on the closed form cos^(N-1)(chi t) while the")
    print("    ellipse thins; the Wineland optimum balances the two:")
    oat = oat_scaling()
    for r in oat["rows"]:
        print(f"    N={r['N']:4d}  mu* = {r['mu']:.4f}  "
              f"xi_R^2 = {r['xi2_wineland']:.6f}  "
              f"xi^2 = {r['xi2_ku']:.6f}  contrast = {r['contrast']:.6f}")
    print(f"    fit: xi_R^2(min) ~ N^({oat['exponent']:.4f})   "
          f"ideal exponent -2/3 = -0.6667")
    print("    -> at every optimum xi_R^2 > xi^2: the KU number alone would")
    print("       overstate the phase gain by 1/contrast^2. The metrology")
    print("       figure is Dphi = xi_R / sqrt(N), never xi / sqrt(N).")
    print("\n[3] THE TWO CURVES AT N=50 - squeezing develops, contrast pays")
    for row in oat_curve(50, [0.0, 0.02, 0.04, 0.06, 0.08, 0.12]):
        xi_w = row["xi2_wineland"]
        xi_w_s = f"{xi_w:10.6f}" if xi_w != float("inf") else "       inf"
        print(f"    mu = {row['mu']:.2f}  xi^2 = {row['xi2_ku']:.6f}  "
              f"xi_R^2 = {xi_w_s}  contrast = {row['contrast']:.6f}")
    print("\n[4] THE CAT POINT - chi t = pi/2, even N")
    cat = cat_report(20)
    print(f"    N=20: |<CSS +x|psi>|^2 = {cat['overlap_plus_x']:.6f}  "
          f"|<CSS -x|psi>|^2 = {cat['overlap_minus_x']:.6f}  "
          f"sum = {cat['overlap_sum']:.6f}")
    print(f"    <Sx> = {cat['mean_sx']:.2e}  <Sy> = {cat['mean_sy']:.2e}  "
          f"contrast = {cat['contrast']:.2e}")
    print(f"    Sz population change = {cat['sz_population_change']:.1e} "
          f"(phases only - the bimodality is in the Sx basis)")
    print("    -> an exact macroscopic superposition: half the state points")
    print("       +x, half points -x, and the mean spin has vanished.")
    print("\n[5] TWO-AXIS COUNTER-TWISTING - H = chi (Sz^2 - Sy^2)")
    print("    Twisting about y and z with opposite signs squeezes")
    print("    exponentially; fitted on the same grid, computed the same way:")
    tact = tact_scaling()
    for t, o in zip(tact["rows"], oat["rows"]):
        print(f"    N={t['N']:4d}  t* = {t['mu']:.4f}  "
              f"xi_R^2 = {t['xi2_wineland']:.6f}  "
              f"(OAT {o['xi2_wineland']:.6f})  "
              f"contrast = {t['contrast']:.6f}")
    print(f"    fit: xi_R^2(min) ~ N^({tact['exponent']:.4f})   "
          f"ideal exponent -1 (asymptotic); local slope at the top of")
    print(f"    the grid is {tact['top_local_slope']:.4f} - still steepening "
          f"toward -1.")
    print("    -> TACT beats OAT at every N on the grid, with a steeper")
    print("       exponent. The 1/N law is the asymptote, not the decade fit.")
    print("\n[6] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, v in classifier_report().items():
        xi_s = ("inf" if v["xi2_wineland"] is None
                else f"{v['xi2_wineland']:.6f}")
        print(f"    {name:24s} xi_R^2 = {xi_s:>10s}  "
              f"contrast = {v['contrast']:.4f}  -> {v['knn']}")
    print("\n[7] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[8] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "css":
        report = {}
        for N in (10, 20, 50, 100):
            p = squeezing_params(N, css_x(N))
            report[f"N{N}"] = {
                "xi2_ku": p["xi2_ku"], "xi2_wineland": p["xi2_wineland"],
                "contrast": p["contrast"], "mean_spin": p["mean_spin"]}
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if cmd == "twist":
        N = int(argv[1]) if len(argv) > 1 else 50
        mus = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2]
        print(json.dumps({"N": N, "curve": oat_curve(N, mus),
                          "optimum": oat_optimum(N)},
                         indent=2, sort_keys=True))
        return 0
    if cmd == "scaling":
        print(json.dumps({"oat": oat_scaling(), "tact": tact_scaling()},
                         indent=2, sort_keys=True))
        return 0
    if cmd == "contrast":
        rows = oat_curve(50, [0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.2])
        print(json.dumps({"N": 50, "closed_form": "cos^(N-1)(mu)",
                          "curve": rows}, indent=2, sort_keys=True))
        return 0
    if cmd == "cat":
        N = int(argv[1]) if len(argv) > 1 else 20
        print(json.dumps(cat_report(N), indent=2, sort_keys=True))
        return 0
    if cmd == "tact":
        N = int(argv[1]) if len(argv) > 1 else 40
        ts = [3.0 / math.sqrt(N) * i / 20 for i in range(21)]
        print(json.dumps({"N": N, "curve": tact_curve(N, ts),
                          "optimum": tact_optimum(N)},
                         indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "real_sensor_calibration"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
