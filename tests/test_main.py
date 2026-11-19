"""Tests for SpinSqueezing-Lab: exact Dicke-basis spin squeezing."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import main


# ---------------------------------------------------------------- CSS / SQL
def test_css_sql_exact():
    for N in (10, 50):
        p = main.squeezing_params(N, main.css_x(N))
        assert abs(p["xi2_ku"] - 1.0) < 1e-9
        assert abs(p["xi2_wineland"] - 1.0) < 1e-9
        assert abs(p["contrast"] - 1.0) < 1e-12
        assert abs(p["expectations"]["sx"] - N / 2.0) < 1e-9


def test_css_any_direction_is_sql():
    p = main.squeezing_params(30, main.css(30, 0.7, 1.1))
    assert abs(p["xi2_ku"] - 1.0) < 1e-9
    assert abs(p["contrast"] - 1.0) < 1e-9


def test_dicke_basis_dimensions():
    assert len(main.dicke_ms(10)) == 11
    assert main.dicke_ms(10)[0] == -5.0
    assert main.dicke_ms(10)[-1] == 5.0
    assert len(main.css_x(7)) == 8  # odd N: half-integer m, still fine


# ---------------------------------------------------------------- OAT
def test_contrast_closed_form():
    N = 40
    psi0 = main.css_x(N)
    for mu in (0.05, 0.1, 0.3):
        p = main.squeezing_params(N, main.oat_evolve(N, psi0, mu))
        assert abs(p["contrast"] - math.cos(mu) ** (N - 1)) < 1e-9


def test_oat_sz_populations_frozen():
    psi0 = main.css_x(40)
    psi_t = main.oat_evolve(40, psi0, 0.37)
    assert max(abs(abs(a) ** 2 - abs(b) ** 2)
               for a, b in zip(psi0, psi_t)) < 1e-15
    assert abs(sum(abs(a) ** 2 for a in psi_t) - 1.0) < 1e-12


def test_oat_optimum_below_sql():
    opt = main.oat_optimum(20)
    assert opt["xi2_wineland"] < 1.0
    assert opt["xi2_ku"] < opt["xi2_wineland"]  # contrast < 1 at optimum
    assert opt["contrast"] < 1.0


def test_oat_scaling_exponent_window():
    oat = main.oat_scaling()
    lo, hi = main.OAT_EXPONENT_WINDOW
    assert lo <= oat["exponent"] <= hi
    assert abs(oat["exponent"] - (-2.0 / 3.0)) < 0.1


def test_oat_optima_improve_with_N():
    rows = main.oat_scaling()["rows"]
    assert all(rows[i + 1]["xi2_wineland"] < rows[i]["xi2_wineland"]
               for i in range(len(rows) - 1))


def test_wineland_is_ku_over_contrast_squared():
    for r in main.oat_scaling()["rows"]:
        assert abs(r["xi2_wineland"]
                   - r["xi2_ku"] / r["contrast"] ** 2) < 1e-6


# ---------------------------------------------------------------- cat point
def test_cat_point_exact_superposition():
    cat = main.cat_report(20)
    assert abs(cat["overlap_plus_x"] - 0.5) < 1e-9
    assert abs(cat["overlap_minus_x"] - 0.5) < 1e-9
    assert abs(cat["overlap_sum"] - 1.0) < 1e-9
    assert abs(cat["mean_sx"]) < 1e-9
    assert cat["contrast"] < 1e-9
    assert cat["sz_population_change"] < 1e-15


def test_cat_point_N40():
    assert abs(main.cat_report(40)["overlap_sum"] - 1.0) < 1e-9


def test_cat_requires_even_N():
    try:
        main.cat_report(15)
    except ValueError:
        return
    raise AssertionError("cat_report must reject odd N")


# ---------------------------------------------------------------- TACT
def test_tact_matrix_structure():
    H = main.tact_matrix(12)
    d = len(H)
    assert all(H[i][j] == H[j][i] for i in range(d) for j in range(d))
    assert all(H[i][j] == 0.0 for i in range(d) for j in range(d)
               if abs(i - j) not in (0, 2))


def test_jacobi_referee():
    H = main.tact_matrix(12)
    evals, V = main.jacobi_eigh(H)
    d = len(H)
    recon = max(abs(sum(V[i][k] * evals[k] * V[j][k] for k in range(d))
                    - H[i][j]) for i in range(d) for j in range(d))
    ortho = max(abs(sum(V[k][i] * V[k][j] for k in range(d))
                    - (1.0 if i == j else 0.0))
                for i in range(d) for j in range(d))
    assert recon < 1e-9
    assert ortho < 1e-9


def test_tact_norm_preserved():
    psi_t = main.tact_evolve(20, main.css_x(20), 0.05)
    assert abs(sum(abs(a) ** 2 for a in psi_t) - 1.0) < 1e-9


def test_tact_beats_oat_every_N():
    oat_rows = main.oat_scaling()["rows"]
    tact_rows = main.tact_scaling()["rows"]
    assert all(t["xi2_wineland"] < o["xi2_wineland"]
               for t, o in zip(tact_rows, oat_rows))


def test_tact_exponent_window_and_steeper():
    tact = main.tact_scaling()
    oat = main.oat_scaling()
    lo, hi = main.TACT_EXPONENT_WINDOW
    assert lo <= tact["exponent"] <= hi
    assert tact["exponent"] < oat["exponent"] - 0.1
    assert -1.1 < tact["top_local_slope"] < -0.85


# ---------------------------------------------------------------- classifier
def test_classifier_scenarios():
    rep = main.classifier_report()
    assert rep["css_at_t0"]["knn"] == "sql-limited"
    assert rep["oat_optimum_N20"]["knn"] == "deep-squeezed"
    assert rep["oat_optimum_N100"]["knn"] == "deep-squeezed"
    assert rep["oat_mild_N50_mu0.02"]["knn"] == "squeezed"
    assert rep["tact_optimum_N40"]["knn"] == "deep-squeezed"
    assert rep["cat_point_N20"]["knn"] == "over-twisted"
    assert all(v["knn"] == v["policy"] for v in rep.values())


def test_classifier_matches_policy_on_grid():
    for features, label in main.TRAINING_GRID:
        assert main.classify_regime(10.0 ** features[0],
                                    features[1]) == label


# ---------------------------------------------------------------- gates/seats
def test_hold_gates():
    assert len(main.HOLD_GATES) == 14
    assert all(v.startswith("HOLD") for v in main.HOLD_GATES.values())
    assert main.gate("real_sensor_calibration") == "HOLD"
    assert main.gate("atomic_clock_product_claim") == "HOLD"
    assert main.gate("live_experiment_direction") == "HOLD"
    assert main.gate("defense_targeting_use") == "HOLD"
    assert main.gate("squeezing_record_claim") == "HOLD"
    assert main.gate("run_simulation") == "PROCEED"


def test_hold_demo_exit_codes():
    assert main.main(["hold-demo"]) == 2
    assert main.main(["hold-demo", "run_simulation"]) == 0


def test_seats_never_override():
    brains = main.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = main.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_demo_deterministic():
    root = Path(__file__).resolve().parent.parent
    out1 = subprocess.run([sys.executable, "main.py", "demo"],
                          cwd=root, capture_output=True).stdout
    out2 = subprocess.run([sys.executable, "main.py", "demo"],
                          cwd=root, capture_output=True).stdout
    assert out1 == out2
    assert b"CHECKS 49/49" in out1
