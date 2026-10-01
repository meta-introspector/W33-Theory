"""Regression for Passes 11250-11254 (frozen certificates plus fast exact recomputations)."""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11250_sp44_identities():
    d = load("w33_pass11250_qubit_identities_proved.json")
    assert d["sp44_order_enumerated"] == d["sp44_order_formula"] == 979200
    assert d["unipotents"] == d["steinberg_q8"] == 65536
    assert d["identity1_all"] and d["identity2_all"]
    assert d["identity2_q4"]["2,2"]["fraction"] == "1/16"


def test_11250_form_counts():
    import w33_pass11250_qubit_identities_proved as M
    for k in range(1, 9):
        assert M.p_chi0(k, 2) == (Fraction(1, 2 ** k) if k % 2 == 0 else 0)        # Pass 11245's lemma
    for k in (2, 4, 6):
        assert M.p_chi0(k, 4) == Fraction(1, 4 ** k)
    assert M.nondeg_symmetric(2, 4) == 48 and M.nondeg_symmetric(2, 2) == 4


def test_11251_exact_identity():
    import w33_pass11251_exact_reversal as E
    d = load("w33_pass11251_exact_reversal.json")
    assert d["proved_F_T_equals_F_min_squared"] and d["F_T_squared_equals_F_min_fourth_exactly"]
    assert d["n_maximisers"] == 243 and d["maximiser_r"] == 0
    t = tuple(d["exact_trace_coeffs"])
    g = E.add(E.add(E.ONE, E.zpow(1)), E.zpow(8))
    assert E.mul(t, E.conj(t)) == E.mul(E.mul(g, g), E.mul(g, g))                  # t tbar = (1 + z + z^8)^4
    # the mechanism: phases {0,0,0,+-1,+-1,+-2} = sumset {-1,0,1}+{-1,0,1}
    s = E.ZERO
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            s = E.add(s, E.zpow(i + j))
    assert s == E.mul(g, g)
    assert d["galois_gap_in_F"] > 10 * d["float_error_bound_used"] / 2


def test_11253_cubic_field():
    import w33_pass11253_cubic_field_levels as L
    d = load("w33_pass11253_cubic_field_levels.json")
    assert d["distinct_levels_one_qutrit"] == 103 and d["distinct_levels_two_qutrit"] == 45
    assert d["cubic_levels"] == 147 and d["rational_levels"] == 1
    assert d["levels_with_F_itself_in_Q_cos_2pi_9"] == 143
    assert d["max_float_vs_exact_F2_deviation"] < 1e-13
    c = 2 * np.cos(2 * np.pi / 9)
    for lv in d["levels"]:
        a, b, e = (float(Fraction(x)) for x in lv["F2_abc"])
        assert abs(a + b * c + e * c * c - lv["F"] ** 2) < 1e-12
    assert L.sqrt_in_field([Fraction(25, 81), Fraction(10, 81), Fraction(1, 81)]) == ["5/9", "1/9", "0"]
    assert L.sqrt_in_field([Fraction(0), Fraction(16, 81), Fraction(20, 81)]) == ["-4/9", "2/9", "4/9"]


def test_11253_ring_of_integers():
    # disc(x^3 - 3x + 1) = -4(-3)^3 - 27 = 81 = disc of Q(cos 2pi/9): Z[c] is the full ring of integers
    assert -4 * (-3) ** 3 - 27 * 1 ** 2 == 81


def test_11254_tm1_breaking():
    d = load("w33_pass11254_tm1_yukawa_breaking.json")
    assert d["equivariant"] and d["models"] == 60
    assert abs(d["eps0_s12_tm1"] - 0.3183716) < 1e-6
    assert d["linear_regime_nu_only"] == 1.0 and d["linear_regime_e_only"] == 1.0
    assert 5 < d["abs_dUe1_per_radian_tilt_nu_only"]["median"] < 30
    assert 5 < d["abs_dUe1_per_radian_tilt_e_only"]["median"] < 40
    assert d["eps0_s23_range"][0] < 0.35 and d["eps0_s23_range"][1] > 0.65      # theta23 free at eps = 0


def test_11254_unique_cubic_equivariant():
    import w33_pass11254_tm1_yukawa_breaking as Y
    rng = np.random.default_rng(0)
    p = rng.normal(size=3)
    for g in Y.G:
        assert np.allclose(Y.kcubic(g @ p), g @ Y.kcubic(p) @ g.T)
    chi0 = np.array([0, 1, 1]) / np.sqrt(2)
    R = np.diag([-1, 1, 1])
    quad = Y.mnu_terms(chi0)[:4]
    assert all(np.allclose(R @ M @ R, M) for M in quad)                            # the accidental reflection
    assert not np.allclose(R @ Y.kcubic(chi0) @ R, Y.kcubic(chi0))                  # broken by K


def test_11252_exact_decisions():
    d = load("w33_pass11252_exact_reversibility.json")
    assert d["two_qutrit_validation"]["disagree"] == 0 and d["two_qutrit_validation"]["agree"] == 204
    c = d["three_qutrit_candidates"]
    assert c["(T(x)T(x)T) SUM23 SUM12"]["reversible"] is False
    assert c["(T(x)I(x)T) SUM23 SUM12"]["reversible"] is True
    assert c["(T(x)T)SUM (x) I"]["reversible"] is False
    assert d["three_qutrit_conjugated_control_reversible"] == "6/6"
    assert all(v["undecided"] == 0 for v in d["three_qutrit_census"].values())


def test_11252_criterion_fast():
    import w33_pass11227_two_qutrit_cp_mixing as P2
    import w33_pass11252_exact_reversibility as R
    R.WEYL[2] = R.Weyl(2)
    rng = np.random.default_rng(3)
    T, I3 = P2.T, P2.I3
    assert R.decide(np.kron(T, T) @ P2.SUM, 2, rng)[0] is False
    assert R.decide(np.kron(I3, T) @ P2.SUM, 2, rng, twirl=True)[0] is True
    assert R.decide(R.random_clifford(2, rng), 2, rng, twirl=True)[0] is True           # no Clifford T-violation
