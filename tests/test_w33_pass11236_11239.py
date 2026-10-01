"""Regression for Passes 11236-11239 (frozen certificates plus fast recomputations)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def load(name):
    return json.loads((ROOT / "data" / name).read_text())


def test_11236_no_complete_pattern():
    d = load("w33_pass11236_cubic_flavour.json")
    assert d["n_pmns_complete_hits"] == 0 and d["n_cabibbo_hits"] == 0 and d["n_distinct_modulus_patterns"] == 6
    for k in ("magic_e_x_clifford_nu", "clifford_e_x_magic_nu", "magic_e_x_magic_nu"):
        assert d[k]["pmns_hits"] == 0 and d[k]["cabibbo_hits"] == 0
    assert d["depth1_finite_order_bases"] == 148


def test_11236_T_eigenbasis_is_stabilizer():
    import w33_pass11213_cubic_t_violation as P
    import w33_pass11236_cubic_flavour as F
    # T is diagonal: its eigenbasis is the Z eigenbasis, so <T> and <Z> give the same U_e
    assert F.basis_key(F.eigbasis(P.T1)) == F.basis_key(F.eigbasis(P.Z1))


def test_11237_levels():
    d = load("w33_pass11237_fidelity_levels.json")
    lv = {round(x["level"], 6): x for x in d["levels"]}
    assert lv[0.84403]["F2_minpoly"] == [81, -81, 18, -1]
    assert lv[0.810095]["F2_minpoly"] == [-59049, 39366, -405, 1]
    assert lv[0.84403]["F2_in_Q_cos_2pi_9"] and lv[0.810095]["F2_in_Q_cos_2pi_9"]
    assert d["F_min_squared_is_the_two_qutrit_level_0_712386014"]
    x = QUANTUM ** 2
    assert abs(81 * x ** 3 - 81 * x ** 2 + 18 * x - 1) < 1e-12


def test_11238_sectors():
    d = load("w33_pass11238_one_sector_law.json")
    assert d["control_sector"]["verdicts"] == {"reversible": 13}
    assert d["target_sector"]["verdicts"]["violating"] == 12
    assert set(d["target_sector"]["violator_fidelities"]) == {"0.712386014"}
    assert d["both_sectors"]["verdicts"]["violating"] == 20
    assert abs(QUANTUM ** 2 - 0.712386014) < 1e-9


def test_11239_three_qutrit_bounds():
    d = load("w33_pass11239_three_qutrit_t.json")
    c = d["candidates"]
    assert abs(c["(T(x)T)SUM (x) I"]["lower_bound"] - QUANTUM) < 1e-9
    assert abs(d["factorised"]["(T(x)T(x)T) SUM12"]["product_bound"] - QUANTUM) < 1e-9
    assert c["(I(x)I(x)T) SUM23 (x) control"]["certified_reversible"]
    # only the factorised one-cubic-gate control is certified; the rest are lower bounds below 1
    assert [k for k, v in c.items() if v["certified_reversible"]] == ["(I(x)I(x)T) SUM23 (x) control"]
    # the generic search is weak: it misses the exact factorised bound
    assert c["(T(x)T(x)T) SUM12"]["random"] < 0.5 < QUANTUM


def test_11239_factorisations():
    import w33_pass11239_three_qutrit_t as M
    U2 = np.kron(M.T, M.T) @ M.P2.SUM
    assert np.allclose(M.kron3(M.T, M.T, M.T) @ M.sum_gate(3, 0, 1), np.kron(U2, M.T))
    W2 = np.kron(M.I3, M.T) @ M.P2.SUM
    assert np.allclose(M.kron3(M.I3, M.I3, M.T) @ M.sum_gate(3, 1, 2), np.kron(M.I3, W2))
    # a diagonal tick is reversed by complex conjugation alone: tr(T^* T)/3 = 1
    assert abs(np.trace(M.T.conj() @ M.T) / 3 - 1) < 1e-12
