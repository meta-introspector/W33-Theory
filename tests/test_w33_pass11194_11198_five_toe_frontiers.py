"""Focused regression for the five exact Pass 11194--11198 frontier audits."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass11194_11198_five_toe_frontiers as P  # noqa: E402
import w33_qutrit_clifford_phase_displacement_lift as PHASE  # noqa: E402

D = json.loads((ROOT / "data/w33_pass11194_11198_five_toe_frontiers.json").read_text())


def test_frozen_packet():
    assert D["checks"]["all_five_pass"]
    assert set(D["passes"]) == {"11194", "11195", "11196", "11197", "11198"}
    assert all(all(row["checks"].values()) for row in D["passes"].values())


def test_reversible_measure_and_hecke_traces_live():
    y = P.pass11194_reversible_yukawa_measure()
    assert y["solution"] == "a=b" and y["cross_edges"] == 40
    assert y["stationary_sector_masses"] == ["1/9", "8/9"]
    assert P.hecke_trace(1) == 0
    assert P.hecke_trace(2) == Fraction(15, 4)
    assert all(P.hecke_trace(n) for n in range(2, 20))


def test_signed_incidence_firewall_live():
    s = P.pass11196_signed_incidence_firewall()
    assert s["unsigned_rank"] == s["signed_rank"] == 21
    assert s["centered_gram_spectrum"] == {"6": 20, "0": 7}
    assert s["canonical_sign_profile"] == {"plus": 22, "minus": 23}


def test_heat_limit_live():
    rows = [P.w3q_heat_row(q) for q in (3, 9, 27, 81, 243)]
    errors = [abs(r["error_from_exp_minus_1"]) for r in rows]
    assert all(r["diameter"] == 2 for r in rows)
    assert all(errors[i + 1] < errors[i] for i in range(len(errors) - 1))


def test_phase_join_frozen_words_live():
    c = D["passes"]["11195"]
    words = {
        name: tuple(tuple(x) for x in row["word"])
        for name, row in c["representative_words"].items()
    }
    normal = words["L"] + words["p"] + words["h_star"] + words["p"] + words["R"]
    assert PHASE.word_matrix(normal) == PHASE.word_matrix(words["SUM"])
    scalar, residual = P._phase_scalar(PHASE.word_unitary(words["SUM"]), PHASE.word_unitary(normal))
    assert np.isclose(scalar, 1) and residual < 1e-12

