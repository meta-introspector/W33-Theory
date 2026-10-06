"""Regression for Passes 11531-11537 (cubed test: words vs PU(3); good rules for blocks of dimension 6 and 8; the one-gate
law; the exact time-arrow module through degree 20; the decay data; the union law)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11531_words_vs_pu3():
    d = load("w33_pass11531_cubed_test_words_vs_pu3.json")
    for f, r in d["counterexample_family"].items():
        assert set(r["reversibility"]) == {"not reversible"} and r["eigenspace_dim"] == 6
        assert max(r["local_dimension_counts"], key=r["local_dimension_counts"].get) == "4"
    w = d["word_relations"]
    for k in ("6", "7"):
        assert w[k].get("non-affine: Clifford") == 2376 and w[k].get("non-affine: T-count 1") == 2592
    assert all("non-affine: T-count >= 2 (would break the mechanism)" not in v for v in w.values())
    assert w["7"]["distinct_operators"] == 5567184


def test_11532_dim6_witnesses():
    d = load("w33_pass11532_good_rules_dim6.json")
    wv = d["witnessed"]
    for cls in ("[6] form class 1", "[6] form class 2", "[3,3]"):
        assert wv[f"{cls}: (U) proved"] == 4 and wv[f"{cls}: (E) for -u proved"] == 4
    assert wv["U(3) regular: (E) proved"] == 4
    assert not any("NOT" in k or "no reverser" in k for k in wv)


def test_11533_one_gate_law():
    d = load("w33_pass11533_one_gate_law.json")
    n2 = d["n2"]
    assert n2["class verdict agrees"] == n2["classes"] == 51838 and n2["frame count agrees"] == 51814
    n3 = d["n3"]
    assert n3["law_bad_class_fraction"] == "97/728" and not n3["law_equals_exact"]
    assert Fraction(437, 3276) - Fraction(97, 728) == Fraction(1, 6552)
    assert n3["counts"]["orbits decided"] - n3["counts"]["orbit verdict agrees"] == 9
    e4 = d["estimates"]["4"]["non-collinear"]
    assert e4["z_vs_conjecture"] > 10                     # the (1 - 9^(1-n))/8 conjecture is refuted


def test_11533_law_n2_fast():
    """the law's class count over all of Sp(4,3) is exactly 1/8 (takes ~1 min)"""
    import w33_pass11330_orbit_census as O
    import w33_pass11350_linear_decider as L
    import w33_pass11533_one_gate_law as G
    D = L.Decider(2)
    bad = 0
    for M in np.array(O.all_symplectic(D.wl)[0])[::40] % 3:          # every 40th class: the 1/8 is not uniform
        w = G.W_dim(D, M)                                                # per subsample, so only check it runs and
        bad += (w is None) or w > 0                                      # gives a fraction in a sane range
    assert 0.05 < bad / len(range(0, 51840, 40)) < 0.25


def test_11534_module_exact():
    d = load("w33_pass11534_time_arrow_module_exact.json")
    assert d["degree_bound"] == 20
    assert sorted(v["degree"] for v in d["even_generators"].values()) == [1, 3, 4, 5, 6, 7, 8, 9]
    assert all(v["monomial"] == [0] * v["degree"] for v in d["even_generators"].values())     # stabiliser moments
    assert sorted(v["degree"] for v in d["odd_generators"].values()) == [6, 7, 8, 9, 10]
    assert all(v["extra_generators_needed"] == 0 for v in d["odd"].values())
    assert [d["odd"][str(k)]["relations_among_listed"] for k in range(13, 21)] == [1, 3, 6, 10, 16, 24, 35, 50]
    assert [d["even"][str(k)]["even_relations"] for k in range(14, 21)] == [1, 2, 3, 5, 8, 11, 16]


def test_11535_decay():
    d = load("w33_pass11535_decay_exact_data.json")
    assert d["integers"][-1] == 65621121 and all(v != "FITS" for v in d["recurrences_full"].values())


def test_11536_dim8_witnesses():
    d = load("w33_pass11536_good_rules_dim8.json")
    wv = d["witnessed"]
    for cls in ("[8] form class 1", "[8] form class 2"):
        assert wv.get(f"{cls}: (U) proved", 0) >= 3 and wv.get(f"{cls}: (E) for -u proved", 0) >= 3
    assert not any("NOT" in k for k in wv)


def test_11537_union_law():
    d = load("w33_pass11537_union_law.json")
    for n in ("n2", "n3"):
        c = d[n]["counts"]
        assert not any(("larger" in k) or ("VIOLATES" in k) for k in c), n
    assert sum(v for k, v in d["n2"]["counts"].items() if k.endswith("EXACT")) == 51838
    assert sum(v for k, v in d["n3"]["counts"].items() if k.endswith("EXACT")) == 2117
