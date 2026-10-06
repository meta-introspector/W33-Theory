"""Regression for Passes 11511-11515 (the good rules for every n via the unipotent part; the cubed-phase test exhaustively
up to symmetry; exact P10; six qutrits to 4000 classes; the five MUB-chirality generators and the degree-13 syzygy)."""

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11511_exhaustion_and_coverage():
    d = load("w33_pass11511_good_rules_all_n.json")
    ex = d["exhaustion_U"]
    assert ex["Sp(2,3)"]["U_holds"] == ex["Sp(2,3)"]["unipotent"] == 9
    assert ex["Sp(4,3)"]["U_holds"] == ex["Sp(4,3)"]["unipotent"] == 6561
    c = d["n2"]["counts"]
    assert not any("NOT covered" in k for k in c)
    for cell in ("M z1 = -z1", "M^2 z1 = -z1"):
        assert c[f"{cell}: decided"] == c[f"{cell}: decided all frames reversible"]
    for k, v in d["exhaustion_E"].items():
        assert v["failures"] == 0 and v["vectors"] > 0, k
    for n in ("n2_blocks", "n3_blocks"):
        assert not any("outside" in k for k in d[n]["counts"]), n


def test_11511_lagrangian_lemma_fast():
    """the semisimple case: for M1 = I on F3^4 and any a, an anti-symplectic involution with -1-space a Lagrangian
    through a satisfies Im(A^-1 - I) perp a"""
    import w33_pass11511_good_rules_all_n as G
    import w33_pass11350_linear_decider as L
    import w33_pass11330_orbit_census as O
    D2 = L.Decider(2)
    S4 = np.array(O.all_symplectic(D2.wl)[0]) % 3
    A4 = (S4 @ np.diag([1, 2, 1, 2])) % 3
    ok, nrev = G.U_holds(np.eye(4, dtype=np.int64), A4[:5000], D2.wl.Om)
    assert ok and nrev > 0


def test_11512_cubed_exhaustive():
    d = load("w33_pass11512_cubed_phase_exhaustive.json")
    c = d["counts"]
    assert c["orbit representatives (L, non-affine f mod affine)"] == 816
    sc = d["strong"]["counts"]
    assert sc["unitary-containing eigenspaces"] == 1212 and sc["... NOT reversible"] == 2
    ce = d["counterexample"][0]
    assert ce["unitarity_error"] < 1e-12 and ce["prefilter_passes"] and ce["max_overlap"] < 3 - 1e-3
    assert ce["pass11252_decider"] is False
    U = np.array(ce["U_real"]) + 1j * np.array(ce["U_imag"])
    assert np.allclose(U @ U.conj().T, np.eye(3), atol=1e-12)


def test_11513_exact_depth10():
    d = load("w33_pass11513_exact_depth10.json")
    lv = d["levels"]
    assert lv["8"]["matches_known"] and lv["9"]["matches_known"]
    assert lv["8"]["distinct_operators"] == 32445792
    assert d["dedupe_checks"]["8"]["classes_at_1e6"] == d["dedupe_checks"]["8"]["classes_at_1e5"]
    assert all(d["denominator_law"])
    P10 = Fraction(lv["10"]["P"])
    assert 0 < P10 < Fraction(11519231, 201326592)


def test_11514_six_qutrits():
    d = load("w33_pass11514_six_qutrits_4000.json")
    assert d["classes"] >= 3000 and d["estimate"]["stderr"] < 0.0065


def test_11515_generators_and_syzygy():
    d = load("w33_pass11515_odd_generators_and_syzygy.json")
    for k in range(6, 11):
        g = d["generators"][str(k)]
        assert g["exponents"] == [k - 3, 2, 1] and len(g["mubs"]) == 3
    s = d["syzygies"]
    assert all(s[str(k)]["relations"] == 0 for k in range(6, 13)) and s["13"]["relations"] == 1
    assert d["predicted_relations_from_molien"]["13"] == 1
    kept, dropped = s["13"]["gap"]
    assert kept > 1e-8 and dropped < 1e-12
