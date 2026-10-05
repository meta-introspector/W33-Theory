"""Regression for Passes 11456-11460 (four-qutrit closure; the magic-axis theorem; the generic blind spot; deep decay;
five qutrits)."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def load(name):
    return json.load(open(ROOT / "data" / name))


def test_11456_closure():
    d = load("w33_pass11456_n4_closure.json")
    f = d["fixed_axis_cell"]
    assert sum(f["classes"].values()) == 549 and f["group_order"] == 20056328248320
    assert f["classes"]["decided by the decider"] == 232
    assert "FAILS" not in d["same_line_cell_n4"]["verdicts"]


def test_11457_theorem():
    d = load("w33_pass11457_magic_axis_theorem.json")
    assert len(set(round(x, 6) for x in d["tau_moduli"])) == 3
    assert all(s["classes"] == s["support_is_Im_M_minus_I_uniform"] for s in d["support_lemma"])
    n2 = d["n2"]
    assert n2["sound (all theorem frames violating)"] == n2["classes"] == 376 and n2["Rad dim 0"] == 352
    assert d["n4_unreached_classes"]["count"] == 317
    tm = d["extension_route_twisted_moduli"]
    assert np.allclose(tm["0"], 1) and len(set(round(x, 6) for x in tm["1"])) == 3
    m = d["magnitude_obstruction_n2"]
    assert m["Rad = 0 | all F' frames certified"] == 16 and not any(k.startswith("Rad = 0 |") and "all" not in k for k in m)


def test_11457_theorem_frames_fast():
    import w33_pass11350_linear_decider as L
    import w33_pass11457_magic_axis_theorem as TH
    D = L.Decider(2)
    M = np.eye(4, dtype=np.int64)
    M[0, 1] = 1                                    # a transvection fixing z1? check the hypotheses handle it gracefully
    out = TH.theorem_frames(M % 3, D.z1, D.wl.Om, D.wl.labels.astype(np.int64))
    assert out is None or out[0].dtype == bool


def test_11458_generic():
    d = load("w33_pass11458_generic_component.json")
    assert all(r["dF_rank"] == 4 and r["residual_norms"][-1] < 1e-12 and r["J8"] > 1e-3 for r in d["points"])
    assert min(d["continuous_symmetry_singular_values"]) > 0.1
    mins = d["min_eigvec_Delta6_after"]
    assert min(mins) < 1e-12 < max(mins)           # achiral eigenvector at some points, not at others


def test_11459_deep():
    d = load("w33_pass11459_deep_decay.json")
    assert all(v["samples"] > 9_000_000 for v in d["depths"].values())
    assert all(0.70 < r["rate"] < 0.74 for r in d["local_rates"])                     # constant ~0.72 (tight cut)
    rho88 = max(np.roots([2592, -1332, -585, 140]).real)
    assert all(abs(r["rate"] - rho88) < 3 * r["stderr"] for r in d["local_rates"][1:])
    s30 = d["threshold_sensitivity"]["30"]
    assert s30["false_fraction_at_2999"] > 0.5 and s30["counts"]["3-1e-6"] == s30["counts"]["3-1e-10"]
    assert all(0.68 < r["rate"] < 0.72 for r in d["no_identity_coset_walk"]["local_rates"])


def test_11460_five_qutrits():
    d = load("w33_pass11460_five_qutrits.json")
    e = d["n5"]["estimate"]
    assert e["shares"] == {"collinear": "2460/7381", "noncollinear": "19683/29524"}
    assert e["P_estimate"] > 0.125 + 2 * e["stderr"]
    assert "FAILS" not in d.get("n5_fixed_axis_F_prime", {})
