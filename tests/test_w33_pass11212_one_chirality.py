"""Regression for Pass 11212: one chirality (chiral representatives recomputed; member sweeps frozen)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def _data():
    return json.loads((ROOT / "data" / "w33_pass11212_one_chirality.json").read_text())


def test_representatives_recomputed():
    import w33_pass11212_one_chirality as C
    d = _data()
    reps = {int(k): v for k, v in d["representatives"].items()}
    pair256 = sorted(k for k, v in reps.items() if v["size"] == 256)
    for o in pair256:
        B = np.array(reps[o]["B"], np.int64)
        chi = C.M.chirality(C.M.spectrum(B))
        q = C.H.twist(B)
        q = 1 if q == 1 else -1
        im, _ = C.value(B, C.G.CHIRAL, C.UNITS["CHIRAL"])
        assert chi == 54 * q and im == q
    chiral6912 = sorted(k for k, v in reps.items() if v["size"] == 6912 and v["tau_image"] != k)
    a, b = chiral6912
    ma = C.tagged_multiset(np.array(reps[a]["B"], np.int64), C.G.CHIRAL, C.UNITS["CHIRAL"])
    mb = C.tagged_multiset(np.array(reps[b]["B"], np.int64), C.G.CHIRAL, C.UNITS["CHIRAL"])
    assert ma != C.conj(ma) and mb == C.conj(ma)
    assert {C.M.chirality(C.M.spectrum(np.array(reps[k]["B"], np.int64))) for k in (a, b)} == {18, -18}


def test_frozen():
    d = _data()
    assert d["orbitals"] == 20 and d["tau_swapped_equals_reversal"]
    assert d["pair_256_members_checked"] == 512 and d["chi_equals_54Q"] and d["imag_equals_Q"]
    assert not d["pattern_6912_constant_on_orbitals"] and d["pattern_6912_nonzero_on_achiral_orbital"]
    assert d["tagged_constant_on_orbitals"]
    assert all(d["tagged_chiral_iff_tau_swapped"].values())
    assert d["tagged_chiral_on"]["CHIRAL"] == [256, 256, 6912, 6912]
    assert d["tagged_classifies"]["CHIRAL"] == 20
