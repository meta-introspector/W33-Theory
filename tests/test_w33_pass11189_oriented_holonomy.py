"""Regression for Pass 11189: oriented holonomy and the symplectic twist classify all 20 three-qutrit split
relations including their time direction (frozen certificate; the labels need the ~25 minute orbit computation)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11189_oriented_holonomy.json").read_text())
    assert d["holonomy_plus_twist_classifies_all_20"] and not d["holonomy_classifies_all_20"]
    rp = d["reversal_pairs"]
    assert rp["12-13"]["separated"] and rp["12-13"]["shortest_walk_half_length"] == 2
    assert not rp["3-4"]["separated"] and rp["3-4"]["twist"]["separated"]
    tau = {o["orbital"]: o["tau_image"] for o in d["orbitals"]}
    rev = {o["orbital"]: o["reverse"] for o in d["orbitals"]}
    assert tau[3] == 4 and tau[12] == 13 and tau[7] == 7 and tau[9] == 9
    assert all(tau[o] == o for o in tau if rev[o] == o)


def test_sl_class_separates_n_from_minus_n():
    import w33_pass11189_oriented_holonomy as H
    N = np.array([[0, 1], [0, 0]])
    assert H.slclass(N) != H.slclass((-N) % 3)
    u = np.array([[1, 1], [0, 1]])
    assert H.slclass(u) != H.slclass(np.array([[1, 2], [0, 1]]))
