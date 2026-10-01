import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis" / "w33_20261001_g26_e8_mirror_separation.py"
SPEC = importlib.util.spec_from_file_location("mirror_sep", SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)


def test_payload_replays_frozen_certificate():
    got = M.payload()
    frozen = json.loads((ROOT / "data/w33_20261001_g26_e8_mirror_separation.json").read_text())
    assert got == frozen


def test_divisor_lattice_index_and_isolation():
    d = M.payload()
    assert d["valuation_matrix_determinant"] == 5
    assert d["divisor_isolation"]["2P_minus_J"] == [5, 0]
    assert d["divisor_isolation"]["3J_minus_P"] == [0, 5]
    ids = d["normalized_polynomial_identities"]
    assert ids["SIC_fifth_power"]["degree_each_side"] == 78
    assert ids["stabilizer_fifth_power"]["degree_each_side"] == 99
    assert all(d["checks"].values())
