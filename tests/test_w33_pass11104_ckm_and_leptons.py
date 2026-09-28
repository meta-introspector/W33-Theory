"""Regression for Pass 11104: CKM texture and charged leptons of the 12 survivors."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11104_ckm_and_leptons as P  # noqa: E402


def test_alignment_and_textures():
    r = P.summarize()
    assert r["models"] == 12
    assert r["pairs"] == {"same_doublet=True|tau_heavy=True": 36, "same_doublet=False|tau_heavy=True": 72}
    assert r["each_Hu_has_one_aligned_Hd"]
    for g, tab in r["ckm_exponents_aligned"].items():
        for k in tab:
            E = eval(k)
            assert E[0][0] == 0 and E[0][1] == 2 and E[0][2] == 2 and E[1][0] == 2 and E[2][0] == 2   # V_tb ~ 1, rest ~ eps^2
    for k in r["lepton_exponents_aligned"]:
        for e in eval(k).values():
            assert e[0] == 0 and e[1] == e[2] > 0            # single heavy tau, mu and e at the same order
    assert r["tree_level_neutrino_dirac"] == 12
