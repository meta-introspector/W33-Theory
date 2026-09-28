"""Regression for Pass 11094: every theta-sector tachyon of the 87 Z6-I twins is coloured or has charge +-1/2."""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11094_tachyons_break_electromagnetism as P  # noqa: E402

C = json.loads(P.OUT.read_text())
S = C["summary"]


def test_no_tachyon_is_neutral_and_colourless():
    assert S["models"] == 87 and S["python_neutral_colourless_tachyonic_momenta"] == 0
    assert set(S["python_colourless_charges"]) == {"1/2", "-1/2"}
    assert S["python_models_with_coloured_tachyon"] == 27


def test_two_engines_agree():
    assert S["orbifolder_models_every_tachyon_field_sm_charged"] == 87
    assert S["orbifolder_models_every_tachyon_field_without_neutral_component"] == 87
    assert set(S["orbifolder_tachyon_reps"]) == {"(1,1)_1/2", "(1,1)_-1/2", "(1,2)_0", "(3,1)_1/6", "(-3,1)_1/6",
                                                  "(3,1)_-1/6", "(-3,1)_-1/6"}


def test_recompute_one_model_weight_by_weight():
    orb = json.loads(P.ORB.read_text())["17"]
    m = {f"{x['index']:02d}": x for x in json.loads(P.P.MODELS.read_text())["models"]}["17"]
    tY = [F(x) for x in orb["tY"]]
    col = [[F(x) for x in r] for r in orb["colour_roots"]]
    su2 = [[F(x) for x in r] for r in orb["su2_roots"]]
    st = [P.classify(p, tY, col, su2) for _, _, p in P.tachyon_momenta([F(x) for x in m["V"]], [F(x) for x in m["W5"]])]
    assert len(st) == C["models"]["17"]["tachyonic_momenta"] and not any(s["neutral_colourless"] for s in st)
