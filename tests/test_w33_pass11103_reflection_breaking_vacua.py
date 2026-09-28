"""Regression for Pass 11103: no alignment of the SM-neutral scalar VEVs splits charm from up in the 12 survivors."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11103_reflection_breaking_vacua as P  # noqa: E402


def test_up_never_splits_down_splits_only_through_vectorlike_dbar():
    res = P.summarize()["results"]
    for var in ("A", "B"):
        r = res[var]
        assert r["models"] == 12 and r["up_split_cases"] == 0
        assert r["down_split_models"] == r["models_with_6_dbar"] and len(r["down_split_models"]) == 8
        # VEVs aligned at the Higgs point never split anything (Schur)
        for g in ("1.0", "4.0"):
            assert set(r["table"][f"g={g}|only_p0"]) == {"up|gap=0.0", "down|gap=0.0"}
        # small torus area: nothing splits for any alignment
        for w, _ in P.GRID:
            assert set(r["table"][f"g=1.0|{w}"]) == {"up|gap=0.0", "down|gap=0.0"}


def test_frozen_rows_complete():
    for path in P.FROZEN.values():
        d = json.loads(path.read_text())
        assert sorted(map(int, d)) == [2, 10, 13, 14, 15, 35, 53, 57, 69, 77, 78, 102]
        assert all(len(r["res"]["up"]) == 3 and len(r["res"]["down"]) == 3 for r in d.values())
