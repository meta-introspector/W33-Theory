"""Regression for Pass 11182: the paper's ticks under quantum mereology (frozen profiles)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen():
    d = json.loads((ROOT / "data" / "w33_pass11182_paper_ticks_mereology.json").read_text())
    assert d['factorisations'] == 110565
    assert d['ticks']['clock_K'] == {'local': 108, 'perfect': 8748, 'partial': 101709, 'order': 3}
    assert d['ticks']['perfect_tick_VKV']['local'] == 0 and d['ticks']['perfect_tick_VKV']['order'] == 9
    assert d['ticks']['f9_gate']['local'] == 27
    c = d['clock_splits']
    assert c['all_planes_fixed'] == 108 and c['containing_lightcone'] == 0 and c['containing_transverse'] == 3
    assert c['position_momentum_aligned'] == 4 and c['aligned_q_orthogonal'] == 4 and c['q_orthogonal_frames'] == 4
