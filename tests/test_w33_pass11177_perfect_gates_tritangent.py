"""Regression for Pass 11177: perfect two-qutrit gates <=> collinear factorisations in GQ(4,2) (tritangent planes)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11177_perfect_gates_tritangent as P  # noqa: E402


def test_dictionary():
    r = P.summarize()
    assert r['factorisations'] == 45 and r['degree'] == [12] and r['srg_lambda'] == [3] and r['srg_mu'] == [3]
    assert r['frames'] == 27 and r['frames_partition_40'] and r['pairs_per_frame_count'] == {1: 270}
    assert r['census'] == {'02|equal': 576, '12|noncollinear': 18432, '20|equal': 576, '21|noncollinear': 18432,
                           '22|collinear': 13824}
    assert r['H_order'] == 1152 and r['H_orbit_on_collinear'] == 12 and r['H_stabiliser_of_neighbour'] == 96
    assert r['gap_rank3_45'] and r['gap_rank20_110565']
