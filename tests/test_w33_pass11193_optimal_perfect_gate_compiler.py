"""Regression for Pass 11193's optimal W33 perfect-gate compiler."""
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRODUCER = ROOT / "analysis/w33_pass11193_optimal_perfect_gate_compiler.py"
CERT = ROOT / "data/w33_pass11193_optimal_perfect_gate_compiler.json"


def load():
    return json.loads(CERT.read_text(encoding="utf-8"))


def test_producer_replays():
    subprocess.run([sys.executable, str(PRODUCER)], cwd=ROOT, check=True,
                   capture_output=True, text=True, timeout=120)


def test_optimal_double_coset_normal_form():
    x = load()["normal_form"]
    assert x["optimal_depth_counts"] == {"0": 1152, "1": 13824, "2": 36864}
    assert x["optimal_projective_depth_counts"] == {"0": 576, "1": 6912, "2": 18432}
    assert x["all_group_elements_reconstructed"] == 51840
    assert x["maximum_perfect_gate_depth"] == 2
    assert x["mean_perfect_gate_depth"] == "76/45"
    assert len(x["depth_table_sha256"]) == 64


def test_hecke_fusion_and_sum_word():
    x = load()
    g = x["hecke_fusion"]["graph"]
    assert g["identity"] == "A^2 = 12 I + 3 A + 3 (J-I-A)"
    assert g["spectrum"] == {"-3": 24, "3": 20, "12": 1}
    assert g["uniform_two_step_probabilities"] == {
        "local": "1/12", "partial": "2/3", "perfect": "1/4"
    }
    assert x["hecke_fusion"]["p_H_p_middle_factor_counts"] == {
        "local": 96, "partial": 768, "perfect": 288
    }
    assert x["SUM_compilation"]["optimal_perfect_gate_depth"] == 2
    assert x["SUM_compilation"]["product_verified"] is True


def test_e6_yukawa_sector_walk_is_exactly_lumpable():
    x = load()["E6_yukawa_sector_walk"]
    assert x["sector_sizes"] == {"1.10.10": 5, "16.16.10": 40}
    assert x["perfect_neighbour_quotient"] == [[4, 8], [1, 11]]
    assert x["transition_matrix"] == [["1/3", "2/3"], ["1/12", "11/12"]]
    assert x["stationary_distribution"] == ["1/9", "8/9"]
    assert x["eigenvalues"] == ["1", "1/4"]
    assert "no Yukawa amplitudes" in x["scope"]


def test_universality_boundary_and_visible_surfaces():
    x = load()["universal_computation_boundary"]
    assert "Clifford" in x["boundary"]
    assert x["status_with_existing_E6_cubic_and_Fourier"].startswith("APPROXIMATELY_UNIVERSAL")
    docs = (ROOT / "docs/index.html").read_text(encoding="utf-8")
    tail = (ROOT / "analysis/W33_SHARED_FRONTIER_TAIL.tex").read_text(encoding="utf-8")
    assert docs.count('id="pass11193-perfect-gate-compiler"') == 1
    assert tail.count("PASS11193_OPTIMAL_W33_PERFECT_GATE_COMPILER_INSERT") == 1
