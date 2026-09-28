"""Regression for Pass 11099: positive 10D SO(16)xSO(16) vacuum energy; the A8 models run away at large volume."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11099_lambda10_so16xso16 as L  # noqa: E402
import w33_pass11099_vacuum_energy_so16_a8 as P  # noqa: E402

C = json.loads(P.OUT.read_text())


def test_massless_term_is_n_boson_minus_n_fermion():
    assert L.exact_level_matched()[0] == -2112 == 8 * (8 + 240) - (8 * 256 + 8 * 256)


def test_integrand_modular_invariant_and_integral_negative():
    tau = complex(0.07, 1.21)
    z = L.Z(tau)
    assert abs(L.Z(-1 / tau) - z) / abs(z) < 1e-10 and abs(L.Z(tau + 1) - z) / abs(z) < 1e-10
    up, low, _ = L.integrate(80, 40)
    assert abs((up + low) - C["lambda10"]["I"]) < 1e-6 and up + low < 0          # Lambda_10 = -I/2 > 0


def test_large_volume_weight_and_bose_fermi():
    assert P.volume_sector_weight() == C["large_volume"] and C["large_volume"]["weight_of_so16_on_T6"] == "1/3"
    bf = C["massless_bose_fermi"]
    assert bf["models"] == 104 and bf["degenerate"] == 0 and bf["nB_minus_nF_max"] < 0
