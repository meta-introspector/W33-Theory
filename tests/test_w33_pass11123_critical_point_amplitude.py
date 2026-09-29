"""Regression for Pass 11123: critical radii sqrt3 and phi^2/sqrt3; no massless vectors; no contact quartic."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11123_critical_point_amplitude as P  # noqa: E402


def test_critical_point():
    r = P.summarize()
    radii = sorted(round(v, 6) for v in r['critical_radii'].values())
    assert radii == [round(r['phi2_over_sqrt3'], 6)] * 4 + [round(r['sqrt3'], 6)] * 2
    assert not r['massless_vector_from_pairs'] and r['t_u_symmetric']
    assert r['amplitude_lead'] == '2*a/(b*(a + b))' and r['amplitude_const'] == '-3*a**2/(b*(a + b))'
    f, n = r['klt_formula_vs_numeric']
    assert abs(f - n) < 1e-6
    assert all(v > 0.2 for v in r['radion_slopes'].values())
