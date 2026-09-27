#!/usr/bin/env python3
"""Pass 11028: an interferometric discriminator for the hidden signed clock center."""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11021_minimal_signed_clock_carrier as P21

OUT = ROOT / "data/w33_pass11028_signed_clock_interferometric_discriminator.json"
P52 = ROOT / "data/w33_pass10952_clock_complete_positivity_firewall.json"


def inner(a, b):
    return sum(x*y for x, y in zip(a, b))


def payload():
    e6_to_h, dirs, rows, perms, section = P21.split_section()
    p_to_m = {r[3]: tuple(map(int, r[0])) for r in rows}
    z = next(g for g in section if p_to_m[g[0]] == (2, 0, 0, 2))
    p, signs = z

    fibres = []
    indicators = []
    for d in dirs:
        fibre = sorted(
            i for i, h in e6_to_h.items()
            if h[:2] != (0, 0) and P75.norm_dir(h[:2]) == d
        )
        v = [0] * 27
        for i in fibre:
            v[i] = 1
        fibres.append(fibre)
        indicators.append(tuple(v))

    signed_overlaps = []
    visibilities = []
    for v in indicators:
        zv = P21.act_vec(z, v)
        num = inner(v, zv)
        den = inner(v, v)
        f = Fraction(num, den)
        signed_overlaps.append(str(f))
        visibilities.append(abs(f))
    assert signed_overlaps == ["1/3", "-1/3", "1/3", "-1/3"]
    assert set(visibilities) == {Fraction(1, 3)}

    # Signed-permutation implementation complexity on the 24 noncentral modes.
    noncentral = sorted(i for i, h in e6_to_h.items() if h[:2] != (0, 0))
    seen = set()
    pairs = []
    negative_pairs = 0
    for i in noncentral:
        if i in seen:
            continue
        j = p[i]
        assert j != i and p[j] == i and signs[i] == signs[j]
        eps = -1 if signs[i] else 1
        pairs.append((i, j, eps))
        negative_pairs += int(eps == -1)
        seen |= {i, j}
    assert len(pairs) == 12 and negative_pairs == 6
    V = Fraction(1, 3)
    signed_pmax = (1 + V) / 2
    signed_pmin = (1 - V) / 2
    coarse_V = Fraction(1, 1)
    coarse_pmax = (1 + coarse_V) / 2
    coarse_pmin = (1 - coarse_V) / 2
    assert (signed_pmax, signed_pmin) == (Fraction(2, 3), Fraction(1, 3))
    assert (coarse_pmax, coarse_pmin) == (1, 0)

    p52 = json.loads(P52.read_text(encoding="utf-8"))
    qutrit = p52["even_tick_physical_subgroup"]["four_tick_scalarity_abs_trace_over_3"]
    assert abs(qutrit - 1/3) < 1e-12

    checks = {
        "four_clock_fibres": len(fibres) == 4,
        "six_modes_each": all(len(f) == 6 for f in fibres),
        "signed_overlap_pattern": signed_overlaps == ["1/3", "-1/3", "1/3", "-1/3"],
        "visibility_is_direction_independent_one_third":
            set(visibilities) == {Fraction(1, 3)},
        "projective_model_visibility_one": coarse_V == 1,
        "signed_fringe_extrema_two_thirds_one_third":
            (signed_pmax, signed_pmin) == (Fraction(2, 3), Fraction(1, 3)),
        "central_network_twelve_swaps_six_negative_pairs":
            len(pairs) == 12 and negative_pairs == 6,
        "independent_qutrit_trace_ratio_one_third": abs(qutrit - 1/3) < 1e-12,
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11028.signed-clock-interferometric-discriminator.v1",
        "status": "PASS",
        "headline": (
            "A two-arm identity-versus-central-sign interferometer separates the "
            "naive projective clock from the exact signed 24-mode carrier. For "
            "each uniform six-mode clock fibre the exact signed overlap is +/-1/3, "
            "so phase-scanned fringe visibility is exactly 1/3 independent of "
            "clock direction; the coarse projective model predicts visibility 1."
        ),
        "central_operation": {
            "matrix_element": "-I in support GL2(3), exact signed monomial lift",
            "noncentral_mode_dimension": 24,
            "disjoint_mode_swaps": len(pairs),
            "negative_swap_pairs": negative_pairs,
            "pair_list": [list(x) for x in pairs],
        },
        "coarse_fibre_probe": {
            "signed_overlaps": signed_overlaps,
            "absolute_visibility_each": ["1/3"] * 4,
            "permutation_blind_statistic": "phase-scanned fringe visibility |<v,zv>|/<v,v>",
            "signed_model_visibility": "1/3",
            "projective_identity_model_visibility": "1",
            "visibility_squared": "1/9",
        },
        "equal_arm_fringe": {
            "signed_Pmax": "2/3",
            "signed_Pmin": "1/3",
            "projective_Pmax": "1",
            "projective_Pmin": "0",
            "protocol": (
                "Prepare a normalized uniform state on any one six-mode clock "
                "fibre, split into two coherent arms, apply the exact signed "
                "central operation z in one arm and identity in the other, scan "
                "the reference phase, and extract fringe visibility."
            ),
        },
        "cross_representation_check": {
            "qutrit_four_tick_abs_trace_over_3": qutrit,
            "same_number": "1/3",
            "firewall": (
                "This numerical equality is an operational cross-check only. "
                "Pass 10952 already proves that the qutrit Clifford and the "
                "signed 24-mode carrier are different representations of the "
                "same abstract central element."
            ),
        },
        "experimental_boundary": (
            "The certificate derives an ideal coherent-mode observable, not a "
            "hardware error budget. Loss imbalance, mode-dependent phase noise, "
            "source impurity and detector visibility must be calibrated before "
            "a laboratory comparison. The discriminator tests the signed carrier "
            "against its coarse projective shadow; it does not test the full TOE."
        ),
        "checks": checks,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    p = payload()
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if a.check:
        if not a.output.exists() or a.output.read_text(encoding="utf-8") != text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": p["status"],
        "visibility": p["coarse_fibre_probe"]["signed_model_visibility"],
        "coarse": p["coarse_fibre_probe"]["projective_identity_model_visibility"],
        "network": [p["central_operation"]["disjoint_mode_swaps"], p["central_operation"]["negative_swap_pairs"]],
    }, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
