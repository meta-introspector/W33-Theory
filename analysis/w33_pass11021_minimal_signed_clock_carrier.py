#!/usr/bin/env python3
"""Pass 11021: minimal exact signed clock carrier.

Using the deterministic split section from Pass 11020, close the four coarse
clock fibre indicators and the three clock augmentation generators under the
exact signed GL2(3) action. Both closures have dimension 24: precisely the
noncentral H27 coordinates.

The support orbits refine 27 as 1+2+8+16, with each six-point clock fibre
splitting 2+4 across the noncentral 8- and 16-orbits.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10975_e6_cubic_reynolds_clock_weld as P75
import w33_pass11020_signed_clock_extension_split as P78

OUT = ROOT / "data" / "w33_pass11021_minimal_signed_clock_carrier.json"
P11020 = ROOT / "data" / "w33_pass11020_signed_clock_extension_split.json"
def act_vec(element, v):
    p, s = element
    out = [0] * len(v)
    for i, a in enumerate(v):
        if a:
            out[p[i]] += (-1 if s[i] else 1) * a
    return tuple(out)


def orbit_span_rank(section, vectors):
    images = []
    orbit_sizes = []
    for v in vectors:
        orbit = {act_vec(g, v) for g in section}
        orbit_sizes.append(len(orbit))
        images.extend(orbit)
    matrix = sp.Matrix.hstack(*[sp.Matrix(v) for v in images])
    return matrix.rank(), orbit_sizes, len(set(images))


def coordinate_orbits(perms, n):
    unseen = set(range(n))
    out = []
    while unseen:
        i = min(unseen)
        orbit = {p[i] for p in perms}
        out.append(tuple(sorted(orbit)))
        unseen -= orbit
    return tuple(sorted(out, key=lambda x: (len(x), x)))


def mon_trace(element, subset):
    p, s = element
    return sum((-1 if s[i] else 1) for i in subset if p[i] == i)
def split_section():
    (
        signed, e6_to_h, dirs, rows, triads, perms,
        A, kernel_basis, particular,
    ) = P78.signed_setup()
    g0, g1 = P78.find_support_generators(perms)
    _, _, section, _, _, _, _ = P78.find_split_pair(
        g0, g1, particular, kernel_basis
    )
    assert len(section) == 48
    return e6_to_h, dirs, rows, perms, section


def fit_quadratic_graph(points):
    fits = []
    for coeff in itertools.product(range(3), repeat=6):
        aa, ab, bb, la, lb, z = coeff
        if all(
            (aa*a*a + ab*a*b + bb*b*b + la*a + lb*b + z - c) % 3 == 0
            for a, b, c in points
        ):
            fits.append(coeff)
    return fits


def payload():
    parent = json.loads(P11020.read_text(encoding="utf-8"))
    assert parent["extension"]["splits"] is True

    e6_to_h, dirs, rows, perms, section = split_section()
    support_orbits = coordinate_orbits(perms, 27)
    orbit_sizes = [len(o) for o in support_orbits]
    assert orbit_sizes == [1, 2, 8, 16]
    center = sorted(i for i, h in e6_to_h.items() if h[:2] == (0, 0))
    noncentral = sorted(i for i, h in e6_to_h.items() if h[:2] != (0, 0))
    assert len(center) == 3 and len(noncentral) == 24

    noncentral_orbits = [o for o in support_orbits if set(o) <= set(noncentral)]
    assert sorted(map(len, noncentral_orbits)) == [8, 16]
    orbit8 = next(o for o in noncentral_orbits if len(o) == 8)
    orbit16 = next(o for o in noncentral_orbits if len(o) == 16)

    fibres = []
    for d in dirs:
        fibre = sorted(
            i for i, h in e6_to_h.items()
            if h[:2] != (0, 0) and P75.norm_dir(h[:2]) == d
        )
        assert len(fibre) == 6
        fibres.append(fibre)

    intersections = [
        [len(set(f) & set(orbit8)), len(set(f) & set(orbit16))]
        for f in fibres
    ]
    assert intersections == [[2, 4]] * 4

    clock_indicators = []
    for fibre in fibres:
        v = [0] * 27
        for i in fibre:
            v[i] = 1
        clock_indicators.append(tuple(v))
    clock_rank, clock_orbit_sizes, clock_image_count = orbit_span_rank(
        section, clock_indicators
    )
    assert clock_rank == 24

    augmentation = [
        tuple(a - b for a, b in zip(clock_indicators[i], clock_indicators[3]))
        for i in range(3)
    ]
    aug_rank, aug_orbit_sizes, aug_image_count = orbit_span_rank(
        section, augmentation
    )
    assert aug_rank == 24

    center_basis = []
    for i in center:
        v = [0] * 27
        v[i] = 1
        center_basis.append(tuple(v))
    center_rank, _, _ = orbit_span_rank(section, center_basis)
    assert center_rank == 3

    # The 8-point support orbit is a unique quadratic phase graph in this gauge.
    orbit8_addresses = [e6_to_h[i] for i in orbit8]
    fits = fit_quadratic_graph(orbit8_addresses)
    assert fits == [(0, 1, 0, 1, 1, 0)]

    chars24 = [mon_trace(g, noncentral) for g in section]
    chars3 = [mon_trace(g, center) for g in section]
    chars27 = [mon_trace(g, range(27)) for g in section]
    trace24 = Counter(chars24)
    trace3 = Counter(chars3)
    trace27 = Counter(chars27)
    inv24 = sum(chars24) // 48
    inv3 = sum(chars3) // 48
    inv27 = sum(chars27) // 48
    comm24 = sum(x*x for x in chars24) // 48
    comm3 = sum(x*x for x in chars3) // 48
    comm27 = sum(x*x for x in chars27) // 48
    assert (inv24, inv3, inv27) == (2, 2, 4)
    assert (comm24, comm3, comm27) == (19, 5, 34)

    checks = {
        "parent_split_section": parent["extension"]["splits"] is True,
        "support_orbits_1_2_8_16": orbit_sizes == [1, 2, 8, 16],
        "coordinate_split_24_plus_3": len(noncentral) == 24 and len(center) == 3,
        "clock_fibres_all_size6": all(len(f) == 6 for f in fibres),
        "each_clock_fibre_is_2_plus_4": intersections == [[2, 4]] * 4,
        "four_clock_indicators_close_to_rank24": clock_rank == 24,
        "augmentation_alone_closes_to_rank24": aug_rank == 24,
        "center_closes_to_rank3": center_rank == 3,
        "unique_quadratic_graph": fits == [(0, 1, 0, 1, 1, 0)],
        "noncentral_invariant_dimension2": inv24 == 2,
        "noncentral_commutant_dimension19": comm24 == 19,
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11021.minimal-signed-clock-carrier.v1",
        "status": "PASS",
        "headline": (
            "Under the exact split signed GL2(3) symmetry, the four coarse clock "
            "fibre indicators and even the three-dimensional clock augmentation "
            "module generate exactly the full 24-dimensional noncentral H27 sector. "
            "The remaining three central coordinates form a separate invariant sector."
        ),
        "coordinate_decomposition": {
            "full_dimension": 27,
            "support_orbit_sizes": orbit_sizes,
            "central_dimension": len(center),
            "noncentral_dimension": len(noncentral),
            "signed_invariant_split": "27 = 24 noncentral + 3 central",
        },
        "clock_fibres": {
            "count": len(fibres),
            "coordinates_per_fibre": 6,
            "noncentral_support_orbits": [8, 16],
            "intersection_with_8_and_16": intersections,
            "law": "each clock direction resolves internally as 6 = 2 + 4",
        },
        "minimal_closure": {
            "four_indicator_orbit_sizes": clock_orbit_sizes,
            "four_indicator_distinct_images": clock_image_count,
            "four_indicator_span_rank": clock_rank,
            "augmentation_generator_orbit_sizes": aug_orbit_sizes,
            "augmentation_distinct_images": aug_image_count,
            "augmentation_span_rank": aug_rank,
            "center_span_rank": center_rank,
        },
        "quadratic_phase_orbit": {
            "orbit_size": len(orbit8),
            "addresses": [list(e6_to_h[i]) for i in orbit8],
            "polynomial_basis": ["a^2", "ab", "b^2", "a", "b", "1"],
            "unique_coefficients_mod3": list(fits[0]),
            "equation": "c = a*b + a + b (mod 3)",
            "complement_noncentral_orbit_size": len(orbit16),
            "interpretation": (
                "one phase lift over each nonzero (a,b) lies on the quadratic graph; "
                "the other two phase lifts form the 16-point orbit"
            ),
        },
        "character_diagnostics": {
            "noncentral24": {
                "trace_histogram": {str(k): int(v) for k, v in sorted(trace24.items())},
                "invariant_dimension": inv24,
                "commutant_dimension": comm24,
            },
            "center3": {
                "trace_histogram": {str(k): int(v) for k, v in sorted(trace3.items())},
                "invariant_dimension": inv3,
                "commutant_dimension": comm3,
            },
            "full27": {
                "trace_histogram": {str(k): int(v) for k, v in sorted(trace27.items())},
                "invariant_dimension": inv27,
                "commutant_dimension": comm27,
            },
        },
        "what_changed": (
            "The exact signed clock field need not be all 27 cubic coordinates, but "
            "it cannot be only four or three amplitudes either. The minimal invariant "
            "carrier generated by the clock quotient is canonically the 24 noncentral "
            "H27 coordinates."
        ),
        "boundary": (
            "This is a finite representation-closure theorem in the current signed "
            "H27 gauge. It does not identify the 24 coordinates with physical particle "
            "species, spacetime dimensions, or measured degrees of freedom, and it "
            "does not derive a coarse-graining dynamics."
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
        "support_orbits": p["coordinate_decomposition"]["support_orbit_sizes"],
        "clock_span": p["minimal_closure"]["four_indicator_span_rank"],
        "augmentation_span": p["minimal_closure"]["augmentation_span_rank"],
        "phase_graph": p["quadratic_phase_orbit"]["equation"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
