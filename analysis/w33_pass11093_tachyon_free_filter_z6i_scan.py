#!/usr/bin/env python3
"""Pass 11093: tachyon-freedom and the Standard Model exclude each other in the Z6-I W(3,3) scan.

Pass 11092: a fixed point of the non-SUSY Z6-I twin (twist (0,1/6,1/6,2/3)) is tachyonic iff the SUSY parent has a
massless oscillator-excited left mover there.  Applied here to the whole Z6-I W(3,3) scan (58 base shifts of the
Holotrade 5b3f3ad scan, random Wilson lines W5 = W6):

  A. exact (pure Python): the fixed point WITHOUT Wilson line has local shift V alone.  For 47 of the 58 base shifts it is
     already tachyonic, so EVERY model built on them is tachyonic whatever its Wilson lines;
  B. frozen C++ evidence (patched orbifolder, analysis/orbifolder_n0_drivers/nsscan.cpp): 6000 random Wilson-line draws on
     each of the 11 surviving base shifts.  The lemma is applied in C++ to the SUSY parent (theta-sector states with a
     nonzero oscillator number); the C++ test reproduces the Python verdict 87/87 on the known SM parents and the
     Font--Hernandez tachyon-free control.  Every tachyon-free twin is then built and put through the SM test;
  C. which parent states make the SM parents tachyonic: in all 705 SM parents the oscillator-excited theta-sector states
     are the exotics labelled v, w, x (and conjugates) -- never quark, lepton or Higgs fields.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11092_nonsusy_z6i_twins_are_tachyonic as P  # noqa: E402

BASES = ROOT / "data" / "w33_pass11093_z6i_base_shifts.json"
EVID = ROOT / "data" / "w33_pass11093_scan_evidence.json"
OUT = ROOT / "data" / "w33_pass11093_tachyon_free_filter_z6i_scan.json"


def base_fixed_point(V):
    """theta-sector fixed point without Wilson line: (tachyonic left states, massless left states, min 1/2 p_sh^2)"""
    osc = P.osc_counts(P.OSC_THETA, F(3, 4))
    lv = P.half_norm_counts(V, F(3, 4))
    tach = sum(c * osc[F(7, 12) - h] for h, c in lv.items() if F(7, 12) - h in osc)
    mass = sum(c * osc[F(3, 4) - h] for h, c in lv.items() if F(3, 4) - h in osc)
    return dict(tachyonic_left_states=tach, massless_left_states=mass, min_half_norm=str(min(lv)) if lv else None)


def poisson_zero(expected):
    return math.exp(-expected)


def main():
    bases = json.loads(BASES.read_text())["models"]
    A = {}
    for b in bases:
        A[b["label"]] = base_fixed_point([F(x) for x in b["V"]])
        print(b["label"], A[b["label"]], flush=True)
    survivors = sorted(k for k, v in A.items() if v["tachyonic_left_states"] == 0)
    ev = json.loads(EVID.read_text())
    scan = ev["scan"]
    tot = {k: sum(v[k] for v in scan.values()) for k in next(iter(scan.values()))}
    # expected overlap of SM parents and tachyon-free models if the two were independent (per base)
    expected = sum(v["parent_sm"] * v["tachyon_free"] / v["ok"] for v in scan.values())
    summary = dict(
        base_shifts=len(bases), tachyonic_at_wilson_line_free_fixed_point=len(bases) - len(survivors),
        surviving_base_shifts=survivors,
        scan_totals=tot,
        sm_parents_among_tachyon_free=tot["tachyon_free_with_sm_parent"],
        expected_overlap_if_independent=round(expected, 2),
        poisson_probability_of_zero=poisson_zero(expected),
        sm_parent_excited_theta_labels=ev["scan2_parent_sm_excited_labels"],
        sm_parents_checked_for_excited_states=ev["scan2_parent_sm_models"],
        sm_parents_without_excited_theta_state=ev["scan2_parent_sm_models_without_excited_theta_state"],
        tachyon_free_inequivalent_spectra=ev["tachyon_free_inequivalent_models"],
        tachyon_free_with_an_su3_and_an_su2_factor=ev["tachyon_free_with_su3_and_su2_factor"],
    )
    assert set(survivors) == {k.split("/")[-1] for k in scan}, (survivors, list(scan))
    OUT.write_text(json.dumps(dict(pass_id=11093, base_fixed_points=A, summary=summary, scan=scan,
                                   tachyon_free_gauge_groups=ev["tachyon_free_gauge_groups"]), indent=1, sort_keys=True))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
