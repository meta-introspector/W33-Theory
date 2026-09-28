#!/usr/bin/env python3
"""Pass 11095: an independent code checks the N=0 engine, and the W(3,3) A8 shift has tachyon-free non-supersymmetric
three-generation completions in the SO(16)xSO(16) string.

A. Cross-check against the non-SUSY orbifolder (Escalante-Notario, Perez-Martinez, Ramos-Sanchez, Vaudrevange,
   arXiv:2504.20137, github.com/StringsIFUNAM/nonSUSYorbifolder).  Its 19 sample models (Z2W x Z_N, all Z_N geometries it
   ships) are converted to orbifolder-1.2.1 format and run through our patched engine: massless spectra IDENTICAL,
   19/19, after PATCH 6 (order of a twist on spinors).  Their code already contains the right-moving oscillator sign fix
   ("//hacking here!!! wrong trafo sign for R-moving oscillator excitations in original Orbifolder"): our PATCH 3 is a
   rediscovery and is credited to them.  Their code cannot load the E8xE8 twins of Passes 11092-11096 (its right-mover
   lattice is hard-coded for the Witten-twist structure).  Tachyon control: SO(16)xE8 on T6/Z3 (standard embedding) --
   both codes find exactly one tachyon, (10,1,1).
B. SO(16)xSO(16) completions of the W(3,3) A8 Kac pair V1 = (1/6^7, 5/6; 0^7, 2/3) of T6/Z3: Witten shifts
   V0 = (a; b), a, b in (1/2)(norm-4 vectors of E8) (one Weyl orbit each, 2160), with V0.V1 in Z (exact enumeration here).
C. Frozen scan (non-SUSY orbifolder, analysis/orbifolder_n0_drivers/nsoscan.cpp): random Wilson lines on one V0 per
   loadable gauge class, tachyon test and its SM-like test; every SM-like model re-verified with OUR engine: tachyons at
   the only tachyonic level (M^2/8 = -1/2, Witten sector), anomalies, and the net chiral content with THEIR hypercharge
   applied to OUR fermion spectrum.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
XCHECK = ROOT / "data" / "w33_pass11095_nonsusy_orbifolder_crosscheck.json"
SCAN = ROOT / "data" / "w33_pass11095_a8_scan_evidence.json"
OUT = ROOT / "data" / "w33_pass11095_so16xso16_a8_standard_models.json"

V1A = (F(1, 6),) * 7 + (F(5, 6),)
V1B = (F(0),) * 7 + (F(2, 3),)


def e8_norm(n2):
    out = []
    for half in (False, True):
        rng = [F(k, 2) for k in range(-5, 6, 2)] if half else [F(k) for k in range(-2, 3)]
        for v in itertools.product(rng, repeat=8):
            if sum(x * x for x in v) == n2 and sum(v) % 2 == 0:
                out.append(v)
    return out


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def root_components(a, V, roots):
    kept = [r for r in roots if dot(r, a).denominator == 1 and dot(r, V).denominator == 1]
    seen, comps = set(), []
    for r in kept:
        if r in seen:
            continue
        stack, n = [r], 0
        seen.add(r)
        while stack:
            x = stack.pop()
            n += 1
            for y in kept:
                if y not in seen and dot(x, y) != 0:
                    seen.add(y)
                    stack.append(y)
        comps.append(n)
    names = {2: "A1", 6: "A2", 12: "A3", 20: "A4", 30: "A5", 42: "A6", 56: "A7", 24: "D4", 84: "D7", 112: "D8"}
    return " x ".join(names.get(c, f"[{c}]") for c in sorted(comps, reverse=True))


def witten_shift_classes():
    roots = e8_norm(2)
    half = [tuple(x / 2 for x in v) for v in e8_norm(4)]
    assert len(roots) == 240 and len(half) == 2160
    fa = Counter(dot(a, V1A) % 1 for a in half)
    fb = Counter(dot(b, V1B) % 1 for b in half)
    admissible = sum(fa[x] * fb[(-x) % 1] for x in fa)
    ga = {}
    for a in half:
        ga.setdefault(dot(a, V1A) % 1, Counter())[root_components(a, V1A, roots)] += 1
    gb = {}
    for b in half:
        gb.setdefault(dot(b, V1B) % 1, Counter())[root_components(b, V1B, roots)] += 1
    classes = Counter()
    for x, ca in ga.items():
        for na, ka in ca.items():
            for nb, kb in gb.get((-x) % 1, {}).items():
                classes[f"{na} | {nb}"] += ka * kb
    return dict(half_norm4_vectors=len(half), admissible_V0=admissible, total_pairs=len(half) ** 2,
                gauge_classes=dict(classes))


def main():
    B = witten_shift_classes()
    print(json.dumps(B, indent=1), flush=True)
    xc = json.loads(XCHECK.read_text())
    sc = json.loads(SCAN.read_text())
    ver = sc["verification"]
    sm = {"(3,2)_1/6": 3, "(-3,1)_-2/3": 3, "(-3,1)_1/3": 3, "(1,2)_-1/2": 3, "(1,1)_1": 3}
    sm_c = {"(-3,2)_1/6": 3, "(3,1)_-2/3": 3, "(3,1)_1/3": 3, "(1,2)_-1/2": 3, "(1,1)_1": 3}   # 3 <-> 3bar relabelled
    three_gen = sum(1 for v in ver.values() if v["net_chiral"] in (sm, sm_c))
    summary = dict(
        crosscheck_models=len(xc["models"]), crosscheck_identical=sum(v["gauge_same"] and v["spectrum_same"] for v in xc["models"].values()),
        tachyon_control=xc["tachyon_control"],
        witten_shifts=B,
        loadable_classes=sc["loadable_classes"], rejected_classes=sc["rejected_classes"],
        scan=sc["scan"],
        sm_like_models=len(ver),
        verified_three_generations=three_gen,
        verified_no_net_fractional_colourless=sum(1 for v in ver.values() if not v["net_fractional_colourless"]),
        verified_with_higgs_scalar=sum(1 for v in ver.values() if v["higgs_doublet_scalars"] > 0),
        our_tachyonic_fields_at_minus_half=sc["our_tachyonic_fields_at_minus_half"],
        their_sm_and_tachyon_free=sc["their_sm_and_tachyon_free"],
        our_anomaly=sc["our_anomaly"],
        higgs_scalar_fields=dict(Counter(v["higgs_doublet_scalars"] for v in ver.values())),
        # the open problem: vector-like but present -- every model carries fractionally charged colourless fermions
        vectorlike_fractional_fermion_states_min=min(v["fractional_colourless_fermion_states"] for v in ver.values()),
        vectorlike_fractional_fermion_states_max=max(v["fractional_colourless_fermion_states"] for v in ver.values()),
        models_without_any_fractional_fermion=sum(1 for v in ver.values() if v["fractional_colourless_fermion_states"] == 0),
    )
    OUT.write_text(json.dumps(dict(pass_id=11095, summary=summary, verification=ver), indent=1, sort_keys=True, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
