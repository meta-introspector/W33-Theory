#!/usr/bin/env python3
"""Pass 11262: the G26/MUB origin of the line-stabilizer S4 Yukawa layer.

The twelve order-three G26 mirrors are the four qutrit MUBs.  Complex-linear
G26 induces A4 on those four bases; adjoining complex conjugation supplies one
odd transposition and closes the action to S4.  This gives an object-level
bridge to the four tetrahedral vacua used by the repository's S4 TM1 model.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11262_g26_mub_s4_yukawa_bridge.json"
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11269_g26_qutrit_dictionary as G26  # noqa: E402
import w33_pass11243_tm1_line_geometry as LG  # noqa: E402
import w33_pass11230_tm1_alignment as TM1  # noqa: E402


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def closure(gens):
    identity = tuple(range(len(gens[0])))
    group = {identity}; frontier = [identity]
    while frontier:
        new = []
        for a in frontier:
            for b in gens:
                c = compose(a, b)
                if c not in group:
                    group.add(c); new.append(c)
        frontier = new
    return group


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def mub_action():
    stab = [v / np.linalg.norm(v) for v in G26.num(G26.STAB)]
    sic = [v / np.linalg.norm(v) for v in G26.num(G26.SIC)]
    components = []
    for i, u in enumerate(stab):
        c = tuple(sorted([i] + [j for j, v in enumerate(stab) if abs(np.vdot(u, v)) < 1e-9]))
        if c not in components:
            components.append(c)
    assert components == [(0,1,2), (3,8,10), (4,6,11), (5,7,9)]

    def ray_index(v):
        overlaps = [abs(np.vdot(q, v)) for q in stab]
        j = int(np.argmax(overlaps)); assert overlaps[j] > 1 - 1e-9
        return j

    def induced(M, anti=False):
        p = []
        for c in components:
            v = stab[c[0]].conj() if anti else stab[c[0]]
            j = ray_index(M @ v)
            p.append(next(k for k, d in enumerate(components) if j in d))
        return tuple(p)

    w = np.exp(2j*np.pi/3)
    refs = [np.eye(3) - 2*np.outer(v, v.conj()) for v in sic]
    refs += [np.eye(3) + (w-1)*np.outer(v, v.conj()) for v in stab]
    linear_gens = sorted(set(induced(M) for M in refs))
    linear = closure(linear_gens)
    conjugation = induced(np.eye(3), anti=True)
    completed = closure(linear_gens + [conjugation])
    assert len(linear) == 12 and all(parity(p) == 1 for p in linear)
    assert conjugation == (0,1,3,2) and parity(conjugation) == -1
    assert len(completed) == 24
    return components, linear_gens, conjugation, linear, completed


def flavon_dictionary(components):
    d4 = np.array(LG.D4, dtype=int)
    chords = []
    for i in range(4):
        for j in range(i+1,4):
            for sign in (1,-1):
                chi = sign*(d4[i]-d4[j])
                dots = [int(d @ chi) for d in d4]
                chords.append({"mub_pair": [i,j], "orientation": sign, "chi": chi.tolist(), "endpoint_dots": dots})
    assert len(chords) == 12
    assert sorted({abs(r["endpoint_dots"][r["mub_pair"][0]]) for r in chords}) == [4]
    return {
        "four_phi_axes": [d.tolist() for d in d4],
        "four_MUB_components": [list(c) for c in components],
        "twelve_oriented_chords": chords,
        "rule": "phi axis i labels MUB i; ±(d_i-d_j) labels an oriented pair of MUBs",
    }


def yukawa_audit():
    p54 = json.loads((ROOT / "data/w33_pass11254_tm1_yukawa_breaking.json").read_text())
    p67 = json.loads((ROOT / "data/w33_pass11267_tm1_selector.json").read_text())
    curves = []
    examples = p54["examples"]
    for row_i, eps in enumerate([r["eps"] for r in examples[0]["scan"]]):
        rows = [m["scan"][row_i] for m in examples]
        curves.append({
            "epsilon": eps,
            "sin2_theta13_three_fixed_models": [r["s13"] for r in rows],
            "cos_delta_three_fixed_models": [r["cos_delta"] for r in rows],
            "TM1_column_breaking_dUe1": [r["dUe1"] for r in rows],
        })
    bosons = [r for r in p67["loop"]["rows"] if r["statistics"] == "boson"]
    natural = [r for r in bosons if abs(r["eps_eff"]) <= 1e-3]
    return {
        "epsilon_does_not_determine_angles": True,
        "evidence_at_epsilon_zero": {
            "sin2_theta23_range_over_60_TM1_fits": p54["eps0_s23_range"],
            "cos_delta_range_over_60_TM1_fits": p54["eps0_cos_delta_range"],
        },
        "three_fixed_model_curves": curves,
        "TM1_alignment_window": "|epsilon| approximately <= 1e-3",
        "bosonic_loop_rows_inside_window": natural,
        "interpretation": (
            "epsilon controls departure from the fixed TM1 column. theta13 and delta additionally depend on "
            "Yukawa coefficients. Odd MUB permutations require the antiunitary completion, so the CP phase "
            "cannot be assigned by the complex-linear G26 vacuum coordinate alone."
        ),
    }


def run():
    components, gens, conj, linear, completed = mub_action()
    out = {
        "schema": "w33.pass11262.g26-mub-s4-yukawa-bridge.v1",
        "status": "PASS_S4_IS_ANTIUNITARY_COMPLETION_OF_G26_A4_MUB_ACTION_AND_YUKAWA_SCOPE_AUDITED",
        "mub_action": {
            "components": [list(c) for c in components],
            "complex_linear_image_order": len(linear),
            "complex_linear_image": "A4",
            "all_linear_permutations_even": all(parity(p) == 1 for p in linear),
            "complex_conjugation_permutation": list(conj),
            "conjugation_is_odd": parity(conj) == -1,
            "completed_image_order": len(completed),
            "completed_image": "S4",
        },
        "flavon_dictionary": flavon_dictionary(components),
        "yukawa_audit": yukawa_audit(),
        "boundary": (
            "The four-MUB/tetrahedron dictionary is canonical up to relabeling. It identifies the symmetry "
            "objects and the role of antiunitarity; it does not identify arbitrary Cartan amplitudes with "
            "real flavon expectation values or remove the independent Yukawa coefficients."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
