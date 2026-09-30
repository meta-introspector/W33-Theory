#!/usr/bin/env python3
"""Pass 11174: the orientation pattern forced on every perfect interacting tick (Pass 11168) is the image of the
TIME-REVERSING swap of GL(2,3), and a single oblique Hesse clock, used as a potential, already makes the tick perfect.

(1) The Lorentz group of the finite light cone comes from the tetracode clock group.  The 27 events are Sym^2(F_3^2)
    (v = (x^2, xy, y^2), Pass 10946), and g in GL(2,3) acts by Sym^2(g), which preserves q = v0 v2 - v1^2 and has
    det Sym^2(g) = det g.  The 48 elements give 24 of the 48 Lorentz maps (g and -g agree; the other 24 are their
    negatives): 12 of det +1 (from SL(2,3)) and 12 of det -1 (the other coset, PGL(2,3) = S4 in all).  The ONLY
    coordinate permutations among them are the identity (from +-I) and the
    light-cone reflection v0 <-> v2, which is Sym^2 of the swap [[0,1],[1,0]] -- determinant -1.
(2) So the forced pattern pi = (v0 <-> v2, v1 fixed) of all 108 perfect kick-tick-kick steps is the Sym^2 image of a
    det -1 element: the coset that the repository's determinant character marks as time-reversing (Pass 10955: the
    extended-Clifford unitary/antiunitary character, the clock CP grading, the D4 half-spin sheet swap).  A perfect
    interacting tick carries each past into the future of the coordinate that TIME REVERSAL of the light cone assigns
    to it, reversed in orientation.
(3) Four-clock potentials.  The four Hesse clocks are the four projective null directions n_1 = (1,0,0), n_2 = (0,0,1)
    (the light-cone axes, rays (1,0), (0,1)) and n_3 = (1,1,1), n_4 = (1,2,1) (the oblique rays (1,1), (2,1)).  Among the
    81 potentials N = sum_i c_i n_i n_i^T, exactly 36 make V(N) K V(N) perfect: c_1, c_2 are irrelevant (they are local
    phase gates on the light-cone qutrits, which sit at both ends of the step) and the condition is on the oblique
    clocks alone,  c_3 != c_4  and  c_3 + c_4 != 2  (from the determinants 1 + N_01 N_12 = -1 and (1 + N_02)^2 != 0 of
    Pass 11168, with N_01 = N_12 = c_3 - c_4, N_02 = c_3 + c_4).  The simplest: one oblique clock, N = (x0 + x1 + x2)^2.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11168_clock_dual_pair as D  # noqa: E402
import w33_pass11169_perfect_gate_local_orbits as O  # noqa: E402

OUT = ROOT / "data" / "w33_pass11174_lightcone_reflection_time_reversal.json"


def sym2(g):
    a, b = g[0]
    c, d = g[1]
    return np.array([[a * a, 2 * a * b, b * b], [a * c, a * d + b * c, b * d], [c * c, 2 * c * d, d * d]]) % 3


def summarize():
    GL = [np.array(m).reshape(2, 2) for m in itertools.product(range(3), repeat=4) if (m[0] * m[3] - m[1] * m[2]) % 3]
    images, dets, perms = {}, {}, []
    for g in GL:
        L = sym2(g)
        dg = int(round(np.linalg.det(g))) % 3
        dL = int(round(np.linalg.det(L))) % 3
        assert np.array_equal((L.T @ D.MQ @ L) % 3, D.MQ)
        images[L.tobytes()] = L
        dets[(dg, dL)] = dets.get((dg, dL), 0) + 1
        if ((L == 0) | (L == 1)).all() and (L.sum(0) == 1).all() and (L.sum(1) == 1).all():
            perms.append(dict(g=g.tolist(), det_g=dg, L=L.tolist()))
    swap_img = sym2(np.array([[0, 1], [1, 0]]))
    K = D.xz_to_party(D.I3, D.MQ, D.Z3, D.I3)
    pis = set()
    for n in itertools.product(range(3), repeat=6):
        Nm = np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])
        S = (D.kick(Nm) @ K @ D.kick(Nm)) % 3
        if D.perfect(S):
            pis.add(O.pattern(S))
    pi = list(pis)[0]
    pi_matrix = np.zeros((3, 3), int)
    for j in range(3):
        pi_matrix[pi[j], j] = 1
    rays = [(1, 0), (0, 1), (1, 1), (2, 1)]
    ns = [np.array([x * x, x * y, y * y]) % 3 for x, y in rays]
    perfect_c = []
    for c in itertools.product(range(3), repeat=4):
        N = sum(ci * np.outer(nv, nv) for ci, nv in zip(c, ns)) % 3
        if D.perfect((D.kick(N) @ K @ D.kick(N)) % 3):
            perfect_c.append(c)
    rule = all(((c[2] != c[3]) and ((c[2] + c[3]) % 3 != 2)) == (c in perfect_c) for c in itertools.product(range(3), repeat=4))
    res = dict(pass_id=11174, lorentz_images=len(images), det_pairs={str(k): v for k, v in dets.items()},
               coordinate_permutations=perms, lightcone_reflection_is_sym2_swap=bool(np.array_equal(swap_img, pi_matrix)),
               vkv_pi=list(pi), four_clock_perfect=len(perfect_c), four_clock_rule_exact=bool(rule),
               single_oblique_clock_perfect=[(0, 0, 1, 0) in perfect_c, (0, 0, 0, 1) in perfect_c],
               single_lightcone_clock_perfect=[(1, 0, 0, 0) in perfect_c, (0, 1, 0, 0) in perfect_c])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
