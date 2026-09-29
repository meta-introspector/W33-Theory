#!/usr/bin/env python3
"""Pass 11153: in time, the magic must sit in the measurements -- the start state cannot carry it (the opposite of space,
where either can).

Temporal CGLMP (one qutrit, Lueders projective measurements at two times, Pass 11144):
  P(a,b|x,y) = <a_x|rho|a_x> |<b_y|a_x>|^2.
  (a) Pauli measurements only (the four Pauli eigenbases, all labellings; Alice's first fixed by Clifford symmetry), ANY
      start: for fixed bases I is linear in rho, so the maximum over starts is the top eigenvalue of a 3x3 operator --
      enumerated exactly;
  (b) ANY bases, stabilizer start: with arbitrary bases every pure start is unitarily equivalent to every other (rotate
      all four bases), so the optimum equals the unrestricted temporal optimum 3.1628065 (Pass 11144) -- the start's
      magic is irrelevant.
Compare space (Passes 11144, 11149): Pauli measurements + magic state 2.6017 (> 2); stabilizer state + magic
measurements 2.8729.  Result for (a): see summarize().
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11153_temporal_magic_budget.json"
W3 = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0)
Z = np.diag([1, W3, W3 * W3])


def pauli_bases():
    out = []
    for op in (Z, X, X @ Z, X @ Z @ Z):
        _, V = np.linalg.eig(op)
        V = V / np.linalg.norm(V, axis=0)
        for perm in itertools.permutations(range(3)):
            out.append(V[:, perm])
    return out


REL = [((0, 0), lambda a, b: a == b, lambda a, b: a == (b - 1) % 3), ((1, 0), lambda a, b: b == (a + 1) % 3, lambda a, b: b == a),
       ((1, 1), lambda a, b: a == b, lambda a, b: a == (b - 1) % 3), ((0, 1), lambda a, b: b == a, lambda a, b: b == (a - 1) % 3)]


def temporal_operator(A, B):
    """3x3 operator O with I = <psi|O|psi> for the temporal CGLMP with bases A[x], B[y]"""
    O = np.zeros((3, 3), complex)
    for (x, y), plus, minus in REL:
        for a in range(3):
            for b in range(3):
                s = (1 if plus(a, b) else 0) - (1 if minus(a, b) else 0)
                if s:
                    O += s * abs(B[y][:, b].conj() @ A[x][:, a]) ** 2 * np.outer(A[x][:, a], A[x][:, a].conj())
    return O


def summarize():
    PB = pauli_bases()
    zb = PB[0]
    best = -9
    for a1, b0, b1 in itertools.product(range(len(PB)), repeat=3):
        O = temporal_operator([zb, PB[a1]], [PB[b0], PB[b1]])
        best = max(best, float(np.linalg.eigvalsh((O + O.conj().T) / 2)[-1]))
    res = dict(pass_id=11153, temporal_pauli_measurements_any_start=best,
               temporal_any_bases_stabilizer_start=3.1628065077657714,
               spatial_pauli_measurements_magic_state=2.6016791318831545, spatial_stabilizer_state_magic_measurements=2.872934051172338)
    res['in_time_magic_must_be_in_measurements'] = best <= 2 + 1e-9
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
