#!/usr/bin/env python3
"""Pass 11149: how much magic does spooky action need?  Where the magic sits -- in the state or in the measurements --
and how much of it (mana) the Bell-violating resources carry.

Setting: the qutrit CGLMP inequality (local bound 2; Pass 11144).  Magic is measured by the mana
M(rho) = log sum_u |W_rho(u)| of the discrete Wigner function on F_3^(2n) (phase-point operators A_u = D_u A_0 D_u^dag,
A_0 the parity; Veitch-Mousavian-Gottesman-Emerson).  Stabilizer states have M = 0.
Computed:
  (a) Pauli measurements only (four Pauli eigenbases, all outcome labellings; Alice's first basis fixed by local-Clifford
      symmetry), optimal state = top eigenvector of each Bell operator: the maximum over ALL two-qutrit states;
  (b) the mana of the optimal state of (a), of the Acin-Durt-Gisin-Latorre state, and of |Omega>;
  (c) the mana of the CGLMP measurement vectors (magic carried by the measurements when the state is |Omega>).
  (d) the minimum mana of a state that reaches 2 + eps with the optimal Pauli-measurement Bell operator
      (analysis/w33_pass11149_scan_min_mana.py; frozen data/w33_pass11149_min_mana_curve.json).
Results:
  * with Pauli measurements ONLY, the best state reaches 2 sqrt2 cos((1/3) arccos(1/(2 sqrt2))) = 2.6016791, the largest
    root of x^3 - 6x - 2 = 0 (the Bell operator's spectrum splits into the roots of x^3 - 6x - 2 and of x^3 - 3x + 1,
    i.e. 2cos(pi/9), 2cos(2pi/9), 2cos(4pi/9)); that state carries mana 0.392 -- magic in the STATE alone suffices;
  * the stabilizer state |Omega> (mana 0) with magic measurements (CGLMP vectors of mana up to 0.437) reaches 2.8729 --
    magic in the MEASUREMENTS alone suffices; the optimal Acin et al. state carries mana 0.173;
  * NO threshold: the minimum mana needed for a violation 2 + eps rises linearly from zero, about 0.63-0.66 mana per unit
    of violation (eps = 0.02, 0.1, 0.3, 0.6 -> mana 0.0128, 0.064, 0.183, 0.378) -- infinitesimal magic buys a small
    spooky violation, the state-side counterpart of Pass 11144's quadratic onset on the measurement side.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11149_magic_for_spookiness.json"
W3 = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0)
Z = np.diag([1, W3, W3 * W3])


def D(u):
    # Heisenberg-Weyl displacement with the symmetric phase for odd d: w^(2^-1 x z) X^x Z^z, 2^-1 = 2 mod 3
    return W3 ** ((2 * u[0] * u[1]) % 3) * np.linalg.matrix_power(X, u[0] % 3) @ np.linalg.matrix_power(Z, u[1] % 3)


A0 = sum(D(u) for u in itertools.product(range(3), repeat=2)) / 3
PTS = list(itertools.product(range(3), repeat=2))
A1 = {u: D(u) @ A0 @ D(u).conj().T for u in PTS}


def wigner1(rho):
    return np.array([np.trace(rho @ A1[u]).real / 3 for u in PTS])


def wigner2(rho):
    return np.array([np.trace(rho @ np.kron(A1[u], A1[v])).real / 9 for u in PTS for v in PTS])


def mana(rho, n):
    w = wigner1(rho) if n == 1 else wigner2(rho)
    return float(np.log(np.abs(w).sum()))


def pauli_bases():
    out = []
    for op in (Z, X, X @ Z, X @ Z @ Z):
        _, V = np.linalg.eig(op)
        V = V / np.linalg.norm(V, axis=0)
        for perm in itertools.permutations(range(3)):
            out.append(V[:, perm])
    return out


def bell_operator(A0_, A1_, B0_, B1_):
    """CGLMP Bell operator (d = 3) for bases given as columns"""
    P = lambda M, a: np.outer(M[:, a], M[:, a].conj())
    def term(Am, Bm, rel_plus, rel_minus):
        T = np.zeros((9, 9), complex)
        for a in range(3):
            for b in range(3):
                s = (1 if rel_plus(a, b) else 0) - (1 if rel_minus(a, b) else 0)
                if s:
                    T += s * np.kron(P(Am, a), P(Bm, b))
        return T
    return (term(A0_, B0_, lambda a, b: a == b, lambda a, b: a == (b - 1) % 3)
            + term(A1_, B0_, lambda a, b: b == (a + 1) % 3, lambda a, b: b == a)
            + term(A1_, B1_, lambda a, b: a == b, lambda a, b: a == (b - 1) % 3)
            + term(A0_, B1_, lambda a, b: b == a, lambda a, b: b == (a - 1) % 3))


def pauli_only_max():
    PB = pauli_bases()
    zb = PB[0]
    best, arg = -9, None
    for a1, b0, b1 in itertools.product(range(len(PB)), repeat=3):
        Bop = bell_operator(zb, PB[a1], PB[b0], PB[b1])
        ev, V = np.linalg.eigh((Bop + Bop.conj().T) / 2)
        if ev[-1] > best + 1e-12:
            best, arg = float(ev[-1]), V[:, -1]
    return best, arg


def cglmp_bases(t=1.0):
    def basis(al, sg):
        return np.array([[np.exp(2j * np.pi * j * (sg * k + al) / 3) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    return [basis(0, 1), basis(t / 2, 1)], [basis(t / 4, -1), basis(-t / 4, -1)]


def summarize():
    best, v = pauli_only_max()
    rho = np.outer(v, v.conj())
    g = (np.sqrt(11) - np.sqrt(3)) / 2
    adgl = np.zeros(9, complex); adgl[[0, 4, 8]] = [1, g, 1]; adgl /= np.linalg.norm(adgl)
    om = np.zeros(9, complex); om[[0, 4, 8]] = 1; om /= np.sqrt(3)
    A, B = cglmp_bases(1.0)
    meas_mana = sorted({round(mana(np.outer(M[:, k], M[:, k].conj()), 1), 9) for M in A + B for k in range(3)})
    A0_, B0_ = cglmp_bases(0.0)
    cubic_root = float(max(np.roots([1, 0, -6, -2]).real))
    trig = float(2 * np.sqrt(2) * np.cos(np.arccos(1 / (2 * np.sqrt(2))) / 3))
    res = dict(pass_id=11149,
               pauli_measurements_best_state_value=best, cubic_x3_6x_2_root=cubic_root, trig_form=trig,
               pauli_best_state_mana=mana(rho, 2),
               omega_mana=mana(np.outer(om, om.conj()), 2), adgl_mana=mana(np.outer(adgl, adgl.conj()), 2),
               omega_with_cglmp_bases_value=float(np.real(om.conj() @ bell_operator(A[0], A[1], B[0], B[1]) @ om)),
               cglmp_measurement_vector_mana=meas_mana,
               min_mana_curve=json.loads((ROOT / "data" / "w33_pass11149_min_mana_curve.json").read_text()),
               pauli_basis_vector_mana=sorted({round(mana(np.outer(M[:, k], M[:, k].conj()), 1), 9) for M in A0_ + B0_ for k in range(3)}))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
