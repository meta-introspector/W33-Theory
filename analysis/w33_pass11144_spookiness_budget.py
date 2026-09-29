#!/usr/bin/env python3
"""Pass 11144: the spookiness budget of W(3,3) -- in space, correlations built from W(3,3)'s own (stabilizer) resources
never violate the qutrit Bell (CGLMP) inequality; a non-Clifford rotation is needed.  In time, one qutrit persisting
from a maximally mixed start reproduces the spatial Bell value exactly, and larger values reflect invasiveness, not
nonlocality.

CGLMP (d = 3, two settings): I = P(A0=B0) + P(B0=A1+1) + P(A1=B1) + P(B1=A0) - P(A0=B0-1) - P(B0=A1) - P(A1=B1-1)
- P(B1=A0-1); local hidden variables give I <= 2 (checked over all 81 deterministic strategies).
Spatial:  P(a,b|x,y) = <psi| Pi^x_a (x) Pi^y_b |psi>.
Temporal: one qutrit prepared in rho, measured projectively in basis x at t1 (Lueders), then in basis y at t2 (any
unitary in between is absorbed into y): P(a,b|x,y) = <a_x|rho|a_x> |<b_y|a_x>|^2.
Results:
  * validation: |Omega> with the CGLMP bases 2.872934 (Collins et al.); the Acin-Durt-Gisin-Latorre state 2.914854
    (the spatial quantum optimum, also reproduced by 6-start optimisation); local bound 2 over all 81 strategies;
  * SPACE, stabilizer only (every two-qutrit stabilizer state -- two local-Clifford orbits, represented by |Omega> and
    |00> -- times every choice of the four Pauli bases with every outcome labelling): max I = 2 exactly -- W(3,3)'s own
    resources give no violation (a non-negative discrete Wigner function on F_3^2 x F_3^2, the point space of W(3,3),
    is a local hidden-variable model; Gross);
  * the stabilizer point sits EXACTLY on the facet: along the path from the Pauli bases (t = 0) to the CGLMP bases
    (t = 1), I(t) = 2 + (pi^2/6) t^2 + O(t^3) (symmetric second difference 1.6449326 vs pi^2/6 = 1.6449341) -- any
    non-Clifford rotation along this path violates, with a quadratic onset;
  * TIME, maximally mixed start: I = 2.872934, the spatial |Omega> value -- one qutrit persisting in time reproduces the
    Bell statistics exactly (Pass 11143 on statistics);
  * TIME, stabilizer only (stabilizer starts, Pauli bases, Lueders updates): max I = 2 exactly -- still classical;
  * TIME, optimum over pure starts and projective bases (Lueders): I = 3.1628065 (40 restarts, reproducible from the
    frozen parameters; no closed form found, NOT sqrt10 = 3.16228), with a strongly biased start (outcome probabilities
    0.758, 0.121, 0.121).  For qubits the temporal CHSH optimum equals the spatial one (2 sqrt2, Fritz; reproduced
    here); for qutrit CGLMP the temporal optimum EXCEEDS the spatial one (3.163 > 2.915);
  * classical invasive models: without memory I <= 2; with any memory of >= 2 states I reaches the algebraic 4.
    So the orderings invert: in space classical 2 < quantum 2.915; in time quantum 3.163 < classical-with-memory 4 --
    measurement disturbance (non-orthogonal post-measurement states) caps the quantum signalling.
Reading: within W(3,3) (stabilizer states, Pauli measurements, Clifford gates) correlations are classical in space AND
in time.  "Spooky action at a distance" -- a no-signalling violation -- needs a resource outside W(3,3) (magic), and the
W(3,3) point is exactly on the edge, violated at second order.  The temporal mirror of the same correlations needs only
persistence; its larger values measure invasiveness, not nonlocality.  cf. qudit CHSH and Wigner negativity
(arXiv:2405.14367), Fritz (NJP 12, 083055), Budroni-Emary (arXiv:1309.3678).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11144_spookiness_budget.json"
W = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0)
Z = np.diag([1, W, W * W])


def cglmp(P):
    """P[x][y] is a 3x3 array P(a,b|x,y)"""
    def pr(x, y, rel):
        return sum(P[x][y][a, b] for a in range(3) for b in range(3) if rel(a, b))
    return (pr(0, 0, lambda a, b: a == b) + pr(1, 0, lambda a, b: b == (a + 1) % 3) + pr(1, 1, lambda a, b: a == b)
            + pr(0, 1, lambda a, b: b == a) - pr(0, 0, lambda a, b: a == (b - 1) % 3) - pr(1, 0, lambda a, b: b == a)
            - pr(1, 1, lambda a, b: a == (b - 1) % 3) - pr(0, 1, lambda a, b: b == (a - 1) % 3))


def lhv_max():
    best = -9
    for a0, a1, b0, b1 in itertools.product(range(3), repeat=4):
        P = [[np.zeros((3, 3)) for _ in range(2)] for _ in range(2)]
        for x, a in ((0, a0), (1, a1)):
            for y, b in ((0, b0), (1, b1)):
                P[x][y][a, b] = 1
        best = max(best, cglmp(P))
    return best


def spatial(psi, A, B):
    P = [[np.zeros((3, 3)) for _ in range(2)] for _ in range(2)]
    for x in range(2):
        for y in range(2):
            for a in range(3):
                for b in range(3):
                    v = np.kron(A[x][:, a], B[y][:, b]).conj() @ psi
                    P[x][y][a, b] = abs(v) ** 2
    return cglmp(P)


def temporal(phi, A, B):
    P = [[np.zeros((3, 3)) for _ in range(2)] for _ in range(2)]
    for x in range(2):
        for y in range(2):
            for a in range(3):
                pa = abs(A[x][:, a].conj() @ phi) ** 2
                for b in range(3):
                    P[x][y][a, b] = pa * abs(B[y][:, b].conj() @ A[x][:, a]) ** 2
    return cglmp(P)


def temporal_mixed(A, B):
    P = [[np.zeros((3, 3)) for _ in range(2)] for _ in range(2)]
    for x in range(2):
        for y in range(2):
            for a in range(3):
                for b in range(3):
                    P[x][y][a, b] = abs(B[y][:, b].conj() @ A[x][:, a]) ** 2 / 3
    return cglmp(P)


def cglmp_bases(t=1.0):
    """Collins et al. bases |k>_x = sum_j w^(j(k + alpha_x)) |j>/sqrt3 (Alice), |l>_y = sum_j w^(j(-l + beta_y)) |j>/sqrt3
    (Bob), alpha = (0, t/2), beta = (t/4, -t/4); t = 1 is optimal, t = 0 gives Pauli (Fourier) bases"""
    def basis(al, sg):
        return np.array([[np.exp(2j * np.pi * j * (sg * k + al) / 3) for k in range(3)] for j in range(3)]) / np.sqrt(3)
    A = [basis(0, 1), basis(t / 2, 1)]
    B = [basis(t / 4, -1), basis(-t / 4, -1)]
    return A, B


def omega():
    return np.array([1 if i == j else 0 for i in range(3) for j in range(3)], complex) / np.sqrt(3)


def unitary(p):
    H = np.zeros((3, 3), complex)
    iu = np.triu_indices(3, 1)
    H[iu] = p[:3] + 1j * p[3:6]
    H = H + H.conj().T + np.diag(p[6:9])
    return expm(1j * H)


def opt(kind, starts=40, seed=0):
    rng = np.random.default_rng(seed)
    def f(p):
        A = [unitary(p[0:9]), unitary(p[9:18])]
        B = [unitary(p[18:27]), unitary(p[27:36])]
        if kind == 'spatial':
            v = p[36:45] + 1j * p[45:54]
            return -spatial(v / np.linalg.norm(v), A, B)
        v = p[36:39] + 1j * p[39:42]
        return -temporal(v / np.linalg.norm(v), A, B)
    n = 54 if kind == 'spatial' else 42
    best = max((-minimize(f, rng.normal(size=n), method='L-BFGS-B', options=dict(maxiter=3000)).fun for _ in range(starts)))
    return float(best)


def pauli_bases():
    out = []
    for op in (Z, X, X @ Z, X @ Z @ Z):
        _, V = np.linalg.eig(op)
        V = V / np.linalg.norm(V, axis=0)
        for perm in itertools.permutations(range(3)):
            out.append(V[:, perm])
    return out


def stabilizer_states_two():
    """all 360 two-qutrit stabilizer states: joint eigenvectors of the 40 Lagrangian groups (form [a,a'] + [b,b'])"""
    import w33_pass11143_space_time_two_quadrangles as Q
    states = []
    for line in Q.lagrangians(Q.form(1)):
        states += Q.stab_states(line)
    return states


def stabilizer_states_one():
    out = []
    for op in (Z, X, X @ Z, X @ Z @ Z):
        _, V = np.linalg.eig(op)
        out += [V[:, i] / np.linalg.norm(V[:, i]) for i in range(3)]
    return out


def _terms(PR):
    """PR[i][j]: 3x3 P(a,b) for Alice basis i, Bob basis j -> the four CGLMP term matrices"""
    n = len(PR)
    T = {k: np.zeros((n, n)) for k in ('00', '10', '11', '01')}
    for i in range(n):
        for j in range(n):
            M = PR[i][j]
            eq = sum(M[a, a] for a in range(3)); up = sum(M[a, (a + 1) % 3] for a in range(3))
            dn = sum(M[a, (a - 1) % 3] for a in range(3))
            T['00'][i, j] = eq - up          # P(A0=B0) - P(A0=B0-1)
            T['10'][i, j] = up - eq          # P(B0=A1+1) - P(B0=A1)
            T['11'][i, j] = eq - up          # P(A1=B1) - P(A1=B1-1)
            T['01'][i, j] = eq - dn          # P(B1=A0) - P(B1=A0-1)
    return T


def _max_over_settings(T):
    I = (T['00'][:, None, :, None] + T['10'][None, :, :, None] + T['11'][None, :, None, :] + T['01'][:, None, None, :])
    return float(I.max())


def spatial_stabilizer_max():
    """local Cliffords act transitively on the 216 entangled and on the 144 product stabilizer states and preserve the
    set of Pauli bases, so one representative of each suffices"""
    PB = pauli_bases()
    reps = [omega(), np.eye(9)[0].astype(complex)]
    best = -9
    for psi in reps:
        PR = [[np.array([[abs(np.kron(A[:, a], B[:, b]).conj() @ psi) ** 2 for b in range(3)] for a in range(3)])
               for B in PB] for A in PB]
        best = max(best, _max_over_settings(_terms(PR)))
    return best


def temporal_stabilizer_max():
    PB = pauli_bases()
    best = -9
    for phi in stabilizer_states_one() + [None]:
        PR = [[np.array([[(abs(A[:, a].conj() @ phi) ** 2 if phi is not None else 1 / 3) * abs(B[:, b].conj() @ A[:, a]) ** 2
                          for b in range(3)] for a in range(3)]) for B in PB] for A in PB]
        best = max(best, _max_over_settings(_terms(PR)))
    return best


def magic_threshold():
    ts = np.linspace(0, 1, 201)
    vals = [spatial(omega(), *cglmp_bases(t)) for t in ts]
    t_star = next((float(t) for t, v in zip(ts, vals) if v > 2 + 1e-9), None)
    return t_star, float(max(vals))


def classical_invasive_max(m):
    """classical system with an m-state memory, invasive measurements: Alice's input x sets outcome a_x and memory s_x;
    Bob's outcome is any function of (y, memory)"""
    best = -9
    for a0, a1 in itertools.product(range(3), repeat=2):
        for s0, s1 in itertools.product(range(m), repeat=2):
            for b in itertools.product(range(3), repeat=4):
                P = [[np.zeros((3, 3)) for _ in range(2)] for _ in range(2)]
                for x, (a, s) in enumerate(((a0, s0), (a1, s1))):
                    for y in range(2):
                        P[x][y][a, b[2 * y + (0 if s == s0 else 1)]] = 1
                best = max(best, cglmp(P))
    return best


def temporal_certificate():
    d = json.loads((ROOT / "data" / "w33_pass11144_temporal_optimum_params.json").read_text())
    p = np.array(d['params'])
    A = [unitary(p[0:9]), unitary(p[9:18])]
    B = [unitary(p[18:27]), unitary(p[27:36])]
    v = p[36:39] + 1j * p[39:42]
    v = v / np.linalg.norm(v)
    return temporal(v, A, B), [float(abs(A[0][:, a].conj() @ v) ** 2) for a in range(3)]


def summarize():
    A, B = cglmp_bases(1.0)
    g = (np.sqrt(11) - np.sqrt(3)) / 2
    adgl = np.zeros(9, complex); adgl[0] = 1; adgl[4] = g; adgl[8] = 1
    tval, tprobs = temporal_certificate()
    f = lambda t: spatial(omega(), *cglmp_bases(t))
    h = 1e-3
    curv = (f(h) + f(-h) - 2 * f(0.0)) / (2 * h * h)
    res = dict(pass_id=11144, lhv_max=lhv_max(),
               spatial_omega_cglmp_bases=spatial(omega(), A, B),
               spatial_optimum_adgl=spatial(adgl / np.linalg.norm(adgl), A, B),
               temporal_mixed_equals_spatial=temporal_mixed([a.conj() for a in A], B),
               temporal_optimum=tval, temporal_optimum_first_outcome_probs=tprobs,
               spatial_stabilizer_max=spatial_stabilizer_max(), temporal_stabilizer_max=temporal_stabilizer_max(),
               classical_invasive_max_1state=classical_invasive_max(1), classical_invasive_max_2states=classical_invasive_max(2),
               pauli_point_value=spatial(omega(), *cglmp_bases(0.0)),
               magic_path_curvature=curv, pi2_over_6=np.pi ** 2 / 6)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(ROOT / "analysis"))
    print(json.dumps(summarize(), indent=1))
