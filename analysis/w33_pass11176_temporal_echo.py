#!/usr/bin/env python3
"""Pass 11176: entanglement across time through a scrambling tick is ALL-OR-NOTHING; a perfect tick erases a qutrit's
entanglement with its past until an echo at a later tick; and a temporal Bell violation survives where temporal
entanglement does not.

Setup.  Three qutrits (the 27 events of the paper's clock).  One qutrit is the system, the other two its environment,
maximally mixed at the start.  The system's n-tick channel is N_n(rho) = Tr_env[U^n (rho (x) I/9) U^-n]; its two-time
pseudo-density operator (Pass 11148) is R_n = (1/3)(id (x) N_n)(SWAP), and the temporal negativity N(n) is the sum of
|negative eigenvalues| of R_n (identity: 1).  The register stays maximally mixed, so the three-time polygamy sum of
Pass 11148 is N(t1,t2) + N(t2,t3) + N(t1,t3) = 2 N(1) + N(2).
THEOREM (all-or-nothing).  For a Clifford tick with symplectic matrix S, N(1) = 1 for qutrit q if its column of blocks is
local (S_iq = 0 for i != q) and N(1) = 0 otherwise.  Proof: only the Paulis of q whose images stay on q survive the
partial trace; they form a subgroup of F_3^2 of order 1, 3 or 9, and with cube-root-of-unity phases the eigenvalues of
R are proportional to 1 + 2 cos(2 pi m / 3) >= 0 unless the whole group survives.  Checked on 120 (qutrit, tick) pairs.
RESULTS.
  * identity: 1, 1, 1 (three-time sum 3);  the clock K (order 3): transverse qutrit 1, light-cone qutrits 0, full revival
    at the third tick (three-time sums 0, 3, 0);
  * the perfect tick V(N) K V(N) with N = (x0 + x1 + x2)^2 (order 9): 0 for every qutrit for eight ticks, then 1 -- total
    temporal amnesia with an echo at the ninth tick; random perfect ticks: 0 after one and two ticks;
  * echo census over all 108 perfect kick-tick-kick ticks: orders 7 to 36; the transverse qutrit's entanglement with its
    past often returns first (e.g. at tick 3 of 12, 4 of 24, 9 of 36) -- recorded in vkv_echo_census;
  * BELL WITHOUT TEMPORAL ENTANGLEMENT: the two-time CGLMP value (d = 3) through one clock tick is 3.1628 for the
    transverse qutrit (the identity value) and 3.0606 > 2 for the light-cone qutrits, whose temporal negativity is 0;
    through a perfect tick it is exactly 0 (the channel is completely depolarising).  In time, a Bell violation does not
    witness entanglement: it witnesses the setting-dependent memory the first measurement leaves behind.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11166_clock_tick_clifford as C  # noqa: E402
import w33_pass11168_clock_dual_pair as D  # noqa: E402

OUT = ROOT / "data" / "w33_pass11176_temporal_echo.json"
W = np.exp(2j * np.pi / 3)
X = np.array(list(itertools.product(range(3), repeat=3)))            # index = 9 x0 + 3 x1 + x2 (C.V order)


def kick_unitary(N):
    ph = np.einsum('ki,ij,kj->k', X, N % 3, X) % 3
    return np.diag(W ** ((2 * ph) % 3))


def fourier_on(q):
    F = np.array([[W ** (a * b) for b in range(3)] for a in range(3)]) / np.sqrt(3)
    ops = [np.eye(3)] * 3
    ops[q] = F
    return np.kron(np.kron(ops[0], ops[1]), ops[2])


def reduced_pdo(U, q):
    """R = (1/3) sum_ij |i><j| (x) N(|j><i|), N = reduced channel of U on qutrit q (others maximally mixed)"""
    R = np.zeros((9, 9), complex)
    for i in range(3):
        for j in range(3):
            E = np.zeros((3, 3))
            E[j, i] = 1
            ops = [np.eye(3) / 3] * 3
            ops[q] = E
            X_in = np.kron(np.kron(ops[0], ops[1]), ops[2])
            Y = (U @ X_in @ U.conj().T).reshape([3] * 6)
            # trace out the two environment qutrits
            env = [k for k in range(3) if k != q]
            Y = np.trace(Y, axis1=env[0], axis2=env[0] + 3)
            # after the first trace the remaining axes shift: recompute positions
            env2 = env[1] - (1 if env[1] > env[0] else 0)
            Y = np.trace(Y, axis1=env2, axis2=env2 + 2)
            Ei = np.zeros((3, 3))
            Ei[i, j] = 1
            R += np.kron(Ei, Y) / 3
    return R


def tneg(U, q):
    R = reduced_pdo(U, q)
    ev = np.linalg.eigvalsh((R + R.conj().T) / 2)
    return float(-ev[ev < 0].sum())


def order(U, maxn=400):
    """smallest n with U^n a scalar multiple of the identity"""
    P = np.eye(27, dtype=complex)
    for n in range(1, maxn + 1):
        P = P @ U
        if np.allclose(P, P[0, 0] * np.eye(27), atol=1e-8):
            return n
    return None


def echo(U, nmax):
    out = []
    P = np.eye(27, dtype=complex)
    for n in range(1, nmax + 1):
        P = P @ U
        out.append([round(tneg(P, q), 6) for q in range(3)])
    return out


def column_local(S, q):
    return all(not S[2 * i:2 * i + 2, 2 * q:2 * q + 2].any() for i in range(3) if i != q)


def echo_times(S, maxn=100):
    """first n >= 1 at which each qutrit's column of S^n is local (temporal negativity returns to 1), and the order of S"""
    P = np.eye(6, dtype=np.int64)
    first = [None] * 3
    for n in range(1, maxn + 1):
        P = (P @ S) % 3
        for q in range(3):
            if first[q] is None and column_local(P, q):
                first[q] = n
        if np.array_equal(P, np.eye(6, dtype=np.int64)):
            return first, n
    return first, None


def channel_kraus(U, q):
    """Kraus operators of the reduced channel on qutrit q (environment maximally mixed)"""
    Us = U.reshape([3] * 6)
    env = [k for k in range(3) if k != q]
    ks = []
    for e_out in itertools.product(range(3), repeat=2):
        for e_in in itertools.product(range(3), repeat=2):
            idx_out, idx_in = [slice(None)] * 3, [slice(None)] * 3
            for k, eo, ei in zip(env, e_out, e_in):
                idx_out[k], idx_in[k] = eo, ei
            ks.append(Us[tuple(idx_out) + tuple(idx_in)] / 3)
    return ks


def temporal_cglmp_through(ks, restarts=12, seed=0):
    """max over state and bases of the d = 3 CGLMP value with A_x, then the channel, then B_y"""
    from scipy.optimize import minimize
    C = np.zeros((2, 2, 3, 3))
    for a in range(3):
        for b_ in range(3):
            C[0, 0, a, b_] += (a == b_) - (a == (b_ - 1) % 3)
            C[1, 0, a, b_] += (b_ == (a + 1) % 3) - (b_ == a)
            C[1, 1, a, b_] += (a == b_) - (a == (b_ - 1) % 3)
            C[0, 1, a, b_] += (b_ == a) - (b_ == (a - 1) % 3)

    def unit(p):
        H = np.zeros((3, 3), complex)
        iu = np.triu_indices(3, 1)
        H[iu] = p[:3] + 1j * p[3:6]
        H = H + H.conj().T
        H[np.diag_indices(3)] = p[6:9]
        w, V = np.linalg.eigh(H)
        return V @ np.diag(np.exp(1j * w)) @ V.conj().T

    def val(p):
        A = [unit(p[0:9]), unit(p[9:18])]
        B = [unit(p[18:27]), unit(p[27:36])]
        v = p[36:39] + 1j * p[39:42]
        v = v / np.linalg.norm(v)
        tot = 0.0
        for x in range(2):
            for a in range(3):
                u = A[x][:, a]
                pa = abs(u.conj() @ v) ** 2
                rho = sum(k @ np.outer(u, u.conj()) @ k.conj().T for k in ks)
                for y in range(2):
                    for b_ in range(3):
                        w_ = B[y][:, b_]
                        tot += C[x, y, a, b_] * pa * float(np.real(w_.conj() @ rho @ w_))
        return -tot
    rng = np.random.default_rng(seed)
    best = -9
    for _ in range(restarts):
        r = minimize(val, rng.normal(size=42), method='BFGS', options=dict(maxiter=3000))
        best = max(best, -r.fun)
    return best


def summarize():
    rng = np.random.default_rng(11176)
    K = C.clock()
    Nclock = np.ones((3, 3), int)                            # (x0 + x1 + x2)^2
    Vn = kick_unitary(Nclock)
    T = Vn @ K @ Vn
    S_T = C.symplectic(T)
    sym_check = bool(np.array_equal(S_T, (D.kick(Nclock) @ D.xz_to_party(D.I3, D.MQ, D.Z3, D.I3) @ D.kick(Nclock)) % 3))
    ident = np.eye(27)
    res = dict(pass_id=11176,
               identity_N1=[round(tneg(ident, q), 6) for q in range(3)],
               clock_order=order(K), clock_echo=echo(K, 3),
               perfect_tick_is_perfect=D.perfect(S_T), perfect_tick_symplectic_matches=sym_check,
               perfect_tick_order=order(T), perfect_tick_echo=echo(T, order(T) or 12))
    # random Clifford ticks: products of kicks and local Fourier transforms
    perf_N1, ctrl_N1, perf_echo2 = [], [], []
    tries = 0
    while (len(perf_N1) < 8 or len(ctrl_N1) < 8) and tries < 400:
        tries += 1
        U = np.eye(27, dtype=complex)
        for _ in range(12):
            n = rng.integers(0, 3, 6)
            Nm = np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])
            U = fourier_on(int(rng.integers(0, 3))) @ kick_unitary(Nm) @ U
        S = C.symplectic(U)
        n1 = [round(tneg(U, q), 6) for q in range(3)]
        if D.perfect(S) and len(perf_N1) < 8:
            perf_N1.append(n1)
            perf_echo2.append([round(tneg(U @ U, q), 6) for q in range(3)])
        elif not D.perfect(S) and len(ctrl_N1) < 8:
            ctrl_N1.append(n1)
    res.update(random_perfect_N1=perf_N1, random_perfect_N2=perf_echo2, random_nonperfect_N1=ctrl_N1,
               perfect_N1_all_zero=all(abs(x) < 1e-9 for row in perf_N1 for x in row) and
               all(abs(x) < 1e-9 for x in res['perfect_tick_echo'][0]))
    # all-or-nothing theorem: temporal negativity 1 iff the qutrit's column is local, else 0 (census of random ticks)
    agree, n_checked = True, 0
    for _ in range(40):
        U = np.eye(27, dtype=complex)
        for _ in range(int(rng.integers(1, 6))):
            n = rng.integers(0, 3, 6)
            Nm = np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])
            U = fourier_on(int(rng.integers(0, 3))) @ kick_unitary(Nm) @ U
        S = C.symplectic(U)
        for q in range(3):
            n_checked += 1
            agree &= abs(tneg(U, q) - (1.0 if column_local(S, q) else 0.0)) < 1e-9
    res.update(all_or_nothing_agrees=bool(agree), all_or_nothing_checked=n_checked)
    # echo times for all 108 perfect kick-tick-kick ticks (symplectic level)
    Kx = D.xz_to_party(D.I3, D.MQ, D.Z3, D.I3)
    from collections import Counter
    ech = Counter()
    for n in itertools.product(range(3), repeat=6):
        Nm = np.array([[n[0], n[1], n[2]], [n[1], n[3], n[4]], [n[2], n[4], n[5]]])
        S = (D.kick(Nm) @ Kx @ D.kick(Nm)) % 3
        if D.perfect(S):
            first, o = echo_times(S)
            ech[(tuple(first), o)] += 1
    res['vkv_echo_census'] = {str(k): v for k, v in sorted(ech.items(), key=lambda kv: str(kv[0]))}
    # two-time Bell value (d = 3 CGLMP) through one clock tick, per qutrit, and through the perfect tick
    res['cglmp_through_clock'] = [round(temporal_cglmp_through(channel_kraus(K, q)), 5) for q in range(3)]
    res['cglmp_through_perfect_tick'] = round(temporal_cglmp_through(channel_kraus(T, 0), restarts=4), 5)
    res['three_time_sum'] = dict(identity=3.0,
                                 clock=[2 * a + b for a, b in zip(res['clock_echo'][0], res['clock_echo'][1])],
                                 perfect_tick=[2 * a + b for a, b in zip(res['perfect_tick_echo'][0], res['perfect_tick_echo'][1])])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
