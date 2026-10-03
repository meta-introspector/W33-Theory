"""Pass 11350: one magic gate is F3 linear algebra -- an exact linear decider for C (T (x) I...), any number of qutrits.

From Pass 11331: U = C T1 (C = W(a) V_M, T1 = T on qutrit 1) is substrate-time-reversible iff some Clifford
E = W(e) V_Q solves  T1 E T1^-1 = lambda C^-1 E C^T.  Both sides are written in the normal form W(frame) V_(symplectic):
  * RHS (verified numerically):  C^-1 E C^T  ~  W(M^-1 (e - a) - P J a) V_P,   P = M^-1 Q J M^-1 J,  J = diag(1,-1,...);
  * LHS (verified numerically):  T1 W(e) V_Q T1^-1  ~  W(g(e) + s^-k r_Q) V_{s^-k Q},  k = e_x1, s the unit shear
    (x1, z1) -> (x1, z1 + x1); r_Q is the frame of T1 V_Q T1^-1 V_Q^dag and g(e) the frame of T1 W(e) T1^-1 V_{s^-k}^dag.
Matching symplectic parts and frames (phases are free):
  (S)  M s^-k Q (J M J) = Q,  Q z1 = z1,  Q symplectic          -- linear in Q plus a symplecticity filter;
  (F)  g(e) + s^-k r_Q = M^-1 (e - a) - s^-k Q J a,  e_x1 = k   -- affine over F3 in (e, a) for fixed (Q, k).
So the reversible frames of a class are a UNION of affine subspaces (one per solution of (S)); the census (Pass 11330)
shows the union is always a single affine subspace.  The decider here enumerates the solutions of (S) (an affine space
over F3 intersected with Sp) and solves (F) by Gaussian elimination; it falls back to 'undecided' if (S) has more than
CAP candidate solutions.  Validation: all 51,840 two-qutrit classes against the exhaustive frame counts of Pass 11330.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11350_linear_decider.json"
CAP = 3 ** 9


# ---------------------------------------------------------------- F3 linear algebra
def rref3(A):
    A = A.copy() % 3
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        if A[r, c] == 2:
            A[r] = (2 * A[r]) % 3
        nz = np.flatnonzero(A[:, c])
        nz = nz[nz != r]
        A[nz] = (A[nz] - np.outer(A[nz, c], A[r])) % 3
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A[:r], piv


def solve_affine(A, b):
    """all x with A x = b over F3: returns (x0, basis of nullspace) or None"""
    n = A.shape[1]
    Ab = np.concatenate([A % 3, (b % 3)[:, None]], axis=1)
    Rr, piv = rref3(Ab)
    if n in piv:
        return None
    x0 = np.zeros(n, dtype=np.int64)
    for i, c in enumerate(piv):
        x0[c] = Rr[i, n]
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-Rr[i, f]) % 3
        basis.append(v)
    return x0, basis


class Decider:
    def __init__(self, n, rng=None):
        self.n = n
        self.wl = R.Weyl(n)
        R.WEYL[n] = self.wl
        self.rng = rng or np.random.default_rng(0)
        self.N2 = 2 * n
        self.T1 = R.local(n, 0, P2.T)
        self.J = np.diag([1, 2] * n)
        self.z1 = np.zeros(self.N2, dtype=np.int64)
        self.z1[1] = 1
        s = np.eye(self.N2, dtype=np.int64)
        s[1, 0] = 1
        self.s_pow = [np.linalg.matrix_power(s, k) % 3 for k in range(3)]
        self.sinv_pow = [self.s_pow[(-k) % 3] for k in range(3)]
        self.g = {}
        self._weil_cache = {}

    def weil(self, M):
        """canonical Weil unitary V_M (V W(p) V^dag = W(Mp)), built as sum_p W(Mp)|j><0|W(p)^dag (a single matrix
        product), checked on the generators"""
        key = tuple(M.ravel())
        if key in self._weil_cache:
            return self._weil_cache[key]
        wl = self.wl
        if not hasattr(self, "_col0"):
            self._col0 = wl.W[:, :, 0]                                   # W(p)|0>, shape (P, D)
        img = ((wl.labels @ M.T) % 3) @ wl.pow3
        V = None
        for j in range(wl.D):
            cols = wl.W[img][:, :, j]                                     # W(Mp)|j>
            V0 = cols.T @ self._col0.conj()                               # sum_p W(Mp)|j><0|W(p)^dag
            nrm = np.linalg.norm(V0)
            if nrm > 1e-6:
                V = V0 / nrm * np.sqrt(wl.D)
                break
        for i in range(2 * self.n):
            e = np.zeros(2 * self.n, dtype=np.int64)
            e[i] = 1
            assert np.allclose(V @ wl.W[wl.index(e)] @ V.conj().T, wl.W[wl.index(M @ e % 3)], atol=1e-8)
        self._weil_cache[key] = V
        return V

    def frame_of(self, U):
        ov = np.abs(np.einsum('pij,ij->p', self.wl.W.conj(), U)) / self.wl.D
        k = int(np.argmax(ov))
        assert np.isclose(ov[k], 1), "not a Weyl operator"
        return self.wl.labels[k].astype(np.int64)

    def gframe(self, e):
        key = tuple(e)
        if key not in self.g:
            k = int(e[0])
            X = self.T1 @ self.wl.W[self.wl.index(e)] @ self.T1.conj().T @ self.weil(self.sinv_pow[k]).conj().T
            self.g[key] = self.frame_of(X)
        return self.g[key]

    def r_frame(self, Q):
        V = self.weil(Q)
        return self.frame_of(self.T1 @ V @ self.T1.conj().T @ V.conj().T)

    def symplectic_solutions(self, M, k):
        """affine solution space of M s^-k Q (J M J) = Q, Q z1 = z1 (Q as a 2n x 2n matrix, row-major vector)"""
        N2 = self.N2
        A = (M @ self.sinv_pow[k]) % 3
        B = (self.J @ M @ self.J) % 3
        I = np.eye(N2, dtype=np.int64)
        # vec(A Q B) = (A kron B^T) vec(Q) for row-major vec
        L = (np.kron(A, B.T) - np.eye(N2 * N2, dtype=np.int64)) % 3
        Z = np.zeros((N2, N2 * N2), dtype=np.int64)
        for i in range(N2):
            Z[i, i * N2:(i + 1) * N2] = self.z1
        A_all = np.concatenate([L, Z])
        b_all = np.concatenate([np.zeros(N2 * N2, dtype=np.int64), self.z1])
        return solve_affine(A_all, b_all)

    def is_symplectic(self, Q):
        Om = self.wl.Om
        return ((Q.T @ Om @ Q - Om) % 3 == 0).all()

    def good_frames(self, M):
        """the set of reversible frames a for the class M (as a boolean array over all 3^{2n} frames), or None if
        some (S) solution space exceeds CAP"""
        N2 = self.N2
        Minv = R._inv_mod3(M)
        labels = self.wl.labels.astype(np.int64)
        good = np.zeros(len(labels), dtype=bool)
        IminusMinv = (np.eye(N2, dtype=np.int64) - Minv) % 3
        for k in range(3):
            sol = self.symplectic_solutions(M, k)
            if sol is None:
                continue
            x0, basis = sol
            if 3 ** len(basis) > CAP:
                return None
            for coeffs in itertools.product(range(3), repeat=len(basis)):
                q = x0.copy()
                for c, v in zip(coeffs, basis):
                    q = (q + c * v) % 3
                Q = q.reshape(N2, N2)
                if not self.is_symplectic(Q):
                    continue
                sQ = (self.sinv_pow[k] @ Q) % 3
                rQ = self.r_frame(Q)
                # (F): g(e) + s^-k r_Q = M^-1 e - M^-1 a - sQ J a, e_x1 = k.  g(e) = e + t_k is checked below.
                t_k = (self.gframe(np.eye(N2, dtype=np.int64)[0] * k) - np.eye(N2, dtype=np.int64)[0] * k) % 3
                # (I - M^-1) e = -(M^-1 + sQ J) a - t_k - s^-k r_Q   with e_x1 = k  -> for each a solvability
                Lmat = (Minv + sQ @ self.J) % 3
                rhs0 = (-t_k - self.sinv_pow[k] @ rQ) % 3
                # e = k e1 + y, y_x1 = 0:  (I - M^-1) y = rhs0 - (I - M^-1) k e1 - L a
                e1 = np.zeros(N2, dtype=np.int64)
                e1[0] = 1
                base = (rhs0 - IminusMinv @ (k * e1)) % 3
                Y = IminusMinv[:, 1:]                                   # columns for y (y_x1 = 0)
                # solvable iff (base - L a) in col(Y): project with a left nullspace basis of Y
                Rr, piv = rref3(Y.T.copy())                               # row space of Y^T
                ns = solve_affine(Y.T, np.zeros(N2 - 1, dtype=np.int64))
                if ns is None:
                    continue
                _, Nb = ns                                                # w with w^T Y = 0
                if not Nb:
                    good[:] = True
                    continue
                Wm = np.array(Nb)                                         # (m, N2)
                vals = (Wm @ ((base[None, :] - labels @ Lmat.T) % 3).T) % 3   # (m, frames)
                good |= (vals == 0).all(axis=0)
        return good

    def check_g_affine(self, trials=20):
        """g(e) = e + t_{e_x1}: verify on random e"""
        rng = np.random.default_rng(1)
        ok = True
        for _ in range(trials):
            e = self.wl.labels[rng.integers(len(self.wl.labels))].astype(np.int64)
            k = int(e[0])
            t_k = (self.gframe(np.eye(self.N2, dtype=np.int64)[0] * k) - np.eye(self.N2, dtype=np.int64)[0] * k) % 3
            ok &= bool(((self.gframe(e) - e - t_k) % 3 == 0).all())
        return ok


def run():
    import w33_pass11330_orbit_census as O
    res = dict(pass_id=11350)
    D = Decider(2)
    res["g_is_e_plus_t_k"] = D.check_g_affine()
    Ms, keys = O.all_symplectic(D.wl)
    counts = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    agree = undecided = 0
    mism = []
    for i, M in enumerate(Ms):
        g = D.good_frames(M)
        if g is None:
            undecided += 1
            continue
        c = int((~g).sum())
        if c == counts[i]:
            agree += 1
        elif len(mism) < 10:
            mism.append((i, c, int(counts[i])))
    res["two_qutrit_classes"] = len(Ms)
    res["agree_with_pass_11330_counts"] = agree
    res["undecided_cap"] = undecided
    res["mismatches"] = mism
    print(res, flush=True)
    # frame-level spot check against the Weyl criterion
    rng = np.random.default_rng(11350)
    spot = 0
    for _ in range(300):
        i = int(rng.integers(len(Ms)))
        g = D.good_frames(Ms[i])
        if g is None:
            continue
        a = int(rng.integers(81))
        U = D.wl.W[a] @ D.weil(Ms[i]) @ D.T1
        spot += (R.decide(U, 2, rng)[0] is True) == bool(g[a])
    res["frame_spot_checks_agree"] = spot
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
