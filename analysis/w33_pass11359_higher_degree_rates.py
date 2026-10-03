"""Pass 11359: contraction of the Clifford-averaged cubic gate on higher-degree PU(3) tensor spaces (matrix-free).

Pass 11354 found the second modulus rho = 3/8 on V^(x)p (x) Vbar^(x)q for p + q <= 7.  The reversible fraction of
Pass 11312 decays about 0.72 per gate, so higher degrees must carry slower modes.  Here, matrix-free:
  * pi(C) acts on a tensor with p + q legs by C on p legs and conj(C) on q legs (det-1 lifts);
  * the Clifford-invariant subspace is the range of P = avg_C pi(C) applied to m random vectors, m >= its dimension,
    which is known exactly from characters: dim = avg_C tr(C)^p conj(tr(C))^q;
  * pi(T) is diagonal; rho is the largest modulus < 1 of B^dag pi(T) B (eigenvalues of modulus 1 = SU(3)-invariants).
Degrees: (4,4), (6,3), (9,0), (5,5) (dimension up to 3^10 = 59,049).  Scope: finitely many degrees.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402

OUT = ROOT / "data" / "w33_pass11359_higher_degree_rates.json"


def det1(U):
    return U / np.linalg.det(U) ** (1 / 3)


def apply_pi(C, X, p, q):
    """X has shape (3,)*(p+q) + (m,)"""
    Cc = C.conj()
    for leg in range(p + q):
        M = C if leg < p else Cc
        X = np.moveaxis(np.tensordot(M, X, axes=([1], [leg])), 0, leg)
    return X


def invariant_basis(p, q, Cs, rng):
    dimc = np.mean([np.trace(C) ** p * np.conj(np.trace(C)) ** q for C in Cs])
    k = int(round(dimc.real))
    m = k + 6
    shape = (3,) * (p + q) + (m,)
    X = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    S = np.zeros_like(X)
    for C in Cs:
        S += apply_pi(C, X, p, q)
    S /= len(Cs)
    Mtx = S.reshape(-1, m)
    U, s, _ = np.linalg.svd(Mtx, full_matrices=False)
    r = int(np.sum(s > 1e-8 * s[0]))
    return U[:, :r], k


def run():
    rng = np.random.default_rng(11359)
    Cs = [det1(C) for C in P1.clifford1()]
    T = det1(P1.T1)
    tdiag = np.diag(T)
    res = dict(pass_id=11359, reps={})
    for p, q in ((3, 3), (4, 4), (6, 3), (9, 0), (5, 5)):
        B, k = invariant_basis(p, q, Cs, rng)
        # pi(T) diagonal: phases t[i1]...t[ip] conj(t[j1])...
        ph = np.ones((3,) * (p + q), complex)
        for leg in range(p + q):
            v = tdiag if leg < p else tdiag.conj()
            sh = [1] * (p + q)
            sh[leg] = 3
            ph = ph * v.reshape(sh)
        A = B.conj().T @ (ph.reshape(-1, 1) * B)
        ev = np.linalg.eigvals(A)
        mods = np.sort(np.abs(ev))[::-1]
        n1 = int(np.sum(np.abs(mods - 1) < 1e-7))
        rest = mods[n1:]
        res["reps"][f"({p},{q})"] = dict(dim=3 ** (p + q), clifford_invariants_character=k, basis_rank=int(B.shape[1]),
                                         su3_invariants=n1, rho=float(rest[0]) if len(rest) else 0.0,
                                         top_moduli_below_one=[round(float(x), 6) for x in rest[:6]])
        print((p, q), res["reps"][f"({p},{q})"], flush=True)
    res["sup_rho"] = max(v["rho"] for v in res["reps"].values())
    res["empirical_reversible_decay"] = 0.72
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "reps"}, indent=1))


if __name__ == "__main__":
    main()
