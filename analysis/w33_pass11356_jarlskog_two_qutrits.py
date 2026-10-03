"""Pass 11356: the substrate Jarlskog invariant J_6 on two qutrits -- completeness against exact verdicts.

J_2t(U) = avg_C |tr(U^dag C U C^dag)|^(2t) - avg_C |tr(U^dag C U^T C^dag)|^(2t), averaged over the full two-qutrit
Clifford group (51,840 symplectic classes x 81 Pauli frames; Pass 11355 for the definition; J_2 == 0 for any 2-design, and J_4
vanishing is checked here on the named ticks, not derived from a design property -- see Pass 11357).  Tested against exact verdicts:
  * one-gate ticks W(a) V_M (T (x) I): verdicts from the F3-linear decider (Pass 11350), half drawn from bad classes;
  * multi-gate words (2-4 cubic gates, random qutrits): verdicts from the Weyl criterion (Pass 11252);
  * named ticks: (T (x) T) SUM (violating), (I (x) T) SUM (reversible), (I (x) T) SUM (I (x) T^2) (violating, F_min^2).
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11356_jarlskog_two_qutrits.json"
TOL = 1e-9
_S = {}


def _init():
    R.WEYL[2] = R.Weyl(2)
    _S["reps"] = np.load(P2.CACHE)


def J(U, t=3):
    reps = _S["reps"]
    PA = P2.PA
    Ud = U.conj().T
    tot_a = tot_b = 0.0
    for s in range(0, len(reps), 6000):
        V = reps[s:s + 6000]
        for W, acc in ((U, "a"), (U.T, "b")):
            B = np.einsum('rij,jk,rlk->ril', V, W, V.conj())                 # V W V^dag
            # tr(U^dag P B P^dag) for all Paulis P: (P B P^dag)_{kj} -> einsum
            tr = np.einsum('jk,akl,rlm,amj->ra', Ud, PA, B, PA.conj().transpose(0, 2, 1), optimize=True)
            val = float((np.abs(tr) ** (2 * t)).sum())
            if acc == "a":
                tot_a += val
            else:
                tot_b += val
    n = len(reps) * 81
    return tot_a / n - tot_b / n


def _one_gate(args):
    seed, want_bad = args
    import w33_pass11350_linear_decider as L
    import w33_pass11330_orbit_census as O
    D = _S.setdefault("D", L.Decider(2))
    if "Ms" not in _S:
        _S["Ms"], _ = O.all_symplectic(D.wl)
        _S["counts"] = np.load(ROOT / "data" / "w33_pass11330_class_counts.npy")
    rng = np.random.default_rng(seed)
    pool = np.flatnonzero(_S["counts"] > 0) if want_bad else np.arange(51840)
    i = int(rng.choice(pool))
    M = _S["Ms"][i]
    g = D.good_frames(M)
    if g is None:
        return None
    a = int(rng.choice(np.flatnonzero(~g))) if (want_bad and (~g).any()) else int(rng.integers(81))
    U = D.wl.W[a] @ D.weil(M) @ D.T1
    j = J(U)
    return dict(kind="one-gate", violating=not bool(g[a]), J6=j)


def _multi(seed):
    rng = np.random.default_rng(seed)
    k = int(rng.integers(2, 5))
    U = R.random_clifford(2, rng)
    for _ in range(k):
        U = R.random_clifford(2, rng) @ R.local(2, int(rng.integers(2)), P2.T) @ U
    v = R.decide(U, 2, rng)[0]
    return dict(kind=f"{k}-gate", violating=v is False, J6=J(U))


def run(n_one=120, n_multi=120):
    res = dict(pass_id=11356)
    with Pool(11, initializer=_init) as pool:
        one = [r for r in pool.map(_one_gate, [(100 + i, i % 2 == 0) for i in range(n_one)], chunksize=2) if r]
        multi = pool.map(_multi, [5000 + i for i in range(n_multi)], chunksize=2)
    T, I3 = P2.T, P2.I3
    _init()
    named = {"(T(x)T)SUM": np.kron(T, T) @ P2.SUM, "(I(x)T)SUM": np.kron(I3, T) @ P2.SUM,
             "(I(x)T)SUM(I(x)T^2)": np.kron(I3, T) @ P2.SUM @ np.kron(I3, T @ T)}
    res["named"] = {k: dict(J6=J(U), J4=J(U, 2), J2=J(U, 1)) for k, U in named.items()}
    rows = one + multi
    mism = [r for r in rows if r["violating"] != (r["J6"] > TOL)]
    res["checked"] = len(rows)
    res["violators"] = sum(r["violating"] for r in rows)
    res["mismatches"] = len(mism)
    res["mismatch_examples"] = mism[:5]
    res["max_J6_reversible"] = max((r["J6"] for r in rows if not r["violating"]), default=0.0)
    res["min_J6_violator"] = min((r["J6"] for r in rows if r["violating"]), default=None)
    print(res, flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
