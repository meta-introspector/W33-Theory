"""Pass 11353: the substrate arrow is visible in level statistics -- reversible magic ticks are orthogonal-class, violating
ones unitary-class, and the symplectic class cannot occur.

Floquet ticks on n qutrits (D = 3^n): W = random Clifford+T circuit (layers of random Clifford words with a cubic gate
on a random qutrit).  Ensembles:
  * violating: U = W (a long magic circuit; reversible with probability -> 0, Pass 11332);
  * reversible: U = W W^T (then K U K = conj(U) = U^-1: reversed by complex conjugation, a Clifford anti-unitary);
  * Clifford only: U = random Clifford (finite order: degenerate spectrum, no level repulsion).
Statistic: the spacing ratio r = min(s_i, s_{i+1}) / max(...) of eigenphases; random-matrix values (Atas et al., PRL 110,
084101 (2013)): Poisson 0.3863, COE 0.5307, CUE 0.5996, CSE 0.6744.
Substrate remark (elementary): an anti-unitary Theta on an odd-dimensional space cannot square to -1 (Kramers pairs
need even dimension), so any substrate time reversal has Theta^2 = +1: reversible magic dynamics is orthogonal-class,
never symplectic-class.  Checked: Theta^2 for the anti-unitary Clifford reversals found by Pass 11252 at n = 2.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11353_level_statistics.json"
RMT = dict(poisson=0.3863, COE=0.5307, CUE=0.5996, CSE=0.6744)


def rbar(U):
    ph = np.sort(np.angle(np.linalg.eigvals(U)))
    s = np.diff(np.concatenate([ph, [ph[0] + 2 * np.pi]]))
    a, b = s, np.roll(s, -1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.minimum(a, b) / np.maximum(a, b)
    r = r[np.isfinite(r)]
    return float(r.mean()), float(np.mean(s < 1e-9))


def magic_circuit(n, rng, layers):
    U = np.eye(3 ** n, dtype=complex)
    for _ in range(layers):
        U = R.random_clifford(n, rng, length=4 * n) @ U
        U = R.local(n, int(rng.integers(n)), P2.T) @ U
    return U


def theta_squared_check(rng, trials=40):
    """for reversible two-qutrit ticks, recover the anti-unitary Clifford reversal V K and record V conj(V)"""
    import w33_pass11235_two_qutrit_t_spectrum as S35
    reps = np.load(P2.CACHE)
    out = []
    for _ in range(trials):
        U = magic_circuit(2, rng, 1)
        vals = S35.fidelity_fast  # noqa: F841
        Pa = P2.PA
        Y = np.einsum('aji,jk,akl->ail', Pa.conj(), U, Pa)
        best = None
        for s in range(0, len(reps), 4000):
            C = reps[s:s + 4000]
            B = np.einsum('rij,jk,rlk->ril', C, U.conj(), C.conj())
            tr = np.abs(np.einsum('rjk,akj->ra', B, Y)) / 9
            idx = np.argwhere(tr > 1 - 1e-9)
            if len(idx):
                r, a = idx[0]
                best = Pa[a] @ C[r]                       # V = P_a C_r with V U* V^dag ~ U^-1
                break
        if best is None:
            continue
        V = best
        sq = V @ V.conj()
        lam = np.trace(sq) / 9
        out.append(complex(lam) if np.allclose(sq, lam * np.eye(9)) else None)
    return out


def run(samples=40):
    rng = np.random.default_rng(11353)
    res = dict(pass_id=11353, rmt_reference=RMT, ensembles={})
    for n, layers in ((4, 12), (5, 15)):
        for name in ("violating", "reversible", "clifford"):
            rs, deg = [], []
            for _ in range(samples):
                if name == "clifford":
                    U = R.random_clifford(n, rng, length=60 * n)
                else:
                    W = magic_circuit(n, rng, layers)
                    U = W if name == "violating" else W @ W.T
                    if name == "reversible":
                        assert np.allclose(U.conj(), np.linalg.inv(U), atol=1e-8)
                r, d = rbar(U)
                rs.append(r)
                deg.append(d)
            res["ensembles"][f"n={n}:{name}"] = dict(r_mean=float(np.mean(rs)), r_stderr=float(np.std(rs) / np.sqrt(len(rs))),
                                                    degenerate_spacing_fraction=float(np.mean(deg)))
            print(n, name, res["ensembles"][f"n={n}:{name}"], flush=True)
    th = theta_squared_check(rng)
    res["theta_squared_values"] = sorted({str(np.round(x, 6)) if x is not None else "non-scalar" for x in th})
    res["theta_squared_checked"] = len(th)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "ensembles"}, indent=1, default=str))


if __name__ == "__main__":
    main()
