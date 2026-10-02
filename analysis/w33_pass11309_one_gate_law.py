"""Pass 11309: the exact law for one magic gate on two qutrits -- a bad set of symplectic classes times a linear form.

Pass 11266 (one qutrit, exhaustive): C T violates substrate time reversal iff the symplectic part of C is a unit shear
fixing Z (3 of 24) and the Pauli frame shifts X (6 of 9): probability exactly (3/24)(2/3) = 1/12.  For two qutrits it
sampled 0.0912 and found the naive fixed-axis rule fails.

STRUCTURE (found here, then used): for C = P_a V_M (Pauli frame a in F_3^4, symplectic class M), the verdict for
U = C (T (x) I) depends on a only through ONE linear form: either no frame violates (M 'good') or exactly the 54 frames
with l_M(a) + c_M != 0 violate, for an AFFINE form (so 0, 54 or 81 frames).  Consequences:
  * a class is bad iff one of the 5 frames 0, e_1..e_4 violates (a nonconstant affine form cannot vanish on all five);
  * every frame of every bad class is then decided, giving the exact count of violating Cliffords.
Verified: the 0-or-54 dichotomy on every frame of a random sample of classes; |Bad| from all 51840 classes.
The bad set is then tested against simple symplectic predicates involving the magic axis z_1.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11309_one_gate_law.json"
_S = {}
UNIT = [27, 9, 3, 1]          # PA index of the unit frames (a1,a2,a3,a4) = e_i  (PA built as product over (a,b,c,d))


def _init():
    R.WEYL[2] = R.Weyl(2)
    _S["reps"] = np.load(P2.CACHE)
    _S["T1"] = np.kron(P2.T, P2.I3)


def _bad(chunk):
    reps, T1 = _S["reps"], _S["T1"]
    rng = np.random.default_rng(chunk[0])
    out = []
    for r in chunk:
        v = [R.decide(P2.PA[a] @ reps[r] @ T1, 2, rng)[0] is False for a in [0] + UNIT]
        out.append((r, v))
    return out


def _full(r):
    reps, T1 = _S["reps"], _S["T1"]
    rng = np.random.default_rng(r)
    return r, [R.decide(P2.PA[a] @ reps[r] @ T1, 2, rng)[0] is False for a in range(81)]


def symplectic_part(wl, C):
    M = np.zeros((4, 4), int)
    for j in range(4):
        e = np.zeros(4, int)
        e[j] = 1
        img = C @ wl.W[wl.index(e)] @ C.conj().T
        M[:, j] = wl.labels[int(np.argmax(np.abs(np.einsum('pij,ij->p', wl.W.conj(), img))))]
    return M


def run():
    rng = np.random.default_rng(11309)
    res = dict(pass_id=11309)
    with Pool(11, initializer=_init) as pool:
        sample = [int(x) for x in rng.choice(51840, 300, replace=False)]
        full = pool.map(_full, sample, chunksize=4)
        counts = sorted({sum(v) for _, v in full})
        res["dichotomy_sample_classes"] = len(sample)
        res["violating_frames_per_class_values"] = counts
        # affine-form check: violating frames = {a : l(a) + c != 0} for some linear l and constant c (l = 0 allowed)
        lin_ok = True
        labels = np.array([(a // 27, a // 9 % 3, a // 3 % 3, a % 3) for a in range(81)])
        forms = [(l, c) for l in labels for c in range(3)]
        for _, v in full:
            viol = set(map(tuple, labels[np.array(v)]))
            lin_ok &= any(set(map(tuple, labels[(labels @ l + c) % 3 != 0])) == viol for l, c in forms)
        res["violators_are_complement_of_an_affine_hyperplane"] = bool(lin_ok)
        print(res, flush=True)
        chunks = [list(range(s, min(s + 480, 51840))) for s in range(0, 51840, 480)]
        allv = [x for part in pool.map(_bad, chunks) for x in part]
    bad = [r for r, v in allv if any(v)]
    with Pool(11, initializer=_init) as pool:
        fullbad = dict(pool.map(_full, bad, chunksize=8))           # every frame of every bad class
    n_viol = {r: sum(v) for r, v in fullbad.items()}
    res["bad_classes"] = len(bad)
    res["bad_class_frame_counts"] = {str(k): sum(1 for x in n_viol.values() if x == k) for k in sorted(set(n_viol.values()))}
    labels = np.array([(a // 27, a // 9 % 3, a // 3 % 3, a % 3) for a in range(81)])
    forms = [(l, c) for l in labels for c in range(3)]
    res["all_bad_classes_affine"] = all(any(set(map(tuple, labels[(labels @ l + c) % 3 != 0])) ==
                                            set(map(tuple, labels[np.array(v)])) for l, c in forms)
                                        for v in fullbad.values())
    def is_affine(S):
        """S (set of frames) closed under x + y - z, i.e. an affine subspace"""
        S = [np.array(s) for s in S]
        keys = {tuple(s) for s in S}
        return all(tuple((x + y - z) % 3) in keys for x in S for y in S for z in S)
    codim = {}
    for r, v in fullbad.items():
        good = [tuple(labels[a]) for a in range(81) if not v[a]]
        ok = (len(good) == 0) or is_affine(good)
        dim = {0: "empty", 1: 0, 3: 1, 9: 2, 27: 3}.get(len(good), "?")
        codim[str(dim)] = codim.get(str(dim), 0) + 1
        if not ok:
            codim["NOT_AFFINE"] = codim.get("NOT_AFFINE", 0) + 1
    res["nonviolating_frames_affine_subspace_dims"] = codim
    res["law"] = "C(T(x)I) violates iff the Pauli frame avoids an affine subspace A_M of F_3^4 (A_M = F_3^4 for good classes)"
    total = sum(n_viol.values())
    res["violating_cliffords_exact"] = total
    res["probability_exact"] = str(Fraction(total, 51840 * 81))
    res["probability_float"] = total / (51840 * 81)
    print(res["bad_classes"], res["bad_class_frame_counts"], res["probability_exact"], res["probability_float"], flush=True)
    bad54 = {r for r, k in n_viol.items() if k == 54}
    bad81 = {r for r, k in n_viol.items() if k == 81}
    # predicates on the symplectic part, magic axis z1 = (0,1,0,0) (x1, z1, x2, z2)
    wl = R.Weyl(2)
    reps = np.load(P2.CACHE)
    Ms = [symplectic_part(wl, reps[r]) for r in range(51840)]
    badset = set(bad)
    z1 = np.array([0, 1, 0, 0])
    x1 = np.array([1, 0, 0, 0])
    preds = {
        "M z1 = z1": lambda M: tuple(M @ z1 % 3) == tuple(z1),
        "M z1 = +-z1": lambda M: tuple(M @ z1 % 3) in (tuple(z1), tuple(2 * z1 % 3)),
        "<z1, M z1> = 0": lambda M: wl.form(z1, M @ z1 % 3) == 0,
        "M^-1 z1 = z1": lambda M: tuple(R._inv_mod3(M) @ z1 % 3) == tuple(z1),
        "<x1, M z1> = 0": lambda M: wl.form(x1, M @ z1 % 3) == 0,
        "M z1 has zero qutrit-2 part": lambda M: tuple((M @ z1 % 3)[2:]) == (0, 0),
        "M z1 in span(z1, qutrit 2)": lambda M: (M @ z1 % 3)[0] == 0,
    }
    pr = {}
    for name, f in preds.items():
        S = {r for r in range(51840) if f(Ms[r])}
        pr[name] = dict(size=len(S), bad_subset_of=badset <= S, equals_bad=S == badset, overlap=len(S & badset),
                        equals_bad81=S == bad81, equals_bad54=S == bad54, overlap81=len(S & bad81))
    res["predicates"] = pr
    print(pr, flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
