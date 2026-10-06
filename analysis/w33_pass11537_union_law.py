"""Pass 11537: the union law -- one-gate reversibility is linear algebra on the fixed vectors orthogonal to the magic axis.

Setting: U = W(a) V_M T1 on n qutrits; (S) M s^-k Q (J M J) = Q, Q z1 = z1, Q symplectic (Pass 11350); A_Q = Q J.
Let K0 = { w : M w = w, omega(w, z1) = 0 } and N_Q = (I - A_Q^-1) K0.

THEOREM (necessity, every n).  If frame a is reversible through (Q, k), then a is omega-orthogonal to N_Q.
PROOF.  Pair the frame condition (F) with omega(w, .) for w in K0: M and the shear s^k fix w (w_x1 = omega(w, z1) = 0),
omega(w, t_k) = 0 (t_k in span z1), r_Q = 0 (Q z1 = z1); what remains is omega(w, a) + omega(w, A a) = 0, i.e.
omega((I - A^-1) w, a) = 0.  QED.  So the reversible frames lie in the UNION over all solutions Q of N_Q^perp, and a class
with no solution of (S) has no reversible frame.

UNION LAW (observed): the reversible frames are EXACTLY the union over Q of N_Q^perp.
It contains the magic-axis law (Pass 11498: v = z1 + M z1 is anti-fixed by every A), the one-gate W-law of Pass 11533 (a
universally anti-fixed v), the good-cell criterion of Pass 11499 (there K0 = ker(M - I)), and it explains the exceptions of
the W-law (Pass 11533): their violating sets are affine pieces cut out by two vectors of K0, i.e. intersections of the N_Q^perp.
WHAT SUFFICIENCY NEEDS.  (F) is solvable iff its right-hand side pairs correctly with every w in K' = {w : (I - M) w in
span z1}.  Since K^perp = Im(M - I), either some fixed vector has omega(w, z1) != 0 (z1 not in Im(M - I)), or K' = K + span(w0)
with (I - M) w0 = z1; in both cases one extra, k-dependent affine condition appears, which the union over (Q, k) must absorb.

Checked: every class at n = 2 (exhaustive) and every orbit of Pass 11373 at n = 3 whose solution spaces are enumerable.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11498_magic_axis_law_all_n as T  # noqa: E402

OUT = ROOT / "data" / "w33_pass11537_union_law.json"


def null(A):
    A = np.asarray(A) % 3
    s = L.solve_affine(A, np.zeros(A.shape[0], np.int64))
    return [np.asarray(b) % 3 for b in s[1]] if s and s[1] else []


def union_prediction(D, M, cap=3 ** 9):
    """the union law's reversible frames, or None if a solution space is too large to enumerate"""
    N2 = D.N2
    z = D.z1
    Om = D.wl.Om
    lab = D.wl.labels.astype(np.int64)
    K0 = null(np.concatenate([(M - np.eye(N2, dtype=np.int64)) % 3, (z @ Om)[None, :] % 3]))
    sols = []
    for k in range(3):
        s = T.solutions(D, M, k, cap)
        if s is None:
            return None
        sols += s
    pred = np.zeros(len(lab), bool)
    for Q in sols:
        Ai = L.R._inv_mod3((Q @ D.J) % 3) % 3
        ok = np.ones(len(lab), bool)
        for w in K0:
            ok &= (lab @ Om @ ((w - Ai @ w) % 3)) % 3 == 0
        pred |= ok
        if pred.all():
            break
    return pred


def check(D, items, label):
    st, mass = Counter(), Counter()
    for M, size in items:
        M = np.asarray(M) % 3
        g = D.good_frames(M)
        if g is None:
            st["decider undecided"] += 1
            continue
        pred = union_prediction(D, M)
        if pred is None:
            st["solution space beyond cap"] += 1
            continue
        cell = GEO.cell(M, D.z1, D.wl.Om)
        if (pred == g).all():
            key = f"{cell}: EXACT"
        elif not (g & ~pred).any():
            key = f"{cell}: union strictly larger (sufficiency fails)"
        else:
            key = f"{cell}: VIOLATES the necessity theorem"
        st[key] += 1
        mass[key] += size
    tot = sum(mass.values())
    out = dict(set=label, counts=dict(st), checked_mass_fraction={k: v / tot for k, v in mass.items()})
    print(out, flush=True)
    return out


def per_solution_n2(step=7):
    """is the union law exact for each SINGLE solution (Q, k)?  i.e. are the frames reversible through (Q, k) exactly
    N_Q^perp?  (every 7th class of Sp(4,3); all solutions of each)"""
    import w33_pass11252_exact_reversibility as R
    import w33_pass11330_orbit_census as O
    D = L.Decider(2)
    z, Om, N2 = D.z1, D.wl.Om, 4
    lab = D.wl.labels.astype(np.int64)

    def frames_single(M, k, Q):
        Minv = R._inv_mod3(M)
        IminusMinv = (np.eye(N2, dtype=np.int64) - Minv) % 3
        sQ = (D.sinv_pow[k] @ Q) % 3
        rQ = D.r_frame(Q)
        e1 = np.zeros(N2, np.int64)
        e1[0] = 1
        t_k = (D.gframe(e1 * k) - e1 * k) % 3
        Lmat = (Minv + sQ @ D.J) % 3
        base = ((-t_k - D.sinv_pow[k] @ rQ) % 3 - IminusMinv @ (k * e1)) % 3
        Y = IminusMinv[:, 1:]
        _, Nb = L.solve_affine(Y.T, np.zeros(N2 - 1, dtype=np.int64))
        if not Nb:
            return np.ones(len(lab), bool)
        Wm = np.array(Nb)
        return ((Wm @ ((base[None, :] - lab @ Lmat.T) % 3).T) % 3 == 0).all(axis=0)

    st = Counter()
    for M in (np.array(O.all_symplectic(D.wl)[0]) % 3)[::step]:
        K0 = null(np.concatenate([(M - np.eye(N2, dtype=np.int64)) % 3, (z @ Om)[None, :] % 3]))
        cell = GEO.cell(M, z, Om)
        for k in range(3):
            sols = T.solutions(D, M, k, 3 ** 8)
            for Q in (sols or []):
                Ai = L.R._inv_mod3((Q @ D.J) % 3) % 3
                pred = np.ones(len(lab), bool)
                for w in K0:
                    pred &= (lab @ Om @ ((w - Ai @ w) % 3)) % 3 == 0
                act = frames_single(M, k, Q)
                key = "per-solution EXACT" if (pred == act).all() else (
                    "frames of (Q,k) a proper subset of N_Q-perp" if not (act & ~pred).any() else "NECESSITY VIOLATED")
                st[f"{cell}: {key}"] += 1
    print(dict(st), flush=True)
    return dict(st)


def run():
    import w33_pass11330_orbit_census as O
    import w33_pass11373_three_qutrit_exact_fraction as X
    L.CAP = 3 ** 11
    res = dict(pass_id=11537)
    D2 = L.Decider(2)
    res["n2"] = check(D2, [(M, 1) for M in O.all_symplectic(D2.wl)[0]], "n=2 all classes")
    D3 = L.Decider(3)
    res["n3"] = check(D3, X.read_orbits(3), "n=3 all orbits")
    return res


def main():
    if "perq" in sys.argv:
        res = json.load(open(OUT))
        res["n2_per_solution"] = per_solution_n2()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    json.dump(run(), open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
