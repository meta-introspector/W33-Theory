"""Pass 11498: THE MAGIC-AXIS LAW FOR EVERY NUMBER OF QUTRITS -- a five-line proof in the linear decider's own variables.

Setting (Passes 11331, 11350): one magic gate, U = W(a) V_M T1 on n qutrits, z1 the magic axis.  U is reversible iff for
some k in F3 a symplectic Q solves
    (S)  M s^-k Q (J M J) = Q,   Q z1 = z1,
and the frame condition
    (F)  g(e) + s^-k r_Q = M^-1 (e - a) - s^-k Q J a,   e_x1 = k
has a solution e.  Here J = diag(1, -1, ...), s is the unit shear (x1, z1) -> (x1, z1 + x1), g(e) = e + t_k is the frame of
T1 W(e) T1^-1 V_{s^-k}^dag and r_Q the frame of T1 V_Q T1^-1 V_Q^dag (Pass 11350, validated on all 51,840 two-qutrit
classes).

THEOREM.  If M^2 z1 = z1 and v = z1 + M z1 != 0, every frame a with omega(v, a) != 0 is violating (Pass 11420's law F').
PROOF.
 (1) omega(M z1, z1) = omega(M^2 z1, M z1) = omega(z1, M z1) = -omega(M z1, z1), so omega(M z1, z1) = 0.  Since
     omega(y, z1) = y_x1, both z1 and M z1 have x1-coordinate 0: s fixes v, and omega(v, z1) = 0.
 (2) A := Q J is anti-symplectic and (S) reads A = M s^-k A M, i.e. A M = s^k M^-1 A.  With J z1 = -z1, Q z1 = z1 and
     M^-1 z1 = M z1:  A z1 = -z1,  A M z1 = s^k M^-1 A z1 = -s^k M z1 = -M z1.  Hence A v = -v.
 (3) Q z1 = z1 means V_Q W(z1) V_Q^dag = W(z1); T1 is diagonal on qutrit 1, a function of W(z1); so V_Q commutes with T1
     and r_Q = 0.  And t_k lies in span(z1) (T X^k T^-1 = X^k x a quadratic phase on qutrit 1).
 (4) Apply omega(v, .) to (F).  M and s^k fix v, so omega(v, M^-1 y) = omega(v, s^-k y) = omega(v, y); omega(v, t_k) = 0
     by (1); r_Q = 0 by (3).  What remains is  omega(v, a) + omega(v, A a) = 0.  A anti-symplectic with A v = -v gives
     omega(v, A a) = -omega(A^-1 v, a) = omega(v, a).  So 2 omega(v, a) = 0, i.e. omega(v, a) = 0.        QED.

Computed here: each ingredient ((1) on random classes, (2) on every enumerated solution of (S), (3) r_Q and t_k) for
n = 1..4, and the conclusion against the decider's reversible frames (n = 2 all classes, n = 3 orbit representatives,
n = 4 Stab(z1) classes the decider reaches).
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402

OUT = ROOT / "data" / "w33_pass11498_magic_axis_law_all_n.json"


def decider(n):
    return L.Decider(n) if n <= 3 else F4.SparseDecider(n)


def lemma_tk(D):
    N2 = D.N2
    e1 = np.zeros(N2, np.int64)
    e1[0] = 1
    out = []
    for k in range(3):
        t = (D.gframe(e1 * k) - e1 * k) % 3
        out.append(bool((t[np.arange(N2) != 1] == 0).all()))
    return all(out)


def lemma_rQ(D, n, rng, want=30, tries=40000):
    st = Counter()
    gens = GEO.gen_mats(n) if n > 1 else None
    for _ in range(tries):
        if st["tested"] >= want:
            break
        Q = GEO.random_symplectic(rng, gens) % 3 if n > 1 else np.array([[1, 0], [int(rng.integers(3)), 1]])
        if ((Q @ D.z1 - D.z1) % 3).any():
            continue
        st["tested"] += 1
        st["r_Q = 0"] += bool((D.r_frame(Q) % 3 == 0).all())
    return dict(st)


def solutions(D, M, k, cap):
    sol = D.symplectic_solutions(M, k)
    if sol is None:
        return []
    x0, basis = sol
    if 3 ** len(basis) > cap:
        return None
    out = []
    for co in itertools.product(range(3), repeat=len(basis)):
        q = x0.copy()
        for c, b in zip(co, basis):
            q = (q + c * np.asarray(b)) % 3
        Q = q.reshape(D.N2, D.N2)
        if D.is_symplectic(Q):
            out.append(Q)
    return out


def check_classes(D, Ms, cap=3 ** 8, verdicts=True):
    z = D.z1
    lab = D.wl.labels.astype(np.int64)
    st = Counter()
    for M in Ms:
        M = np.asarray(M) % 3
        if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
            continue
        v = (z + M @ z) % 3
        st["classes"] += 1
        st["omega(Mz1, z1) = 0"] += int((M @ z) @ D.wl.Om @ z % 3 == 0)
        allsol, capped = [], False
        for k in range(3):
            s = solutions(D, M, k, cap)
            if s is None:
                capped = True
                break
            allsol += s
        if capped:
            st["solution space beyond cap"] += 1
        else:
            st["solutions enumerated"] += len(allsol)
            st["solutions with QJv = -v"] += sum(bool(((Q @ D.J @ v + v) % 3 == 0).all()) for Q in allsol)
        if verdicts:
            g = D.good_frames(M)
            if g is not None:
                st["decided"] += 1
                st["decided: every reversible frame has omega(v,a) = 0"] += bool(((lab[g] @ D.wl.Om @ v) % 3 == 0).all())
    return dict(st)


def linear_check(D, Ms):
    """step (2) on the WHOLE affine solution space of (S) (symplecticity is not used by the proof): Q J v = -v holds for
    every Q = x0 + sum c_i b_i iff x0 J v = -v and b_i J v = 0 -- no enumeration, any n"""
    z = D.z1
    N2 = D.N2
    st = Counter()
    for M in Ms:
        M = np.asarray(M) % 3
        if ((M @ M @ z - z) % 3).any() or not ((z + M @ z) % 3).any():
            continue
        v = (z + M @ z) % 3
        st["classes"] += 1
        Jv = (D.J @ v) % 3
        for k in range(3):
            sol = D.symplectic_solutions(M, k)
            if sol is None:
                st["k with empty (S)"] += 1
                continue
            x0, basis = sol
            st["affine solution spaces"] += 1
            st["dimension sum"] += len(basis)
            ok = bool(((x0.reshape(N2, N2) @ Jv + v) % 3 == 0).all()) and all(
                ((np.asarray(b).reshape(N2, N2) @ Jv) % 3 == 0).all() for b in basis)
            st["spaces with QJv = -v throughout"] += ok
    return dict(st)


def run_linear():
    import w33_pass11330_orbit_census as O
    import w33_pass11457_magic_axis_theorem as TH
    out = {}
    D2 = L.Decider(2)
    out["n2"] = linear_check(D2, O.all_symplectic(D2.wl)[0])
    D3 = L.Decider(3)
    out["n3"] = linear_check(D3, [m for m, _ in TH.n3_reps()])
    D4 = F4.SparseDecider(4)
    out["n4"] = linear_check(D4, [m for m, _ in TH.n4_classes()])
    print(out, flush=True)
    return out


def run():
    import w33_pass11330_orbit_census as O
    import w33_pass11457_magic_axis_theorem as TH
    rng = np.random.default_rng(11498)
    res = dict(pass_id=11498, lemmas={})
    for n in (1, 2, 3, 4):
        D = decider(n)
        res["lemmas"][str(n)] = dict(t_k_in_span_z1=lemma_tk(D), r_Q=lemma_rQ(D, n, rng, want=30 if n < 4 else 4))
        print(n, res["lemmas"][str(n)], flush=True)
    L.CAP = 3 ** 11
    D2 = L.Decider(2)
    res["n2_all_classes"] = check_classes(D2, O.all_symplectic(D2.wl)[0])
    print(res["n2_all_classes"], flush=True)
    D3 = L.Decider(3)
    res["n3_orbit_reps"] = check_classes(D3, [m for m, _ in TH.n3_reps()])
    print(res["n3_orbit_reps"], flush=True)
    D4 = F4.SparseDecider(4)
    res["n4_stab_z1_classes"] = check_classes(D4, [m for m, _ in TH.n4_classes()], cap=3 ** 6, verdicts=False)
    print(res["n4_stab_z1_classes"], flush=True)
    return res


def main():
    if "linear" in sys.argv:
        res = json.load(open(OUT))
        res["linear_check_whole_solution_space"] = run_linear()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
