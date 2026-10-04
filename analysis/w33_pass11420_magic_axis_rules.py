"""Pass 11420: the universal magic-axis rules -- what is proved, and what is verified at n = 2, 3, 4.

Pass 11373 found four rules exact at n = 2 and 3 (and Pass 11421 tests them at n = 4): with z1 the magic axis,
    Mz1 = z1  => violating frames exist;      Mz1 = -z1 => every frame reversible;
    same line, M^2 z1 = z1 => violating;       same line, M^2 z1 = -z1 => every frame reversible.

PROVED here (elementary):
  L1. Mz1 = +-z1  <=>  V_M T1 V_M^dag = T1^(+-1).  (V_M Z1 V_M^dag = W(Mz1) and T1 = f(Z1), f(w^k) = zeta^(k^3) is odd
      in k, so f(Z1^-1) = T1^-1.)
  L2. If Mz1 = z1 then U^3 is Clifford for every frame a.  (V_M commutes with T1; for a_x1 = 0 everything commutes and
      T1^3 = Z1; for a_x1 = c != 0, U cycles the Z1-eigenspaces k -> k+c -> k+2c and the cubic phases add to
      sum_j (k + jc)^3 = 3k^3 + 6kc^2 = 0 (mod 9) after zeta-exponentiation with c^2 = 1.)  So the magic of a fixed-
      axis tick is a cube root of a Clifford.
  L3. If Mz1 = -z1 and dim ker(M + I) = 1, there is an anti-symplectic involution sigma with sigma M sigma = M^-1 and
      sigma z1 = -z1, hence a symplectic solution Q = sigma J of (S) at k = 0.  (Wonenburger: some anti-symplectic
      involution tau reverses M; tau z1 lies in ker(M^-1 + I) = ker(M + I), so tau z1 = +-z1; use -tau if needed.)
  L4. C_H(M) fixes z1 and M z1, hence preserves omega(a, z1 + M z1): the frame orbits used for the slow n = 3
      classes respect the predicate of F'.

VERIFIED (computer):
  F'. CONJECTURE F' -- if M^2 z1 = z1 (the two 'bad' cells: Mz1 = z1, and Mz1 on a line through z1 with M^2 z1 = z1),
     then v = z1 + M z1 is FIXED by M and EVERY frame with omega(a, v) != 0 violates (for Mz1 = z1, v = -z1 and the
     condition is a_x1 != 0; for Mz1 = -z1, v = 0, consistent with that cell being all good):
       n = 2: all classes of both cells; n = 3: every orbit of both cells (Pass 11373; slow ones by orbit).
  dim ker(M + I) = 1 is counted on the reversed cell (where L3 applies unconditionally).
"""

from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11373_three_qutrit_exact_fraction as X  # noqa: E402

OUT = ROOT / "data" / "w33_pass11420_magic_axis_rules.json"


def kernel_dim_mod3(A):
    Rr, piv = L.rref3(A % 3)
    return A.shape[1] - len(piv)


def lemma_checks(n=2):
    """L1 and L2 numerically on every relevant n = 2 class (a sanity check of the proofs)"""
    import w33_pass11330_orbit_census as O
    D = L.Decider(n)
    Ms, _ = O.all_symplectic(D.wl)
    z, T1, W = D.z1, D.T1, D.wl.W
    l1 = l2 = 0
    rng = np.random.default_rng(0)
    for M in Ms:
        s = 1 if ((M @ z - z) % 3 == 0).all() else (-1 if ((M @ z + z) % 3 == 0).all() else 0)
        if not s:
            continue
        V = D.weil(M)
        l1 += np.allclose(V @ T1 @ V.conj().T, T1 if s == 1 else T1.conj().T)
        if s == 1:
            for a in rng.choice(len(W), 4, replace=False):
                U = W[a] @ V @ T1
                U3 = U @ U @ U
                ok = all(np.isclose(np.abs(np.einsum('pij,ij->p', W.conj(), U3 @ W[D.wl.index(e)] @ U3.conj().T)).max()
                                    / D.wl.D, 1) for e in np.eye(2 * n, dtype=int))
                l2 += ok
    return dict(L1_checked_classes=int(l1), L2_frames_with_U3_clifford=int(l2))


def conjecture_f_n2():
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    z, lab = D.z1, D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(0)
    out = {}
    for cell in ("Mz1 = z1", "same line, M^2 z1 = z1"):
        st = Counter()
        for M in Ms:
            if GEO.cell(M, z, D.wl.Om) != cell:
                continue
            g = D.good_frames(M)
            if g is None:
                V = D.weil(M)
                g = np.array([R.decide(D.wl.W[a] @ V @ D.T1, 2, rng)[0] is True for a in range(len(lab))])
            bad = ~g
            xne = (lab @ D.wl.Om @ ((z + M @ z) % 3)) % 3 != 0            # omega(a, z1 + M z1) != 0
            st["bad == {omega(a, z1+Mz1) != 0}" if (bad == xne).all() else
               ("bad contains it" if (bad >= xne).all() else "FAILS")] += 1
        out[cell] = dict(st)
    return out


_S = {}


def _init():
    L.CAP = 3 ** 11
    _S["D"] = L.Decider(3)


def _orbit_verdict(args):
    i, M, a = args
    D = _S["D"]
    v = R.decide(D.wl.W[a] @ D.weil(M) @ D.T1, 3, np.random.default_rng(a))[0]
    assert v is not None
    return i, a, v


def conjecture_f_n3(nproc=8):
    """every n = 3 orbit (Pass 11373) in the two 'bad' rule cells: are all frames with a_x1 != 0 violating?"""
    rows = X.read_orbits(3)
    cen = X.read_centralisers(3)
    D = L.Decider(3)
    L.CAP = 3 ** 11
    lab = D.wl.labels.astype(np.int64)
    out = {}
    jobs, plan = [], {}
    fast = Counter()
    for i, (M, size) in enumerate(rows):
        cell = GEO.cell(M, D.z1, D.wl.Om)
        if cell not in ("Mz1 = z1", "same line, M^2 z1 = z1"):
            continue
        xne = (lab @ D.wl.Om @ ((D.z1 + M @ D.z1) % 3)) % 3 != 0      # omega(a, z1 + M z1) != 0
        g = D.good_frames(M) if i not in cen else None
        if g is not None:
            bad = ~g
            fast[(cell, "bad == {omega(a, z1+Mz1) != 0}" if (bad == xne).all() else
                  ("bad contains it" if (bad >= xne).all() else "FAILS"))] += 1
            continue
        orbits, _ = X.frame_orbits(D, cen[i][1])
        targets = [orb for orb in orbits if xne[orb[0]]]
        assert all(xne[np.array(orb)].all() for orb in targets)          # L4: orbits respect omega(a, z1 + M z1)
        plan[i] = (cell, targets)
        jobs += [(i, M, orb[0]) for orb in targets]
    with Pool(nproc, initializer=_init) as pool:
        verdict = {(i, a): v for i, a, v in pool.imap_unordered(_orbit_verdict, jobs, chunksize=1)}
    slow = Counter()
    for i, (cell, targets) in plan.items():
        slow[(cell, "all frames with omega(a, z1+Mz1) != 0 violate" if all(verdict[(i, orb[0])] is False for orb in targets)
              else "FAILS")] += 1
    out["fast_classes"] = {f"{c} | {k}": v for (c, k), v in fast.items()}
    out["slow_classes_by_frame_orbits"] = {f"{c} | {k}": v for (c, k), v in slow.items()}
    return out


def conjecture_f_n4(per_cell=20, seed=114204):
    """F' at n = 4 with Pass 11421's sparse decider on targeted classes of the two 'bad' cells"""
    import w33_pass11421_four_qutrits as F4
    F4._init()
    D, gens = F4._S["D"], F4._S["gens"]
    z, Om = D.z1, D.wl.Om
    lab = D.wl.labels.astype(np.int64)
    rng = np.random.default_rng(seed)
    out = Counter()
    for kind in ("Mz1 = z1", "same line, M^2 z1 = z1"):
        done = 0
        while done < per_cell:
            if kind == "Mz1 = z1":
                M = GEO.random_symplectic(rng, gens, length=300)
                M = (GEO.map_to((M @ z) % 3, z, Om, rng) @ M) % 3
            else:
                reps = F4._S["line_reps"]["M^2 z1 = z1"]
                A = np.eye(8, dtype=np.int64)
                A[:4, :4] = reps[int(rng.integers(len(reps)))]
                A[4:, 4:] = GEO.random_symplectic(rng, GEO.gen_mats(2), length=100)
                h = GEO.random_symplectic(rng, gens, length=300)
                N = (GEO.map_to((h @ z) % 3, z, Om, rng) @ h) % 3
                M = (N @ A @ R._inv_mod3(N)) % 3
            assert GEO.cell(M, z, Om) == kind
            g = D.good_frames(M)
            if g is None:
                out[f"{kind} | undecided"] += 1
                done += 1
                continue
            bad = ~g
            xne = (lab @ Om @ ((z + M @ z) % 3)) % 3 != 0
            out[f"{kind} | " + ("bad == {omega(a, z1+Mz1) != 0}" if (bad == xne).all() else
                                ("bad contains it" if (bad >= xne).all() else "FAILS"))] += 1
            done += 1
    return dict(out)


def kernel_condition():
    import w33_pass11330_orbit_census as O
    D = L.Decider(2)
    Ms, _ = O.all_symplectic(D.wl)
    st = Counter()
    for M in Ms:
        if ((M @ D.z1 + D.z1) % 3 == 0).all():
            st[f"dim ker(M+I) = {kernel_dim_mod3((M + np.eye(4, dtype=np.int64)) % 3)}"] += 1
    rows = X.read_orbits(3)
    D3 = L.Decider(3)
    st3 = Counter()
    for M, size in rows:
        if ((M @ D3.z1 + D3.z1) % 3 == 0).all():
            st3[f"dim ker(M+I) = {kernel_dim_mod3((M + np.eye(6, dtype=np.int64)) % 3)}"] += size
    return dict(n2_classes=dict(st), n3_classes=dict(st3))


def run():
    res = dict(pass_id=11420)
    res["lemma_checks_n2"] = lemma_checks()
    print(res["lemma_checks_n2"], flush=True)
    res["L3_kernel_condition"] = kernel_condition()
    print(res["L3_kernel_condition"], flush=True)
    res["conjecture_F_n2"] = conjecture_f_n2()
    print(res["conjecture_F_n2"], flush=True)
    res["conjecture_F_n3"] = conjecture_f_n3()
    print(res["conjecture_F_n3"], flush=True)
    res["conjecture_F_n4"] = conjecture_f_n4()
    print(res["conjecture_F_n4"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
