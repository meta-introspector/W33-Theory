"""Pass 11499: the "good" half of Pass 11373's four rules -- a sufficient condition for EVERY frame to be reversible when
M z1 = -z1 or M^2 z1 = -z1, valid for every n, and how far it reaches.

Linear-decider form (Passes 11350, 11498): reversal needs a symplectic Q with (S) M s^-k Q (J M J) = Q, Q z1 = z1, and the
frame condition (F).  Take k = 0 (t_0 = 0; r_Q = 0 because Q z1 = z1), A = Q J.

THEOREM.  Let M z1 = -z1 or M^2 z1 = -z1.  If some symplectic Q solves (S) at k = 0 and A = Q J fixes ker(M - I)
pointwise, then every frame a is reversible.
PROOF.  At k = 0, (F) reads (I - M^-1) e = -(M^-1 + A) a with e_x1 = 0.
 (i) Solvability: Im(I - M^-1) = Im(M - I) = ker(M - I)^perp.  For x in ker(M - I): omega(x, M^-1 a) = omega(M x, a) = omega(x, a)
     and omega(x, A a) = -omega(A^-1 x, a) = -omega(x, a); so (M^-1 + A) a is orthogonal to ker(M - I) for every a.
 (ii) The constraint e_x1 = omega(e, z1) = 0 is automatic.  Pick y with (I - M) y = z1: y = -z1 if M z1 = -z1,
     y = -(z1 + M z1) if M^2 z1 = -z1.  Then omega(e, z1) = omega((I - M^-1) e, y) = omega((M^-1 + A) a, -y), and
     A z1 = -z1, A M = M^-1 A (from (S): A = M A M) make the two terms cancel: omega(e, z1) = 0.          QED.
The hypothesis is LINEAR in Q (Q J x = x for x in a basis of ker(M - I)), so existence is decided by Gaussian elimination
plus a symplecticity filter.

Computed: the hypothesis for every good class at n = 2 (all 2592) and every good orbit at n = 3 (Pass 11373), against the
decider's verdicts (every frame reversible).
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

OUT = ROOT / "data" / "w33_pass11499_good_rules.json"


def cell(D, M):
    z = D.z1
    Mz = M @ z % 3
    if ((Mz + z) % 3 == 0).all():
        return "M z1 = -z1"
    if ((M @ Mz + z) % 3 == 0).all():
        return "M^2 z1 = -z1"
    return None


def hypothesis_Q(D, M, cap=3 ** 10):
    """a symplectic Q solving (S) at k = 0 with Q J x = x on ker(M - I), or None; 'capped' if the space is too big"""
    N2 = D.N2
    A0 = M % 3
    B = (D.J @ M @ D.J) % 3
    rows = [(np.kron(A0, B.T) - np.eye(N2 * N2, dtype=np.int64)) % 3]
    rhs = [np.zeros(N2 * N2, dtype=np.int64)]
    Z = np.zeros((N2, N2 * N2), dtype=np.int64)
    for i in range(N2):
        Z[i, i * N2:(i + 1) * N2] = D.z1
    rows.append(Z)
    rhs.append(D.z1.copy())
    ker = L.solve_affine((M - np.eye(N2, dtype=np.int64)) % 3, np.zeros(N2, np.int64))[1]
    for x in ker:
        x = np.asarray(x) % 3
        y = (D.J @ x) % 3
        Y = np.zeros((N2, N2 * N2), dtype=np.int64)
        for i in range(N2):
            Y[i, i * N2:(i + 1) * N2] = y
        rows.append(Y)
        rhs.append(x)
    sol = L.solve_affine(np.concatenate(rows) % 3, np.concatenate(rhs) % 3)
    if sol is None:
        return None, len(ker)
    x0, basis = sol
    if 3 ** len(basis) > cap:
        return "capped", len(ker)
    for co in itertools.product(range(3), repeat=len(basis)):
        q = x0.copy()
        for c, b in zip(co, basis):
            q = (q + c * np.asarray(b)) % 3
        Q = q.reshape(N2, N2)
        if D.is_symplectic(Q):
            return Q, len(ker)
    return None, len(ker)


def check(D, Ms, label, verdicts=True):
    st = Counter()
    mass = Counter()
    for item in Ms:
        M, size = (item if isinstance(item, tuple) else (item, 1))
        M = np.asarray(M) % 3
        c = cell(D, M)
        if c is None:
            continue
        st[c] += 1
        Q, dk = hypothesis_Q(D, M)
        key = "capped" if isinstance(Q, str) else ("hypothesis holds" if Q is not None else "hypothesis fails")
        st[f"{c}: {key}"] += 1
        mass[f"{c}: {key}"] += size
        mass[c] += size
        if verdicts:
            g = D.good_frames(M)
            if g is not None:
                st[f"{c}: decided"] += 1
                st[f"{c}: decided all frames reversible"] += bool(g.all())
                if key == "hypothesis holds":
                    st[f"{c}: hypothesis holds and all frames reversible"] += bool(g.all())
    frac = {k: v / mass[k.split(":")[0]] for k, v in mass.items() if ":" in k}
    out = dict(set=label, counts=dict(st), mass_fraction_within_cell=frac)
    print(out, flush=True)
    return out


def union_k0(D, M, cap=3 ** 10):
    """frame-level form of the theorem: a is reversible through a k = 0 solution Q whenever omega((A^-1 - I) x, a) = 0 for
    all x in ker(M - I) (A = Q J); step (ii) of the proof never used the fixing hypothesis.  Returns the union over all
    k = 0 solutions, or None if capped."""
    import w33_pass11498_magic_axis_law_all_n as T
    lab = D.wl.labels.astype(np.int64)
    N2 = D.N2
    sols = T.solutions(D, M, 0, cap)
    if sols is None:
        return None
    ker = [np.asarray(x) % 3 for x in L.solve_affine((M - np.eye(N2, dtype=np.int64)) % 3, np.zeros(N2, np.int64))[1]]
    cover = np.zeros(len(lab), bool)
    for Q in sols:
        A = (Q @ D.J) % 3
        Ainv = L.R._inv_mod3(A) % 3 if hasattr(L, "R") else np.round(np.linalg.inv(A)).astype(np.int64) % 3
        ws = [((Ainv @ x - x) % 3) for x in ker]
        ok = np.ones(len(lab), bool)
        for wv in ws:
            ok &= (lab @ D.wl.Om @ wv) % 3 == 0
        cover |= ok
    return cover


def run_union():
    import w33_pass11330_orbit_census as O
    L.CAP = 3 ** 11
    D = L.Decider(2)
    st = Counter()
    for M in np.array(O.all_symplectic(D.wl)[0]) % 3:
        c = cell(D, M)
        if c is None:
            continue
        Q, dk = hypothesis_Q(D, M)
        if Q is not None and not isinstance(Q, str):
            continue
        st[f"{c}: hypothesis fails (dim ker(M-I) = {dk})"] += 1
        cov = union_k0(D, M)
        if cov is None:
            I = np.eye(D.N2, dtype=np.int64)
            identity_solves = bool(((M @ I @ (D.J @ M @ D.J)) % 3 == I).all()) and dk == 0
            st[f"{c}: capped; ker(M-I) = 0 and Q = I solves (S) (hypothesis vacuous, theorem applies)"] += identity_solves
            st[f"{c}: capped and unresolved"] += not identity_solves
            continue
        st[f"{c}: union of k=0 criteria covers every frame"] += bool(cov.all())
        g = D.good_frames(M)
        if g is not None:
            st[f"{c}: union frames all reversible (sound)"] += bool(g[cov].all())
    print(dict(st), flush=True)
    return dict(st)


def run():
    import w33_pass11330_orbit_census as O
    import w33_pass11457_magic_axis_theorem as TH
    import w33_pass11373_three_qutrit_exact_fraction as X
    L.CAP = 3 ** 11
    res = dict(pass_id=11499)
    D2 = L.Decider(2)
    res["n2"] = check(D2, list(np.array(O.all_symplectic(D2.wl)[0]) % 3), "n=2 all classes")
    D3 = L.Decider(3)
    res["n3"] = check(D3, [(M, size) for M, size in X.read_orbits(3)], "n=3 all orbits (Pass 11373)")
    return res


def main():
    if "union" in sys.argv:
        res = json.load(open(OUT))
        res["n2_failures_by_union"] = run_union()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
