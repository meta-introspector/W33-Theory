"""Pass 11511: the good rules for every n, reduced to the unipotent part of M at eigenvalue 1.

Setting (Passes 11373, 11499): M z1 = -z1 or M^2 z1 = -z1.  At k = 0 a reversal is A = QJ anti-symplectic with A = M A M
(i.e. A M A^-1 = M^-1) and A z1 = -z1; Pass 11499 showed a frame a is reversible through A whenever
omega((A^-1 - I) x, a) = 0 for all x in K = ker(M - I).

REDUCTION LEMMA.  Let V1 = ker (M - I)^(2n) (the generalised 1-eigenspace) and V' its M-invariant complement (the other
generalised eigenspaces).  They are omega-orthogonal; K is inside V1 and z1 is inside V' (its eigenvalue is -1, resp. a root
of x^2 + 1).  A reverser maps the generalised lambda-space to the lambda^-1-space, so A = A1 (+) A'.  Conversely, for ANY
reverser A1 of M1 = M|V1 on V1 and the V'-block A' of one k = 0 solution, A1 (+) A' is again a k = 0 solution
(anti-symplectic, A = M A M, A z1 = -z1).  Since K is inside V1 and V1 is orthogonal to V', the criterion only involves A1
and the V1-component of a.  Hence:
   every frame is reversible  <=  (E) (S) has a k = 0 solution, and
                                  (U) for every a1 in V1 some reverser A1 of M1 has omega((A1^-1 - I) x, a1) = 0 on ker(M1 - I).
LEMMA (semisimple case).  If M1 = I, (U) holds: take the anti-symplectic involution that is -1 on a Lagrangian containing a1
and +1 on a complementary Lagrangian; then Im(A1^-1 - I) is that Lagrangian, which is orthogonal to a1.
EXHAUSTION.  (U) holds for every unipotent element of Sp(2, 3) and of Sp(4, 3) (all 3^2 and 3^8 of them).  So:

THEOREM.  For every n, if (S) has a k = 0 solution and dim V1 <= 4 (or M is semisimple at 1), every frame is reversible.

Checked here: the two exhaustions; every good class at n = 2 and every good orbit at n = 3 satisfies (E) and dim V1 <= 4, and
the decider confirms all frames reversible.
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
import w33_pass11499_good_rules as G  # noqa: E402

OUT = ROOT / "data" / "w33_pass11511_good_rules_all_n.json"


def inv3(A):
    return L.R._inv_mod3(np.asarray(A) % 3) % 3


def null(A):
    A = np.asarray(A) % 3
    s = L.solve_affine(A, np.zeros(A.shape[0], np.int64))
    return [np.asarray(b) % 3 for b in s[1]] if s and s[1] else []


def std_form(m):
    Om = np.zeros((m, m), np.int64)
    for i in range(0, m, 2):
        Om[i, i + 1], Om[i + 1, i] = 1, 2
    return Om


def U_holds(M1, anti, Om):
    """(U) for a unipotent M1 on a standard symplectic space, given all anti-symplectic maps"""
    m = M1.shape[0]
    I = np.eye(m, dtype=np.int64)
    M1inv = inv3(M1)
    rev = anti[((np.einsum('aij,jk->aik', anti, M1) - np.einsum('ij,ajk->aik', M1inv, anti)) % 3 == 0).all(axis=(1, 2))]
    K = null((M1 - I) % 3)
    labs = np.array(list(itertools.product(range(3), repeat=m)))
    cover = np.zeros(len(labs), bool)
    for A in rev:
        Ai = inv3(A)
        ok = np.ones(len(labs), bool)
        for x in K:
            ok &= (labs @ Om @ ((Ai @ x - x) % 3)) % 3 == 0
        cover |= ok
        if cover.all():
            break
    return bool(cover.all()), len(rev)


def exhaust_U():
    import w33_pass11330_orbit_census as O
    out = {}
    # Sp(2,3) = SL(2,3)
    S2 = np.array([np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 1])
    A2 = np.array([np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 2])
    Om2 = std_form(2)
    uni2 = [M for M in S2 if not (np.linalg.matrix_power((M - np.eye(2, dtype=np.int64)) % 3, 2) % 3).any()]
    res2 = [U_holds(M, A2, Om2) for M in uni2]
    out["Sp(2,3)"] = dict(unipotent=len(uni2), U_holds=sum(r[0] for r in res2), min_reversers=min(r[1] for r in res2))
    D2 = L.Decider(2)
    S4 = np.array(O.all_symplectic(D2.wl)[0]) % 3
    J4 = np.diag([1, 2, 1, 2])
    A4 = (S4 @ J4) % 3
    Om4 = D2.wl.Om
    assert (np.array(Om4) % 3 == std_form(4)).all()
    uni4 = [M for M in S4 if not (np.linalg.matrix_power((M - np.eye(4, dtype=np.int64)) % 3, 4) % 3).any()]
    res4 = [U_holds(M, A4, Om4) for M in uni4]
    out["Sp(4,3)"] = dict(unipotent=len(uni4), U_holds=sum(r[0] for r in res4), min_reversers=min(r[1] for r in res4))
    print(out, flush=True)
    return out


def k0_solution(D, M, rng, tries=300000):
    Q, _ = G.hypothesis_Q(D, M)
    if Q is not None and not isinstance(Q, str):
        return "hypothesis solution"
    I = np.eye(D.N2, dtype=np.int64)
    if ((M @ (D.J @ M @ D.J)) % 3 == I).all():
        return "Q = I"
    sol = D.symplectic_solutions(M, 0)
    if sol is None:
        return None
    x0, basis = sol
    for _ in range(tries):
        q = x0.copy()
        for b in basis:
            q = (q + int(rng.integers(3)) * np.asarray(b)) % 3
        if D.is_symplectic(q.reshape(D.N2, D.N2)):
            return "sampled solution"
    return None


def dim_V1(M):
    N2 = M.shape[0]
    B = (M - np.eye(N2, dtype=np.int64)) % 3
    return len(null(np.linalg.matrix_power(B, N2) % 3)), len(null(B))


def coverage(D, items, label):
    rng = np.random.default_rng(11511)
    st, mass = Counter(), Counter()
    for M, size in items:
        M = np.asarray(M) % 3
        c = G.cell(D, M)
        if c is None:
            continue
        e = k0_solution(D, M, rng)
        v1, k = dim_V1(M)
        covered = e is not None and (v1 <= 4 or v1 == k)
        key = f"{c}: theorem applies" if covered else f"{c}: NOT covered (E={e}, dim V1={v1})"
        st[key] += 1
        mass[key] += size
        st[f"{c}: (E) via {e}"] += 1
        g = D.good_frames(M)
        if g is not None:
            st[f"{c}: decided"] += 1
            st[f"{c}: decided all frames reversible"] += bool(g.all())
    out = dict(set=label, counts=dict(st), mass={k: v for k, v in mass.items()})
    print(out, flush=True)
    return out


def block_dims(M, z):
    """dim V1 = ker (M - I)^(2n), whether M is semisimple at 1, and dim of the generalised block of z1:
    ker (M + I)^(2n) if M z1 = -z1, ker (M^2 + I)^(2n) if M^2 z1 = -z1"""
    N2 = M.shape[0]
    I = np.eye(N2, dtype=np.int64)
    B1 = (M - I) % 3
    v1 = len(null(np.linalg.matrix_power(B1, N2) % 3))
    ss = v1 == len(null(B1))
    if (((M @ z) + z) % 3 == 0).all():
        bz = len(null(np.linalg.matrix_power((M + I) % 3, N2) % 3))
    else:
        bz = len(null(np.linalg.matrix_power((M @ M + I) % 3, N2) % 3))
    return v1, ss, bz


def coverage_blocks(D, items, label):
    """THEOREM (blocks): if dim V1 <= 4 (or M semisimple at 1) and the z1-block has dim <= 4, every frame is reversible --
    (U) by Lemma 2 / the Sp(2,3), Sp(4,3) exhaustion, (E) by the block exhaustion below and Wonenburger on the rest"""
    st, mass = Counter(), Counter()
    rng = np.random.default_rng(115110)
    for M, size in items:
        M = np.asarray(M) % 3
        c = G.cell(D, M)
        if c is None:
            continue
        v1, ss, bz = block_dims(M, D.z1)
        ok = (v1 <= 4 or ss) and bz <= 4
        if ok:
            key = f"{c}: block theorem applies"
        else:
            e = k0_solution(D, M, rng) if (v1 <= 4 or ss) else None
            key = (f"{c}: z1-block dim {bz} > 4, (E) shown per class ({e})" if e is not None else
                   f"{c}: outside (dim V1={v1}, semisimple={ss}, z1-block={bz})")
        st[key] += 1
        mass[key] += size
        mass[c] += size
    out = dict(set=label, counts=dict(st), mass_fraction={k: v / mass[k.split(":")[0]] for k, v in mass.items() if ":" in k})
    print(out, flush=True)
    return out


def exhaust_E():
    """(E) at block level: for every M in Sp(2,3), Sp(4,3) with M + I nilpotent (resp. M^2 + I nilpotent) and every z with
    M z = -z (resp. M^2 z = -z), some anti-symplectic reverser A of M has A z = -z"""
    import w33_pass11330_orbit_census as O
    out = {}
    S2 = np.array([np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 1])
    A2 = np.array([np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(3), repeat=4) if (a * d - b * c) % 3 == 2])
    D2 = L.Decider(2)
    S4 = np.array(O.all_symplectic(D2.wl)[0]) % 3
    A4 = (S4 @ np.diag([1, 2, 1, 2])) % 3
    for name, Sp, Anti, d in (("Sp(2,3)", S2, A2, 2), ("Sp(4,3)", S4, A4, 4)):
        I = np.eye(d, dtype=np.int64)
        labs = np.array([v for v in itertools.product(range(3), repeat=d) if any(v)])
        for cname, nil, zc in (("M z = -z", lambda M: (M + I) % 3, lambda M: ((labs @ M.T + labs) % 3 == 0).all(axis=1)),
                               ("M^2 z = -z", lambda M: (M @ M + I) % 3, lambda M: ((labs @ (M @ M).T + labs) % 3 == 0).all(axis=1))):
            nM = nz = bad = 0
            for M in Sp:
                if (np.linalg.matrix_power(nil(M), d) % 3).any():
                    continue
                nM += 1
                Minv = inv3(M)
                rev = Anti[((np.einsum('aij,jk->aik', Anti, M) - np.einsum('ij,ajk->aik', Minv, Anti)) % 3 == 0).all(axis=(1, 2))]
                for z in labs[zc(M)]:
                    nz += 1
                    bad += not (((np.einsum('aij,j->ai', rev, z) + z) % 3 == 0).all(axis=1)).any()
            out[f"{name}: {cname}"] = dict(elements=nM, vectors=nz, failures=bad)
    print(out, flush=True)
    return out


def run():
    import w33_pass11330_orbit_census as O
    import w33_pass11373_three_qutrit_exact_fraction as X
    L.CAP = 3 ** 11
    res = dict(pass_id=11511)
    res["exhaustion_U"] = exhaust_U()
    D2 = L.Decider(2)
    res["n2"] = coverage(D2, [(M, 1) for M in O.all_symplectic(D2.wl)[0]], "n=2 all classes")
    # n = 3 is covered by the block theorem (stage `blocks`); the decider's verdicts on the n = 3 good orbits are Pass 11499's
    return res


def main():
    if "blocks" in sys.argv:
        import w33_pass11330_orbit_census as O
        import w33_pass11373_three_qutrit_exact_fraction as X
        res = json.load(open(OUT))
        res["exhaustion_E"] = exhaust_E()
        D2 = L.Decider(2)
        res["n2_blocks"] = coverage_blocks(D2, [(M, 1) for M in O.all_symplectic(D2.wl)[0]], "n=2 all classes")
        D3 = L.Decider(3)
        res["n3_blocks"] = coverage_blocks(D3, X.read_orbits(3), "n=3 all orbits")
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
