"""Pass 11536: the good rules for blocks of dimension 8 -- witnesses in Sp(8, 3) (the method of Pass 11532).

Pass 11511 reduced the good rules to two block statements and proved them in dimensions 2 and 4 by exhaustion:
  (U) for a unipotent u on V1 and every frame a, some anti-symplectic reverser A of u (A u A^-1 = u^-1) has
      omega((A^-1 - I) x, a) = 0 for all x in ker(u - I);
  (E) for M with M + I nilpotent (resp. M^2 + I nilpotent) and every z with M z = -z (resp. M^2 z = -z), some reverser of M
      negates z.
Sp(6, 3) has 3^18 unipotent elements, too many to exhaust.  But each statement is PROVED for a given element by a witness:
the reversers of u are A0 * C(u) (A0 one reverser, C(u) the symplectic centraliser), both obtained by F3 linear algebra
plus a symplecticity filter; random reversers are drawn until their criteria cover all 729 frames (U), resp. some negates z
for every admissible z (E).  This is run on many random elements of each type, tallied by Jordan type.
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

OUT = ROOT / "data" / "w33_pass11536_good_rules_dim8.json"
N2 = 8
I6 = np.eye(N2, dtype=np.int64)
OM = __import__("w33_pass11421_four_qutrits").LightWeyl(4).Om % 3
J6 = np.diag([1, 2] * 4)


def inv3(A):
    return L.R._inv_mod3(np.asarray(A) % 3) % 3


def null(A):
    A = np.asarray(A) % 3
    s = L.solve_affine(A, np.zeros(A.shape[0], np.int64))
    return [np.asarray(b) % 3 for b in s[1]] if s and s[1] else []


def is_symp(Q, sign=1):
    return ((Q.T @ OM @ Q - sign * OM) % 3 == 0).all()


def order(M, cap=4000):
    P = M.copy()
    for k in range(1, cap):
        if (P % 3 == I6).all():
            return k
        P = P @ M % 3
    return None


def unipotent_part(M):
    m = order(M)
    t = 1
    while m % 3 == 0:
        m //= 3
        t *= 3
    a = m * pow(m, -1, t) if t > 1 else 0
    return np.linalg.matrix_power(M, a) % 3 if t > 1 else I6.copy()


def affine_space(rows, rhs):
    sol = L.solve_affine(np.asarray(rows) % 3, np.asarray(rhs) % 3)
    if sol is None:
        return None
    return sol[0] % 3, [np.asarray(b) % 3 for b in sol[1]]


def commutant_rows(Aleft, Bright):
    """vec(X) with Aleft X - X Bright = 0, row-major: (Aleft kron I - I kron Bright^T)"""
    return (np.kron(Aleft, I6) - np.kron(I6, Bright.T)) % 3


def sample(space, rng, test, tries):
    x0, B = space
    for _ in range(tries):
        x = x0.copy()
        for b in B:
            x = (x + int(rng.integers(3)) * b) % 3
        X = x.reshape(N2, N2)
        if test(X):
            return X
    return None


def reversers(u, rng, n=4000):
    """random anti-symplectic reversers of u: A0 * G, G in the symplectic centraliser"""
    uinv = inv3(u)
    rev_space = affine_space(commutant_rows(I6, u) * 0 + (np.kron(I6, u.T) - np.kron(uinv, I6)) % 3, np.zeros(N2 * N2))
    A0 = sample(rev_space, rng, lambda X: is_symp(X, -1), 200000)
    if A0 is None:
        return None
    cen = affine_space(commutant_rows(u, u), np.zeros(N2 * N2))
    out = [A0]
    x0, B = cen
    for _ in range(n * 20):
        if len(out) >= n:
            break
        x = x0.copy()
        for b in B:
            x = (x + int(rng.integers(3)) * b) % 3
        G = x.reshape(N2, N2)
        if is_symp(G):
            out.append(A0 @ G % 3)
    return out


def U_witness(u, rng):
    K = null((u - I6) % 3)
    labs = np.array(np.meshgrid(*[range(3)] * N2, indexing='ij')).reshape(N2, -1).T
    revs = reversers(u, rng)
    if revs is None:
        return None, 0
    cover = np.zeros(len(labs), bool)
    for A in revs:
        Ai = inv3(A)
        ok = np.ones(len(labs), bool)
        for x in K:
            ok &= (labs @ OM @ ((Ai @ x - x) % 3)) % 3 == 0
        cover |= ok
        if cover.all():
            return True, len(revs)
    return False, len(revs)


def E_witness(M, zs, rng):
    revs = reversers(M, rng, n=3000)
    if revs is None:
        return None
    R = np.array(revs)
    return all((((R @ z) + z) % 3 == 0).all(axis=1).any() for z in zs)


def jordan_type(N):
    """dims of ker N^j, j = 1..6 (determines the Jordan type of a nilpotent N)"""
    return tuple(len(null(np.linalg.matrix_power(N, j) % 3)) for j in range(1, N2 + 1))


def form_class(N):
    """for a regular nilpotent N (one Jordan block of size 6): the square class of omega(v, N^5 v), v outside ker N^5"""
    N5 = np.linalg.matrix_power(N, N2 - 1) % 3
    for v in I6:
        if (N5 @ v % 3).any():
            return int(v @ OM @ (N5 @ v) % 3)
    return None


def run(n_samples=4000, per_class=4, seed=11532):
    """DECOMPOSITION.  (U) and (E) pass to orthogonal sums (block-diagonal reversers act blockwise), and every block of
    dimension <= 4 is settled by Pass 11511's exhaustion.  In dimension 6 the only indecomposable cases are:
      in dimension 8 the only NEW indecomposable unipotent class is the regular one, J8 (two classes, told apart by the
      square class of omega(v, N^7 v)); J4 + J4 and smaller pieces split into blocks settled by Passes 11511/11532; for the
      x^2 + 1 type the new indecomposable is s u with u regular unipotent in C(s) = U(4, 9).
    Each is witnessed `per_class` times on random representatives; all other Jordan types are only tallied."""
    rng = np.random.default_rng(seed)
    gens = GEO.gen_mats(4)
    J0 = np.zeros((N2, N2), np.int64)
    for i in range(0, N2, 2):
        J0[i, i + 1], J0[i + 1, i] = 1, 2
    assert is_symp(J0) and ((J0 @ J0 + I6) % 3 == 0).all()
    REG, PAIR = tuple(range(1, N2 + 1)), (2, 4, 6, 8, 8, 8, 8, 8)     # PAIR = regular unipotent of U(4, 9) (x^2+1 type)
    tally, done = Counter(), Counter()
    wit = Counter()
    for t in range(n_samples):
        M = GEO.random_symplectic(rng, gens) % 3
        u = unipotent_part(M)
        jt = jordan_type((u - I6) % 3)
        tally[f"unipotent part {jt}"] += 1
        key = None
        if jt == REG:
            key = f"[8] form class {form_class((u - I6) % 3)}"
        if key and done[key] < per_class:
            done[key] += 1
            ok, n = U_witness(u, rng)
            wit[f"{key}: (U) {'proved' if ok else ('no reverser found' if ok is None else 'NOT witnessed')}"] += 1
            zs = itertools_nonzero(null((u - I6) % 3))
            e = E_witness((-u) % 3, zs, rng)
            wit[f"{key}: (E) for -u {'proved' if e else ('no reverser found' if e is None else 'NOT witnessed')}"] += 1
        if done["U(4) regular"] < per_class and "--no-x2" not in sys.argv:
            g = GEO.random_symplectic(rng, gens) % 3
            s_ = g @ J0 @ inv3(g) % 3
            cen = affine_space(commutant_rows(s_, s_), np.zeros(N2 * N2))
            Gc = sample(cen, rng, is_symp, 20000)
            if Gc is not None:
                Mi = s_ @ unipotent_part(Gc) % 3
                jt2 = jordan_type((Mi @ Mi + I6) % 3)
                tally[f"x^2+1 type {jt2}"] += 1
                if jt2 == PAIR:
                    done["U(4) regular"] += 1
                    zs = itertools_nonzero(null((Mi @ Mi + I6) % 3))
                    e = E_witness(Mi, zs, rng)
                    wit[f"U(4) regular: (E) {'proved' if e else ('no reverser found' if e is None else 'NOT witnessed')}"] += 1
        if (t + 1) % 200 == 0:
            print(t + 1, dict(done), dict(wit), flush=True)
        need = ("[8] form class 1", "[8] form class 2") + (() if "--no-x2" in sys.argv else ("U(4) regular",))
        if all(done[k] >= per_class for k in need):
            break
    return dict(pass_id=11536, samples=t + 1, witnessed=dict(wit), classes_found=dict(done), tally=dict(tally))


def itertools_nonzero(basis):
    import itertools
    out = []
    for co in itertools.product(range(3), repeat=len(basis)):
        if any(co):
            out.append(sum(c * b for c, b in zip(co, basis)) % 3)
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    res = run(int(args[0]) if args else 6000)
    res["x2_plus_1_type"] = ("skipped: the U(4, 9)-regular class did not occur in 800 random samples of the first run; "
                             "open" if "--no-x2" in sys.argv else "searched")
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
