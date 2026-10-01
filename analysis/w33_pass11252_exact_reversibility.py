"""Pass 11252: an exact, complete decision procedure for substrate time reversibility -- and three qutrits decided.

Pass 11239 could only bound three-qutrit reversal fidelities from below: the anti-unitary three-qutrit Clifford group
has ~3.3e12 elements, and generic search is blind (it found 0.41 on an exactly reversible tick).

REDUCTION (exact, no unitaries needed).  Write U = 3^{-n} sum_p c(p) W(p) over the symmetric Weyl operators
W(x,z) = omega^{2xz} X^x Z^z (omega = exp(2 pi i/3)), for which W(p) W(q) = omega^{<p,q>} W(p+q), W(p)^dag = W(-p),
W(x,z)^* = W(x,-z) =: W(Jp), and the Weil representation acts WITHOUT phases: V_M W(p) V_M^dag = W(Mp).  For
V = W(a) V_M, comparing Weyl coefficients in V U^* V^dag = lambda U^dag gives, with L = -M J (anti-symplectic):

    U is substrate-reversible  <=>  there are an anti-symplectic L, b in F_3^{2n} and |mu| = 1 with
                                    c(L p) = mu omega^{<b, L p>} c(p)   for every p.

So reversibility is a combinatorial property of the Weyl coefficients.  Necessary consequences used for pruning:
|c(Lp)| = |c(p)| (magnitudes), and psi(p) = c(Lp)/c(p) is mu times a LINEAR character on every subspace.

ALGORITHM.  Backtracking over the images of a basis (chosen in the smallest magnitude classes), pruning every partial
map on the anti-symplectic form, on the magnitudes of every vector of the partial span, and on the solvability over
F_3 of psi(p)/psi(p0) = omega^{f(p - p0)} with f linear on the partial span.  A map surviving to a leaf satisfies the
criterion on all of F_3^{2n}; the first such map decides 'reversible', and exhausting the search decides 'violating'.
An independent check builds V_M by the Weil twirl V_M ~ sum_p W(Mp) A W(p)^dag and evaluates |tr(V U^* V^dag U)| for
every Weyl coset; the two must agree (asserted).

VALIDATION (independent methods must agree): on two qutrits the verdict is compared with the exhaustive
51840 x 81 search for every word of Passes 11235 and 11238; positive controls: random three-qutrit Cliffords (always
reversible: no Clifford T-violation) and I (x) [(I(x)T)SUM] (reversible by Pass 11239's product reversal).
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
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11252_exact_reversibility.json"
OM = np.exp(2j * np.pi / 3)
X1 = np.roll(np.eye(3), 1, axis=0)
Z1 = np.diag([1, OM, OM ** 2])


class Weyl:
    def __init__(self, n):
        self.n = n
        self.D = 3 ** n
        self.labels = np.array(list(itertools.product(range(3), repeat=2 * n)))      # (x1, z1, x2, z2, ...)
        self.pow3 = 3 ** np.arange(2 * n - 1, -1, -1)
        mats = []
        for lab in self.labels:
            M = np.array([[1.0 + 0j]])
            ph = 0
            for i in range(n):
                x, z = lab[2 * i], lab[2 * i + 1]
                M = np.kron(M, np.linalg.matrix_power(X1, x) @ np.linalg.matrix_power(Z1, z))
                ph += 2 * x * z
            mats.append(OM ** (ph % 3) * M)
        self.W = np.array(mats)
        om = np.zeros((2 * n, 2 * n), dtype=int)
        for i in range(n):
            om[2 * i, 2 * i + 1] = 1
            om[2 * i + 1, 2 * i] = -1
        self.Om = om
        self.J = np.diag([1, -1] * n) % 3

    def index(self, v):
        return int((np.asarray(v) % 3) @ self.pow3)

    def form(self, p, q):
        return int(np.asarray(p) @ self.Om @ np.asarray(q)) % 3

    def coeffs(self, U):
        return np.einsum('pji,ji->p', self.W.conj(), U)          # tr(W(p)^dag U)


def magnitude_classes(c, gap=1e-6):
    """class labels of |c(p)|^2 by gap clustering (equal values never split, unlike rounding at a fixed digit);
    asserts that distinct classes are separated by far more than the float noise"""
    g = np.abs(c) ** 2
    order = np.argsort(g)
    cls = np.empty(len(g), dtype=int)
    k = 0
    cls[order[0]] = 0
    for i in range(1, len(order)):
        d = g[order[i]] - g[order[i - 1]]
        assert d < 1e-9 or d > gap, "ambiguous magnitude gap"
        if d > gap:
            k += 1
        cls[order[i]] = k
    return cls


def anti_symplectic_automorphisms(wl, cls, c=None, tol=1e-7):
    """generator of all anti-symplectic linear L with cls[L p] == cls[p] for all p (lazy; vectorised pruning on every
    vector of the partial span and on the form).  If the Weyl coefficients c are given, also prune on the necessary
    phase condition c(Lp)^3 = nu c(p)^3 with one global nu (= mu^3; nu = 1 when tr U != 0)."""
    n2 = 2 * wl.n
    sizes = Counter(cls)
    order = sorted(range(1, len(wl.labels)), key=lambda i: (sizes[cls[i]], i))
    basis, span = [], {0}
    for i in order:                                  # greedy basis from the smallest magnitude classes
        if i not in span:
            v = wl.labels[i]
            basis.append(v)
            span |= {wl.index(wl.labels[s] + lam * v) for s in span for lam in (1, 2)}
            if len(basis) == n2:
                break
    B = np.array(basis)
    Binv = _inv_mod3(B.T)
    labels = wl.labels
    by_class = {}
    for i, k in enumerate(cls):
        by_class.setdefault(k, []).append(i)
    by_class = {k: labels[np.array(v)] for k, v in by_class.items()}
    lams = [np.array(list(itertools.product(range(3), repeat=k)), dtype=int).reshape(3 ** k, k) for k in range(n2 + 1)]
    OmB = B @ wl.Om                                      # rows: B[i]^T Om

    c3 = None if c is None else np.where(np.abs(c) > tol, c ** 3, 0)
    sup = None if c is None else np.abs(c) > tol
    nu0 = None if c is None or not sup[0] else 1.0 + 0j

    def rec(images, nu):
        k = len(images)
        if k == n2:
            yield (np.array(images).T @ Binv) % 3
            return
        bk = B[k]
        cand = by_class[cls[wl.index(bk)]]
        if k:
            Im = np.array(images)
            # form: <v, image_i> = -<b_k, b_i>
            need = (-(OmB[k] @ B[:k].T)) % 3                       # <b_k, b_i>  -> required -<..>
            have = (cand @ wl.Om @ Im.T) % 3
            cand = cand[(have == need).all(axis=1)]
            if not len(cand):
                return
            L = lams[k]
            src = (L @ B[:k]) % 3                                   # (3^k, n2)
            dst = (L @ Im) % 3
        else:
            src = dst = np.zeros((1, n2), dtype=int)
        S_idx, D_idx = [], []
        for lk in (1, 2):
            s_idx = ((src + lk * bk) % 3) @ wl.pow3                 # (3^k,)
            d_idx = ((dst[None, :, :] + lk * cand[:, None, :]) % 3) @ wl.pow3   # (|cand|, 3^k)
            keep = (cls[d_idx] == cls[s_idx][None, :]).all(axis=1)
            cand = cand[keep]
            S_idx.append(s_idx)
            D_idx = [d[keep] for d in D_idx] + [d_idx[keep]]
            if not len(cand):
                return
        if c3 is None:
            for v in cand:
                yield from rec(images + [v], nu)
            return
        # full necessary phase condition on the partial span: psi(p) = c(Lp)/c(p) = mu omega^{f(p)} with f LINEAR on
        # the span (f(p) = <b, Lp>); equivalently psi(p)/psi(p0) = omega^{f(p - p0)} on the supported span vectors
        coords = lams[k + 1]                                         # (3^{k+1}, k+1): all span vectors
        Bk = np.vstack([B[:k], bk[None, :]]) if k else bk[None, :]
        s_full = ((coords @ Bk) % 3) @ wl.pow3
        msk = sup[s_full]
        A_rows = coords[msk]
        cs = c[s_full[msk]]
        for v in cand:
            Imk = np.vstack([np.array(images), v[None, :]]) if k else v[None, :]
            d_full = ((coords @ Imk) % 3) @ wl.pow3
            psi = c[d_full[msk]] / cs
            kk = np.angle(psi / psi[0]) / (2 * np.pi / 3)
            kr = np.round(kk)
            if not np.allclose(kk, kr, atol=1e-6):
                continue
            if _solvable_mod3((A_rows - A_rows[0]) % 3, kr.astype(int) % 3):
                yield from rec(images + [v], nu)
        return

    yield from rec([], nu0)


def _inv_mod3(A):
    n = len(A)
    M = np.concatenate([A % 3, np.eye(n, dtype=int)], axis=1)
    for c in range(n):
        p = next(r for r in range(c, n) if M[r, c] % 3)
        M[[c, p]] = M[[p, c]]
        M[c] = (M[c] * (1 if M[c, c] == 1 else 2)) % 3
        for r in range(n):
            if r != c and M[r, c]:
                M[r] = (M[r] - M[r, c] * M[c]) % 3
    return M[:, n:]


def weil(wl, Msym, rng):
    """a unitary V with V W(p) V^dag = W(M p) for all p (up to phase)"""
    imgs = np.array([wl.index((Msym @ lab) % 3) for lab in wl.labels])
    for _ in range(5):
        A = rng.normal(size=(wl.D, wl.D)) + 1j * rng.normal(size=(wl.D, wl.D))
        V = np.einsum('pij,jk,plk->il', wl.W[imgs], A, wl.W.conj())
        nrm = np.linalg.norm(V)
        if nrm > 1e-6:
            V = V / nrm * np.sqrt(wl.D)
            if np.allclose(V @ V.conj().T, np.eye(wl.D), atol=1e-8):
                return V
    raise RuntimeError("Weil twirl failed")


def phase_check(wl, c, L, tol=1e-7):
    """is there mu (|mu| = 1) and b in F_3^{2n} with c(L p) = mu omega^{<b, L p>} c(p) for all p?"""
    Lp = (wl.labels @ L.T) % 3
    idx = Lp @ wl.pow3
    cL = c[idx]
    sup = np.abs(c) > tol
    if not np.array_equal(sup, np.abs(cL) > tol):
        return False
    ps = np.flatnonzero(sup)
    psi = cL[ps] / c[ps]
    if not np.allclose(np.abs(psi), 1, atol=1e-6):
        return False
    q = Lp[ps]                                                   # rows: L p
    for j in range(3):
        mu = psi[0] * OM ** (-j)                                 # fixes <b, L p0> = j
        k = np.angle(psi / mu) / (2 * np.pi / 3)
        kr = np.round(k)
        if not np.allclose(k, kr, atol=1e-6):
            continue
        rhs = kr.astype(int) % 3                                 # <b, q> = rhs  with <b, q> = b^T Om q
        A = (q @ wl.Om.T) % 3                                    # row r: (Om q_r)^T, unknown b
        if _solvable_mod3(A, rhs):
            return True
    return False


def _solvable_mod3(A, y):
    M = np.concatenate([A % 3, (y % 3)[:, None]], axis=1).astype(int)
    rows, cols = M.shape
    r = 0
    for c_ in range(cols - 1):
        p = next((i for i in range(r, rows) if M[i, c_] % 3), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        M[r] = (M[r] * (1 if M[r, c_] == 1 else 2)) % 3
        nz = M[:, c_] % 3 != 0
        nz[r] = False
        M[nz] = (M[nz] - M[nz, c_][:, None] * M[r]) % 3
        r += 1
    return not any((M[i, :-1] == 0).all() and M[i, -1] % 3 for i in range(rows))


def decide(U, n, rng=None, limit=5 * 10 ** 6, twirl=False):
    """exact decision.  Returns (verdict, number of magnitude automorphisms examined, best twirl fidelity or None);
    verdict True (reversible), False (violating: ALL magnitude automorphisms examined, none admits the phases),
    None (undecided: enumeration limit hit).  twirl=True also builds V_M and evaluates |tr(V U^* V^dag U)|."""
    rng = rng or np.random.default_rng(0)
    wl = WEYL[n]
    c = wl.coeffs(U)
    cls = magnitude_classes(c)
    best, count = None, 0
    Y = np.einsum('aji,jk,akl->ail', wl.W.conj(), U, wl.W) if twirl else None
    for L in anti_symplectic_automorphisms(wl, cls, c):
        count += 1
        ok = phase_check(wl, c, L)
        if twirl:
            M = (-L @ wl.J) % 3
            VM = weil(wl, M, rng)
            Bm = VM @ U.conj() @ VM.conj().T
            f = float(np.abs(np.einsum('jk,akj->a', Bm, Y)).max() / wl.D)
            best = f if best is None else max(best, f)
            assert ok == (f > 1 - 1e-9), "phase criterion and twirl disagree"
        if ok:
            return True, count, best
        if count >= limit:
            return None, count, best
    return False, count, best


WEYL = {}


def sum_gate(n, ctrl, tgt):
    D = 3 ** n
    M = np.zeros((D, D), dtype=complex)
    for idx in range(D):
        digs = [(idx // 3 ** (n - 1 - k)) % 3 for k in range(n)]
        out = digs.copy()
        out[tgt] = (digs[tgt] + digs[ctrl]) % 3
        M[sum(out[k] * 3 ** (n - 1 - k) for k in range(n)), idx] = 1
    return M


def local(n, q, g):
    ops = [np.eye(3)] * n
    ops[q] = g
    M = np.array([[1.0 + 0j]])
    for o in ops:
        M = np.kron(M, o)
    return M


def random_clifford(n, rng, length=60):
    gates = [local(n, q, g) for q in range(n) for g in (P2.H, P2.S)]
    gates += [sum_gate(n, c, t) for c, t in itertools.permutations(range(n), 2)]
    V = np.eye(3 ** n, dtype=complex)
    for g in rng.integers(len(gates), size=length):
        V = gates[g] @ V
    return V


def random_word(n, d, rng):
    U = random_clifford(n, rng)
    for _ in range(d):
        U = random_clifford(n, rng) @ local(n, int(rng.integers(n)), P2.T) @ U
    return U


def validate_two_qutrit(rng):
    """agreement with the exhaustive 51840 x 81 search"""
    import w33_pass11235_two_qutrit_t_spectrum as S35
    reps = np.load(P2.CACHE)
    agree = disagree = 0
    rv = Counter()
    for d in range(0, 4):
        for _ in range(50):
            U = random_word(2, d, rng)
            brute = S35.fidelity_fast(U, reps) > 1 - 1e-9
            mine, nL, _ = decide(U, 2, rng, twirl=True)
            rv[(d, brute)] += 1
            if brute == mine:
                agree += 1
            else:
                disagree += 1
    T, I3 = P2.T, P2.I3
    named = {
        "(T(x)T)SUM": np.kron(T, T) @ P2.SUM,
        "(I(x)T)SUM": np.kron(I3, T) @ P2.SUM,
        "(I(x)T)SUM(I(x)T^2)": np.kron(I3, T) @ P2.SUM @ np.kron(I3, T @ T),
        "(I(x)T)SUM(I(x)T)": np.kron(I3, T) @ P2.SUM @ np.kron(I3, T),
    }
    named_out = {}
    for k, U in named.items():
        brute = S35.fidelity_fast(U, reps) > 1 - 1e-9
        mine, nL, best = decide(U, 2, rng, twirl=True)
        named_out[k] = dict(exhaustive_reversible=bool(brute), decided_reversible=mine, n_maps_examined=nL)
        if brute == mine:
            agree += 1
        else:
            disagree += 1
    return dict(agree=agree, disagree=disagree, by_depth_and_verdict={f"depth{k[0]}:reversible={k[1]}": v for k, v in rv.items()},
                named=named_out)


def run():
    rng = np.random.default_rng(11252)
    for n in (2, 3):
        WEYL[n] = Weyl(n)
    res = dict(pass_id=11252)
    res["two_qutrit_validation"] = validate_two_qutrit(rng)
    print(res["two_qutrit_validation"], flush=True)
    T, I3 = P2.T, P2.I3
    S12, S23 = sum_gate(3, 0, 1), sum_gate(3, 1, 2)
    kron3 = lambda a, b, c: np.kron(np.kron(a, b), c)  # noqa: E731
    cands = {
        "(T(x)T)SUM (x) I": np.kron(np.kron(T, T) @ P2.SUM, I3),
        "(T(x)T(x)T) SUM12": kron3(T, T, T) @ S12,
        "(I(x)I(x)T) SUM23": kron3(I3, I3, T) @ S23,
        "(T(x)T(x)T) SUM23 SUM12": kron3(T, T, T) @ S23 @ S12,
        "(T(x)I(x)T) SUM23 SUM12": kron3(T, I3, T) @ S23 @ S12,
        "(I(x)I(x)T) SUM23 SUM12": kron3(I3, I3, T) @ S23 @ S12,
        "(T(x)I(x)I) SUM12 SUM23": kron3(T, I3, I3) @ S12 @ S23,
    }
    out = {}
    for k, U in cands.items():
        rev, nL, best = decide(U, 3, rng, twirl=True)
        out[k] = dict(reversible=rev, n_maps_examined=nL, twirl_fidelity_of_found_reversal=best)
        print(k, out[k], flush=True)
    res["three_qutrit_candidates"] = out
    # positive controls with generic magnitude structure: Clifford conjugates C U C^dag of exactly reversible ticks
    # (reversible by C V C^T K).  Plain Cliffords are reversible too (no Clifford T-violation) but their flat magnitude
    # functions have ~|Sp(6,3)| automorphisms, so they are a poor (slow) test of the search.
    base = kron3(I3, I3, T) @ S23
    ctrl = []
    for _ in range(6):
        Cc = random_clifford(3, rng)
        ctrl.append(decide(Cc @ base @ Cc.conj().T, 3, rng)[0])
    res["three_qutrit_conjugated_control_reversible"] = f"{sum(1 for v in ctrl if v)}/{len(ctrl)}"
    print("controls", res["three_qutrit_conjugated_control_reversible"], flush=True)
    # census: random three-qutrit Clifford + k cubic gates
    census = {}
    for d in (1, 2, 3, 4):
        verdicts = [decide(random_word(3, d, rng), 3, rng) for _ in range(100)]
        census[str(d)] = dict(words=len(verdicts), violating=sum(1 for v in verdicts if v[0] is False),
                              reversible=sum(1 for v in verdicts if v[0] is True),
                              undecided=sum(1 for v in verdicts if v[0] is None),
                              max_maps_examined=max(v[1] for v in verdicts))
        print(d, census[str(d)], flush=True)
    res["three_qutrit_census"] = census
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
