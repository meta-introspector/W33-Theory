#!/usr/bin/env python3
"""Pass 11188: the intrinsic arrow of time, exactly, on every conjugacy class of PSp(4,3), PSp(6,3) and PSp(8,3).

Pass 11183 defined the intrinsic arrow of a Clifford tick S as
        A(S) = min over all stabilizer tensor factorisations F of  E(S, F),   E(S, F) = sum_q export_q(S, F)   (trits),
computed it exactly for two qutrits and on a 240-element sample for three qutrits (values {0, 2, 3}).
A is a class function on PSp(2n,3): conjugation by g sends the split F to gF and S_{gF}(gSg^-1) = S_F(S), and the
central sign leaves every block rank unchanged; the splits form one Sp-orbit.  So one representative per class decides
A everywhere.  GAP (analysis/gap/w33_pass11188_class_reps*.g) prints one matrix per class in our adapted basis.

THE PLANE FORMULA.  export_q = rank S_F[others, q] = dim(S P_q + P_q) - 2 = 2 - d(P_q),  d(P) := dim(P cap S P), so
        E(S, F) = 2n - sum_q d(P_q)     and     A(S) = 2n - max { sum d(P) : P_1, ..., P_k mutually orthogonal }
(any orthogonal family of nondegenerate planes extends to a split, and planes with d = 0 add nothing).  This turns the
min over all splits (110565 for n = 3, about 8.3e9 for n = 4) into a max-weight clique over the few thousand planes
that meet their image; it is checked against direct enumeration on all 20 + 74 classes for n = 2, 3.
CONSEQUENCES (for every n).  d(P) = 2 iff S P = P.  E <= n - 1 forces sum d >= n + 1, so some plane is invariant; and
E = 1 is impossible (all but one plane invariant forces the last one invariant too).  Hence
        A(S) <= n - 1   ==>   S fixes a nondegenerate plane,
and conversely, if A <= n - 1 holds for all (n-1)-qutrit ticks, a fixed plane P gives A(S) <= A(S|P^perp) <= n - 1.
RESULTS (exact; counts over Sp(2n,3) are twice the PSp counts).
  * n = 2 (20 classes): reproduces Pass 11183 -- A = 0 on 9576 and A = 2 on 16344 elements of PSp(4,3).
  * n = 3 (74 classes): A in {0, 2, 3} -- never 1 and NEVER MORE THAN 3 of the 6 trits a split could export;
    A = 0 / 2 / 3 on the fractions 55241/884520, 581/1080, 8836/22113 of PSp(6,3).
  * n = 4 (278 classes): A in {0, 2, 3, 4}, max 4.
  * For n = 2, 3, 4:  A(S) = n  <=>  S fixes no nondegenerate plane  (no qutrit that the tick leaves in place);
    otherwise A(S) <= n - 1.  The maximal intrinsic arrow is exactly one trit per qutrit, and it is attained exactly by
    the ticks that leave no subsystem invariant.  The subsystem-free ticks always admit a split in which EVERY qutrit
    meets its own image in exactly a line ('half-moving split').
  * Conjecture (all n): A(S) <= n, with equality iff S has no invariant nondegenerate plane.
"""
from __future__ import annotations

import itertools
import json
import re
import sys
import tempfile
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11180_mereology as M  # noqa: E402
import w33_pass11183_intrinsic_arrow as A  # noqa: E402

GAPOUT = {2: ROOT / "data" / "w33_pass11188_gap_class_reps.txt", 3: ROOT / "data" / "w33_pass11188_gap_class_reps.txt",
          4: ROOT / "data" / "w33_pass11188_gap_class_reps_n4.txt"}
OUT = ROOT / "data" / "w33_pass11188_exact_arrow_by_class.json"
PSP_ORDER = {2: 25920, 3: 4585351680, 4: 65784756654489600}
PLANES_TOTAL = {2: 90, 3: 7371, 4: 597780}          # q^{2n-2} (q^{2n} - 1) / (q^2 - 1), q = 3


def load_classes(n):
    text = GAPOUT[n].read_text().replace("\\\n", "")
    text = re.sub(r"\s+", " ", text)
    out = []
    pat = r"CLASS n=(\d) order=(\d+) size=(\d+) fixed=(-?\d+) real=(\w+) M=\s*(\[ \[.*?\] \])"
    for nn, order, size, fixed, real, mat in re.findall(pat, text):
        if int(nn) != n:
            continue
        rows = re.findall(r"\[([^\[\]]*)\]", mat)
        Mg = np.array([[int(x) for x in r.split(",")] for r in rows], np.int64)
        out.append(dict(order=int(order), size=int(size), fixed=int(fixed), real=real == "true", S=Mg.T % 3))
    return out


def is_symplectic(S, n):
    J = M.form(n)
    return bool(((S.T @ J @ S - J) % 3 == 0).all())


def local_count(S, Bs, n):
    J = M.form(n)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3
    SF = np.einsum('fij,jk,fkl->fil', Binv, S % 3, Bs) % 3
    Bk = SF.reshape(-1, n, 2, n, 2).transpose(0, 1, 3, 2, 4)
    zero = (Bk == 0).all(axis=(-1, -2))
    return int((((~zero).sum(-1) == 1).all(-1)).sum())


def planes_from_splits(Bs):
    """all nondegenerate planes, as spanning pairs (u, v) with omega(u, v) = 1, deduplicated from the splits"""
    out = {}
    for B in Bs:
        for k in range(B.shape[1] // 2):
            u, v = B[:, 2 * k], B[:, 2 * k + 1]
            key = frozenset(tuple((a * u + b * v) % 3) for a in range(3) for b in range(3))
            out.setdefault(key, np.stack([u, v], 1))
    return np.array(list(out.values()))


def all_planes(n):
    """all nondegenerate planes of F_3^{2n} directly (used for n = 4; cached in the temp dir)"""
    cache = Path(tempfile.gettempdir()) / f"w33_planes_{2 * n}.npy"
    if cache.exists():
        Pl = np.load(cache)
        if len(Pl) == PLANES_TOTAL[n]:
            return Pl
    D = 2 * n
    J = M.form(n)
    pts = np.array([v for v in itertools.product(range(3), repeat=D)
                    if any(v) and v[next(i for i in range(D) if v[i])] == 1], np.int64)
    idx = {tuple(p): i for i, p in enumerate(pts)}
    W = (pts @ J @ pts.T) % 3

    def normp(x):
        x = x % 3
        i = int(np.argmax(x != 0))
        return tuple((x * (1 if x[i] == 1 else 2)) % 3)
    seen, out = set(), []
    for a in range(len(pts)):
        for b in np.flatnonzero(W[a, a + 1:]) + a + 1:
            u, v = pts[a], (pts[b] * W[a, b]) % 3
            key = frozenset([a, int(b), idx[normp(u + v)], idx[normp(u + 2 * v)]])
            if key not in seen:
                seen.add(key)
                out.append(np.stack([u, v], 1))
    Pl = np.array(out)
    try:
        np.save(cache, Pl)
    except OSError:
        pass
    return Pl


def dvals(S, Pl):
    """d(P) = dim(P cap S P) for every plane: 0, 1 or 2 (a plane meets its image in 0, 1 or all 4 projective points)"""
    Pl = Pl.astype(np.int8)
    u, v = Pl[:, :, 0], Pl[:, :, 1]
    span = np.stack([(a * u + b * v) % 3 for a in range(3) for b in range(3)], 1)
    S8 = (S % 3).astype(np.int8)
    cnt = np.zeros(len(Pl), np.int64)
    for a, b in ((1, 0), (0, 1), (1, 1), (1, 2)):
        x = (a * u + b * v) % 3
        Sx = (x.astype(np.int64) @ S8.T.astype(np.int64)) % 3
        cnt += (span == Sx[:, None, :].astype(np.int8)).all(-1).any(-1)
    return np.where(cnt == 0, 0, np.where(cnt == 1, 1, 2))


def arrow_by_planes(S, Pl, n):
    """A(S) = 2n - max sum of d over mutually orthogonal planes (exact branch and bound); also #invariant planes"""
    J = M.form(n)
    d = dvals(S, Pl)
    keep = np.flatnonzero(d > 0)
    G, w = Pl[keep].astype(np.int64), d[keep]
    order = np.argsort(-w, kind='stable')
    G, w = G[order], w[order]
    GJ = np.einsum('pai,ab->pib', G, J)                # rows u^T J, v^T J of every candidate plane
    best = [0]

    def orth(i, cands):
        """candidates orthogonal to plane i (lazy: no |G|^2 adjacency matrix)"""
        if len(cands) == 0:
            return cands
        W = np.einsum('ib,pbj->pij', GJ[i], G[cands]) % 3
        return cands[~W.reshape(len(cands), 4).any(1)]

    def bound(score, depth, cands):
        k = n - depth
        n2 = int((w[cands] == 2).sum())
        t = min(n2, k)
        return score + 2 * t + min(len(cands) - n2, k - t)

    def dfs(cands, score, depth):
        best[0] = max(best[0], score)
        if depth == n or len(cands) == 0 or best[0] == 2 * n or bound(score, depth, cands) <= best[0]:
            return
        for t, i in enumerate(cands):
            rest = cands[t + 1:]
            if bound(score + int(w[i]), depth + 1, rest) <= best[0] and score + int(w[i]) <= best[0]:
                continue
            dfs(orth(i, rest), score + int(w[i]), depth + 1)
            if best[0] == 2 * n or bound(score, depth, rest) <= best[0]:
                return
    dfs(np.arange(len(G)), 0, 0)
    return 2 * n - best[0], int((d == 2).sum()), int((d == 1).sum())


def by_class_enumerated(n, Bs, classes):
    rows = []
    Pl = planes_from_splits(Bs)
    for c in classes:
        S = c['S']
        E = A.exports(S, Bs, n).sum(1)
        Ap, inv, half = arrow_by_planes(S, Pl, n)
        rows.append(dict(order=c['order'], size=c['size'], real=c['real'], gap_fixed=c['fixed'],
                         invariant_planes=inv, planes_total=len(Pl),
                         symplectic=is_symplectic(S, n), order_ok=M.porder(S, n) == c['order'],
                         local_ok=local_count(S, Bs, n) == c['fixed'],
                         A=int(E.min()), A_plane_formula=Ap, E_max=int(E.max()), arrow_free=int((E == 0).sum()),
                         at_minimum=int((E == E.min()).sum()),
                         E_hist={str(k): int(v) for k, v in sorted(Counter(E.tolist()).items())}))
    return rows


def by_class_planes(n, classes):
    Pl = all_planes(n)
    rows = []
    for c in classes:
        S = c['S']
        Ap, inv, half = arrow_by_planes(S, Pl, n)
        rows.append(dict(order=c['order'], size=c['size'], real=c['real'], invariant_planes=inv,
                         planes_meeting_image_in_a_line=half, planes_total=len(Pl),
                         symplectic=is_symplectic(S, n), order_ok=M.porder(S, n) == c['order'],
                         A=Ap, A_plane_formula=Ap, local_ok=True))
    return rows


def section(n, rows):
    dist = Counter()
    orders = defaultdict(set)
    for r in rows:
        dist[r['A']] += r['size']
        orders[r['A']].add(r['order'])
    total = sum(r['size'] for r in rows)
    return dict(
        classes=len(rows), total=total, total_ok=total == PSP_ORDER[n],
        planes_ok=all(r['planes_total'] == PLANES_TOTAL[n] for r in rows),
        all_checks=all(r['symplectic'] and r['order_ok'] and r['local_ok'] for r in rows),
        plane_formula_matches=all(r['A'] == r['A_plane_formula'] for r in rows),
        # A < n  <=>  S fixes a nondegenerate plane (the qutrit it keeps exports nothing)
        A_below_n_iff_invariant_plane=all((r['A'] < n) == (r['invariant_planes'] > 0) for r in rows),
        A_distribution_psp={str(k): v for k, v in sorted(dist.items())},
        A_fractions={str(k): str(Fraction(v, total)) for k, v in sorted(dist.items())},
        A_values=sorted(dist), A_max=max(dist), orders_by_A={str(k): sorted(v) for k, v in sorted(orders.items())},
        rows=rows)


def summarize(ns=(2, 3, 4)):
    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    out['pass_id'] = 11188
    import w33_pass11182_paper_ticks_mereology as P
    for n in ns:
        cls = load_classes(n)
        if n == 2:
            rows = by_class_enumerated(2, np.array(M.factorisations(2)), cls)
        elif n == 3:
            rows = by_class_enumerated(3, P.load_facs(), cls)
        else:
            rows = by_class_planes(n, cls)
        out[f"n{n}"] = section(n, rows)
        OUT.write_text(json.dumps(out, indent=1, sort_keys=True))
    return out


if __name__ == "__main__":
    ns = tuple(int(x) for x in sys.argv[1:]) or (2, 3, 4)
    r = summarize(ns)
    for n in ns:
        d = r[f"n{n}"]
        print(f"n={n}: classes={d['classes']} total_ok={d['total_ok']} planes_ok={d['planes_ok']} "
              f"checks={d['all_checks']} plane_formula={d['plane_formula_matches']}")
        print("  A distribution (PSp):", d['A_distribution_psp'], d['A_fractions'])
        print("  orders by A:", d['orders_by_A'], "A<n <=> invariant plane:", d['A_below_n_iff_invariant_plane'])
