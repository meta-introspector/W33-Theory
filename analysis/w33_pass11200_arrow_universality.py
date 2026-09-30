#!/usr/bin/env python3
"""Pass 11200: the arrow law beyond qutrits, and the law of the arrow for large systems.

Pass 11199 proved A(S) = n - c(S) for qutrit Clifford ticks S in Sp(2n,3): the intrinsic arrow (minimal information
exported per tick over all subsystem splits) equals the number of qutrits the tick cannot keep as invariant subsystems.
The proof uses only odd characteristic, so it covers qudits of every odd prime-power dimension q.
UNIVERSALITY (exact, all conjugacy classes, GAP):  A = n - c on every class of
    Sp(4,2), Sp(6,2), Sp(8,2), Sp(10,2)   (qubits, n <= 5 -- characteristic 2, NOT covered by the proof),
    PSp(4,5), PSp(4,7)                    (5- and 7-dits, a check of the odd-q proof),
    PSp(10,3), all 940 classes            (five qutrits: A computed independently of the theorem -- a split attaining
                                           the lower bound -- and compared with the Jordan closed form).
So the law is not special to qutrits; for qubits it is a conjecture verified through five qubits.
BURNSIDE.  Sp(2n,q) is transitive on nondegenerate planes, so the AVERAGE number of planes a tick leaves invariant is
exactly 1 for every n and q (checked from the class data).
LARGE SYSTEMS (qutrits, closed form of Pass 11199 on random ticks up to n = 32): the number c of protected qutrits has
a limiting law with P(c = 0) ~ 0.46, P(c = 1) ~ 0.39, P(c = 2) ~ 0.13, P(c = 3) ~ 0.02 and E[c] ~ 0.73 (exactly
133/180, 0.72532, 0.72552 for n = 2, 3, 4); hence E[A]/n = 1 - E[c]/n -> 1: a typical large Clifford dynamics loses
exactly one trit per qutrit per tick except for a bounded number of protected qutrits, and the maximal arrow A = n
occurs with probability ~ 0.46.
"""
from __future__ import annotations

import itertools
import json
import re
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

GAP_SMALL = ROOT / "data" / "w33_pass11200_gap_class_reps_q.txt"
GAP_BIG = ROOT / "data" / "w33_pass11200_gap_class_reps_q_big.txt"
OUT = ROOT / "data" / "w33_pass11200_arrow_universality.json"


# ---------------- generic F_p symplectic geometry ----------------
def form(n, p):
    J = np.zeros((2 * n, 2 * n), np.int64)
    for k in range(n):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, p - 1
    return J


def points(D, p):
    out = []
    for v in itertools.product(range(p), repeat=D):
        if any(v):
            i = next(k for k in range(D) if v[k])
            if v[i] == 1:
                out.append(v)
    return np.array(out, np.int64)


def all_planes(n, p):
    """every nondegenerate plane, as (u, v) with omega(u, v) = 1 (vectorised per point)"""
    D = 2 * n
    J = form(n, p)
    pts = points(D, p)
    index = {tuple(x): i for i, x in enumerate(pts)}
    W = (pts @ J @ pts.T) % p
    powers = {a: pow(a, -1, p) for a in range(1, p)}
    seen, out = set(), []
    for a in range(len(pts)):
        u = pts[a]
        for b in np.flatnonzero(W[a, a + 1:]) + a + 1:
            v = (pts[b] * powers[int(W[a, b])]) % p
            key = []
            for s in range(p):
                for t in range(p):
                    if (s, t) != (0, 0):
                        x = (s * u + t * v) % p
                        i = int(np.argmax(x != 0))
                        key.append(index[tuple((x * powers[int(x[i])]) % p)])
            key = frozenset(key)
            if key not in seen:
                seen.add(key)
                out.append(np.stack([u, v], 1))
    return np.array(out)


def dvals(S, Pl, p):
    u, v = Pl[:, :, 0], Pl[:, :, 1]
    span = np.stack([(a * u + b * v) % p for a in range(p) for b in range(p)], 1)
    cnt = np.zeros(len(Pl), np.int64)
    for a, b in [(1, 0)] + [(c, 1) for c in range(p)]:
        x = (a * u + b * v) % p
        Sx = (x @ S.T) % p
        cnt += (span == Sx[:, None, :]).all(-1).any(-1)
    return np.where(cnt == 0, 0, np.where(cnt == 1, 1, 2))


def arrow_and_c(S, Pl, n, p):
    """exact A(S) (branch and bound over planes meeting their image) and c(S) (max orthogonal invariant family)"""
    J = form(n, p)
    d = dvals(S, Pl, p)
    keep = np.flatnonzero(d > 0)
    G, w = Pl[keep], d[keep]
    o = np.argsort(-w, kind='stable')
    G, w = G[o], w[o]
    GJ = np.einsum('pai,ab->pib', G, J)

    def orth(i, c, GG, GGJ):
        if len(c) == 0:
            return c
        W = np.einsum('ib,pbj->pij', GGJ[i], GG[c]) % p
        return c[~W.reshape(len(c), 4).any(1)]
    best = [0]

    def bound(s, dep, c):
        k = n - dep
        n2 = int((w[c] == 2).sum())
        t = min(n2, k)
        return s + 2 * t + min(len(c) - n2, k - t)

    def dfs(c, s, dep):
        best[0] = max(best[0], s)
        if dep == n or len(c) == 0 or best[0] == 2 * n or bound(s, dep, c) <= best[0]:
            return
        for t, i in enumerate(c):
            rest = c[t + 1:]
            dfs(orth(i, rest, G, GJ), s + int(w[i]), dep + 1)
            if best[0] == 2 * n or bound(s, dep, rest) <= best[0]:
                return
    dfs(np.arange(len(G)), 0, 0)
    G2 = G[w == 2]
    G2J = np.einsum('pai,ab->pib', G2, J)
    bc = [0]

    def dfs2(c, k):
        bc[0] = max(bc[0], k)
        if bc[0] == n or k + len(c) <= bc[0]:
            return
        for t, i in enumerate(c):
            if k + len(c) - t <= bc[0]:
                return
            dfs2(orth(i, c[t + 1:], G2, G2J), k + 1)
            if bc[0] == n:
                return
    dfs2(np.arange(len(G2)), 0)
    return 2 * n - best[0], bc[0], int((d == 2).sum())


def load_classes(path, q, n):
    text = re.sub(r"\s+", " ", path.read_text().replace("\\\n", ""))
    head = re.search(rf"q={q} n={n} \|PSp\|=(\d+) classes=(\d+)", text)
    out = []
    for qq, nn, order, size, mat in re.findall(r"CLASS q=(\d) n=(\d) order=(\d+) size=(\d+) M=\s*(\[ \[.*?\] \])", text):
        if int(qq) == q and int(nn) == n:
            rows = re.findall(r"\[([^\[\]]*)\]", mat)
            M = np.array([[int(x) for x in r.split(",")] for r in rows], np.int64)
            out.append(dict(order=int(order), size=int(size), S=M.T % q))
    return out, (int(head.group(1)), int(head.group(2))) if head else (None, None)


def class_check(q, n, path):
    cls, (order, ncls) = load_classes(path, q, n)
    Pl = all_planes(n, q)
    J = form(n, q)
    rows, dist, inv_total = [], Counter(), 0
    for c in cls:
        S = c['S']
        assert not ((S.T @ J @ S - J) % q).any()
        A, cc, inv = arrow_and_c(S, Pl, n, q)
        rows.append(dict(order=c['order'], size=c['size'], A=A, c=cc, invariant_planes=inv))
        dist[A] += c['size']
        inv_total += c['size'] * inv
    tot = sum(dist.values())
    return dict(q=q, n=n, classes=len(cls), classes_expected=ncls, group_order=order, total_ok=tot == order,
                planes=len(Pl), law_holds=all(r['A'] == n - r['c'] for r in rows), A_max=max(dist),
                A_distribution={str(k): v for k, v in sorted(dist.items())},
                A_equals_n_fraction=str(Fraction(dist.get(n, 0), tot)),
                burnside_mean_invariant_planes=str(Fraction(inv_total, tot)))


# ---------------- qutrits, n = 5: exact A on all 940 classes of PSp(10,3), independently of the theorem ----------------
import w33_pass11199_arrow_theorem as T  # noqa: E402

def symplectic_split(Nb, Om):
    """columns of Nb span a nondegenerate subspace; return a list of mutually orthogonal nondegenerate planes"""
    vecs = [Nb[:, i] for i in range(Nb.shape[1])]
    planes = []
    B = Nb.copy()
    while B.shape[1] > 0:
        G = (B.T @ Om @ B) % 3
        i, j = map(int, np.argwhere(G != 0)[0])
        u, v = B[:, i], B[:, j]
        P = np.stack([u, v], 1); planes.append(P)
        # project the rest onto P^perp
        Wv = T.nullspace((P.T @ Om) % 3).T                 # P^perp in the full space
        # basis of span(B) cap P^perp
        M = np.concatenate([B, Wv], 1)
        ns = T.nullspace(M)
        inter = (B @ ns[:, :B.shape[1]].T) % 3 if len(ns) else np.zeros((len(Om), 0), np.int64)
        R, piv = T.rref(inter.T) if inter.size else (np.zeros((0, len(Om)), np.int64), [])
        B = R.T % 3
    return planes

def exact_A(S, n, pts):
    Om = T.form_std(n); D = 2 * n; I = np.eye(D, dtype=np.int64)
    fam = []
    for lam in (1, 2):
        E = T.nullspace((S - lam * I) % 3).T                 # eigenspace basis (columns)
        if E.shape[1] == 0: continue
        G = (E.T @ Om @ E) % 3
        rad = T.nullspace(G)                                 # radical coordinates
        # nondegenerate complement: coordinates complementing the radical
        R, piv = T.rref(np.concatenate([rad, np.eye(E.shape[1], dtype=np.int64)], 0)) if len(rad) else (None, None)
        if len(rad):
            Mfull = np.concatenate([rad, np.eye(E.shape[1], dtype=np.int64)], 0).T   # columns: rad..., e_i...
            _, pv = T.rref(Mfull)
            comp = [c - len(rad) for c in pv if c >= len(rad)]
            Nb = E[:, comp]
        else:
            Nb = E
        if Nb.shape[1] >= 2 and ((Nb.T @ Om @ Nb) % 3).any():
            fam += symplectic_split(Nb, Om)
    # W1 = complement of the eigen families
    if fam:
        Pm = np.concatenate(fam, 1); inW = ~(((pts @ Om @ Pm) % 3) != 0).any(1)
    else:
        inW = np.ones(len(pts), bool)
    P1 = pts[inW]
    SX = (P1 @ S.T) % 3
    w = np.einsum('pi,ij,pj->p', P1, Om, SX) % 3
    eig = ((SX - P1) % 3 == 0).all(1) | ((SX + P1) % 3 == 0).all(1)
    sel = np.flatnonzero((w != 0) & ~eig)
    U, V = P1[sel], SX[sel]; S2 = (V @ S.T) % 3
    inv_ = np.zeros(len(U), bool)
    for a in range(3):
        for b in range(3): inv_ |= (((a * U + b * V - S2) % 3) == 0).all(1)
    cand = {}
    for i in np.flatnonzero(inv_):
        Pp = np.stack([U[i], V[i]], 1)
        key = frozenset(tuple((a * Pp[:, 0] + b * Pp[:, 1]) % 3) for a in range(3) for b in range(3))
        cand.setdefault(key, Pp)
    extra = T.max_orthogonal_family(list(cand.values()), n)
    fam += extra
    c = len(fam)
    if c == n - 1: c = n
    need = n - len(fam)
    if need <= 0 or c == n:
        planes = fam if len(fam) == n else fam + [T.nullspace((np.concatenate(fam, 1).T @ Om) % 3).T]
        return c, T.check_split(S, planes, Om)
    Pm = np.concatenate(fam, 1) if fam else np.zeros((D, 0), np.int64)
    inW = ~(((pts @ Om @ Pm) % 3) != 0).any(1) if fam else np.ones(len(pts), bool)
    P2 = pts[inW]
    SX = (P2 @ S.T) % 3
    w = np.einsum('pi,ij,pj->p', P2, Om, SX) % 3
    eig = ((SX - P2) % 3 == 0).all(1) | ((SX + P2) % 3 == 0).all(1)
    good = np.flatnonzero((w != 0) & ~eig)
    X = np.stack([P2[good], SX[good]], 2); XJ = np.einsum('mai,ab->mib', X, Om)
    rng = np.random.default_rng(0)
    def orth(i, cc):
        if len(cc) == 0: return cc
        Wm = np.einsum('ib,pbj->pij', XJ[i], X[cc]) % 3
        return cc[~Wm.reshape(len(cc), 4).any(1)]
    def complement_ok(chosen):
        planes = fam + [X[i] for i in chosen]
        comp = T.nullspace((np.concatenate(planes, 1).T @ Om) % 3).T
        if comp.shape[1] != 2 or not ((comp.T @ Om @ comp) % 3).any(): return None
        if T.intersect_dim(comp, (S @ comp) % 3) >= 1: return planes + [comp]
        return None
    def dfs(cc, chosen, bud):
        if len(chosen) == need - 1:
            return complement_ok(chosen)
        for t, i in enumerate(cc):
            bud[0] -= 1
            if bud[0] < 0: return None
            r = dfs(orth(i, cc[t + 1:]), chosen + [int(i)], bud)
            if r: return r
        return None
    for _ in range(300):
        r = dfs(rng.permutation(len(X)), [], [20000])
        if r:
            return c, T.check_split(S, r, Om)
    return c, None


GAP_N5 = ROOT / "data" / "w33_pass11191_gap_class_reps_n5.txt"
PSP10_3 = 3 ** 25 * (3 ** 2 - 1) * (3 ** 4 - 1) * (3 ** 6 - 1) * (3 ** 8 - 1) * (3 ** 10 - 1) // 2


def qutrit_n5_classes():
    """for every class: c found constructively (eigenspace parts + invariant span(x,Sx) planes) and a split with
    E = 5 - c (which meets the lower bound, so A = 5 - c exactly); compared with the Jordan closed form"""
    import w33_pass11191_arrow_five_six_qutrits as P91
    text = re.sub(r"\s+", " ", GAP_N5.read_text().replace("\\\n", ""))
    cls = []
    for order, size, mat in re.findall(r"CLASS n=5 order=(\d+) size=(\d+) fixed=-1 real=\w+ M=\s*(\[ \[.*?\] \])", text):
        rows = re.findall(r"\[([^\[\]]*)\]", mat)
        cls.append((int(order), int(size), np.array([[int(x) for x in r.split(",")] for r in rows], np.int64).T % 3))
    pts = P91.points(10)
    dist, bad = Counter(), []
    for k, (o, sz, S) in enumerate(cls):
        c, E = exact_A(S, 5, pts)
        cj = T.c_jordan(S, 5)
        if not (E is not None and E == 5 - c and c == cj):
            bad.append(dict(index=k, order=o, c=c, c_jordan=cj, E=E))
        dist[5 - cj] += sz
    tot = sum(dist.values())
    return dict(classes=len(cls), total_ok=tot == PSP10_3, exact_everywhere=not bad, failures=bad[:20],
                A_distribution={str(k): v for k, v in sorted(dist.items())},
                A_fractions={str(k): str(Fraction(v, tot)) for k, v in sorted(dist.items())})


def large_n_statistics(seed=2026):
    import w33_pass11165_three_qutrit_arrow as A3
    import w33_pass11199_arrow_theorem as T
    rng = np.random.default_rng(seed)
    out = {}
    for n, N in ((6, 4000), (8, 4000), (12, 3000), (16, 2000), (24, 1000), (32, 600)):
        cnt = Counter()
        for _ in range(N):
            S, _J = A3.random_symplectic(rng, n=n, steps=30 * n)
            cnt[T.c_jordan(S.astype(np.int64) % 3, n)] += 1
        ec = sum(j * v for j, v in cnt.items()) / N
        out[f"n{n}"] = dict(sampled=N, c_distribution={str(j): v for j, v in sorted(cnt.items())},
                            P_A_equals_n=cnt[0] / N, E_c=ec, E_A_over_n=1 - ec / n)
    return out


def exact_small_qutrits():
    d = json.loads((ROOT / "data" / "w33_pass11199_arrow_theorem.json").read_text())
    out = {}
    for n in ('n2', 'n3', 'n4'):
        dist = {int(k): v for k, v in d['classes'][n]['c_distribution_psp'].items()}
        tot = sum(dist.values())
        out[n] = dict(c_distribution={str(j): str(Fraction(v, tot)) for j, v in sorted(dist.items())},
                      E_c=str(Fraction(sum(j * v for j, v in dist.items()), tot)))
    return out


def summarize(include_big=True):
    t0 = time.time()
    res = dict(pass_id=11200, classes=[])
    cases = [(2, 2, GAP_SMALL), (2, 3, GAP_SMALL), (2, 4, GAP_SMALL), (5, 2, GAP_SMALL)]
    if include_big:
        cases += [(7, 2, GAP_BIG), (2, 5, GAP_BIG)]
    for q, n, path in cases:
        res['classes'].append(class_check(q, n, path))
        OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    # qutrit Burnside check from Pass 11188's class data
    import w33_pass11188_exact_arrow_by_class as P88
    D = json.loads(P88.OUT.read_text())
    res['qutrit_burnside'] = {n: str(Fraction(sum(r['size'] * r['invariant_planes'] for r in D[n]['rows']),
                                              sum(r['size'] for r in D[n]['rows']))) for n in ('n2', 'n3', 'n4')}
    res['qutrit_exact_small'] = exact_small_qutrits()
    res['qutrit_n5_classes'] = qutrit_n5_classes()
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    res['qutrit_large_n'] = large_n_statistics()
    res['law_holds_everywhere'] = (all(c['law_holds'] and c['total_ok'] for c in res['classes'])
                                   and res['qutrit_n5_classes']['exact_everywhere'] and res['qutrit_n5_classes']['total_ok'])
    res['seconds'] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    for c in r['classes']:
        print({k: c[k] for k in ('q', 'n', 'classes', 'total_ok', 'law_holds', 'A_max', 'A_equals_n_fraction',
                                 'burnside_mean_invariant_planes')})
    print('qutrit burnside', r['qutrit_burnside'])
    print('qutrit exact', r['qutrit_exact_small'])
    print('qutrit n=5 classes', {k: v for k, v in r['qutrit_n5_classes'].items()})
    for k, v in r['qutrit_large_n'].items():
        print(k, v)
