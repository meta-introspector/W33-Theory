#!/usr/bin/env python3
"""Pass 11182: which subsystems the paper's own ticks select -- quantum mereology for the clock and the perfect ticks.

For each tick S of three qutrits and each of the 110565 tensor factorisations F (Pass 11180 method, S_F = B_F^-1 S B_F):
local (S maps the planes of F onto planes of F), perfect (all blocks of S_F invertible) or partial.
Results:
  * the clock K of Theorem 4.3 (order 3), its square and its Fourier-dual kick V: local in 108 splits, perfect in 8748,
    partial in the remaining 101709.  In all 108 the clock fixes each plane -- it is literally a product of three
    independent single-qutrit gates there.  None of the 108 contains a light-cone coordinate plane and only 3 contain
    the transverse one; exactly 4 keep positions and momenta separate, and those 4 are the 4 orthogonal frames of the
    light-cone metric q (the only q-orthogonal frames of F_3^3).  So the clock's natural subsystems are the principal
    axes of the Lorentzian form (4 of them) plus 104 splits that mix positions with momenta -- not the light-cone
    coordinates the paper writes it in;
  * the perfect interacting tick V(N) K V(N), N = (x0 + x1 + x2)^2 (order 9): local in NO split -- intrinsically
    entangling (its class is fixed-point-free, Pass 11181); perfect in 4860;
  * the F9 gate K = P^2 + i*1 (order 3): local in 27 splits, perfect in 3564;
  * random ticks (300): 2/3 local in no split (exact value 7922/12285, Pass 11181).
Reading: the clock does not entangle anything -- in the right subsystems; the interaction that makes a tick scramble
perfectly also removes every subsystem split in which it could be seen as local.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11168_clock_dual_pair as D  # noqa: E402
import w33_pass11180_mereology as M  # noqa: E402
import w33_pass11165_three_qutrit_arrow as A3  # noqa: E402

OUT = ROOT / "data" / "w33_pass11182_paper_ticks_mereology.json"
CACHE = Path(r"C:/Users/wiljd/AppData/Local/Temp/facs3.npy")


def load_facs():
    if CACHE.exists():
        B = np.load(CACHE)
        if len(B) == 110565:
            return B
    B = np.array(M.factorisations(3))
    try:
        np.save(CACHE, B)
    except OSError:
        pass
    return B


def ticks():
    K = D.xz_to_party(D.I3, D.MQ, D.Z3, D.I3)
    V = D.kick(D.MQ)
    J3 = np.ones((3, 3), np.int64)
    VKV = (D.kick(J3) @ K @ D.kick(J3)) % 3
    F9 = A3.tetracode3()
    return dict(clock_K=K, clock_K2=(K @ K) % 3, dual_kick_V=V, perfect_tick_VKV=VKV, f9_gate=F9)


def clock_splits(Bs):
    """characterise the splits in which the clock is local"""
    import itertools
    K = ticks()['clock_K']
    J = M.form(3)
    Binv = np.einsum('ij,fkj,kl->fil', -J, Bs, J) % 3
    SF = np.einsum('fij,jk,fkl->fil', Binv, K, Bs) % 3
    Bk = SF.reshape(-1, 3, 2, 3, 2).transpose(0, 1, 3, 2, 4)
    zero = (Bk == 0).all(axis=(-1, -2))
    idx = np.flatnonzero(((~zero).sum(-1) == 1).all(-1))
    all_planes_fixed = int(sum(all(not zero[f, k, k] for k in range(3)) for f in idx))

    def key(u, v):
        return frozenset(tuple((a * u + b * v) % 3) for a in range(3) for b in range(3))
    E = np.eye(6, dtype=np.int64)
    coord = [key(E[0], E[1]), key(E[2], E[3]), key(E[4], E[5])]
    with_transverse = with_lightcone = aligned = orth = 0
    for f in idx:
        B = Bs[f]
        planes = [key(B[:, 2 * k], B[:, 2 * k + 1]) for k in range(3)]
        with_transverse += coord[1] in planes
        with_lightcone += (coord[0] in planes) or (coord[2] in planes)
        xs, ok = [], True
        for k in range(3):
            sp = [(a * B[:, 2 * k] + b * B[:, 2 * k + 1]) % 3 for a in range(3) for b in range(3)]
            xo = [v for v in sp if v.any() and not v[[1, 3, 5]].any()]
            zo = [v for v in sp if v.any() and not v[[0, 2, 4]].any()]
            if not xo or not zo:
                ok = False
                break
            xs.append(xo[0][[0, 2, 4]])
        if ok:
            aligned += 1
            orth += all((xs[a] @ D.MQ @ xs[b]) % 3 == 0 for a in range(3) for b in range(3) if a != b)
    V = [np.array(v) for v in itertools.product(range(3), repeat=3) if any(v)]
    lines = {tuple((v * pow(int(v[next(k for k in range(3) if v[k])]), -1, 3)) % 3) for v in V}
    nonnull = [np.array(l) for l in lines if (np.array(l) @ D.MQ @ np.array(l)) % 3]
    frames = sum(1 for c in itertools.combinations(nonnull, 3)
                 if all((a @ D.MQ @ b) % 3 == 0 for a, b in itertools.combinations(c, 2)))
    return dict(local_splits=len(idx), all_planes_fixed=all_planes_fixed, containing_transverse=with_transverse,
                containing_lightcone=with_lightcone, position_momentum_aligned=aligned, aligned_q_orthogonal=orth,
                q_orthogonal_frames=frames)


def summarize(n_random=300):
    Bs = load_facs()
    out = {}
    for name, S in ticks().items():
        p = M.profile(S, Bs, 3)
        p['order'] = M.porder(S % 3, 3)
        out[name] = p
    rng = np.random.default_rng(11182)
    rand = defaultdict(Counter)
    free = 0
    for _ in range(n_random):
        S, _J = A3.random_symplectic(rng)
        p = M.profile(S, Bs, 3)
        o = M.porder(S % 3, 3)
        rand[o][(p['local'], p['perfect'])] += 1
        free += p['local'] == 0
    res = dict(pass_id=11182, factorisations=len(Bs), ticks=out, clock_splits=clock_splits(Bs),
               random_sample=n_random, random_local_none_fraction=free / n_random,
               random_by_order={str(o): {f"{k[0]}|{k[1]}": v for k, v in sorted(c.items())} for o, c in sorted(rand.items())})
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    for k, v in r['ticks'].items():
        print(k, v)
    print('random: local in none', r['random_local_none_fraction'])
    for o, c in r['random_by_order'].items():
        print(o, c)
