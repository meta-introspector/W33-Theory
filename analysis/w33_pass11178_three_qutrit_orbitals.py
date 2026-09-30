#!/usr/bin/env python3
"""Pass 11178: the relations between three-qutrit tensor factorisations -- a rank-20 action of Sp(6,3) -- and how much of
it a gate's Choi-entanglement signature sees.

GAP (analysis/gap/w33_pass11178_orbital_signatures.g; frozen output data/w33_pass11178_gap_orbital_signatures.txt):
Sp(6,3) acts on the 110565 factorisations F3^6 = P1 + P2 + P3 (orthogonal nondegenerate planes) with rank 20.  For a
representative F' = (Q1,Q2,Q3) of each suborbit, block (i,j) is the projection of Q_j onto P_i along the other two
planes, in symplectic bases: 0 zero, 1 rank one, 2 invertible with det +1, 3 invertible with det -1.  For a gate S with
S(F0) = F', these are exactly its party blocks S_ij, so the code matrix up to relabelling rows and columns (S3 x S3) is
the gate's Choi-entanglement signature (which input-output pairs are entangled, fully or by one trit, and with which
orientation).
Results:
  * the 20 orbitals give 15 signature classes: 13 orbitals are determined by their signature; one class merges two
    orbitals (6912 + 6912), and one merges five (256 + 256 + 2304 + 6912 + 6912 = 16640) -- the class 'one
    orientation-preserving invertible block in each row and column, every other block rank one';
  * the perfect relation (every block invertible, one orientation reversal per row and column) is the single orbital of
    size 3456 = 110565 * 128/4095 (Pass 11163);
  * the orbital sizes over |Sp(6,3)| reproduce the frequencies of a direct sample of Sp(6,3) (recomputed here).
Two qutrits (Pass 11177): rank 3, and the signature (ranks of S_AA, S_BA) is complete.  Three qutrits: rank 20, and the
signature misses 5 of the 20 relations -- the finer invariant separating them is open.
"""
from __future__ import annotations

import itertools
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GAPOUT = ROOT / "data" / "w33_pass11178_gap_orbital_signatures.txt"
OUT = ROOT / "data" / "w33_pass11178_three_qutrit_orbitals.json"
PERMS = list(itertools.permutations(range(3)))


def canon(M):
    M = np.array(M)
    return min(tuple(M[list(r)][:, list(c)].flatten()) for r in PERMS for c in PERMS)


def gap_classes():
    t = GAPOUT.read_text().replace('\n', ' ')
    rows = re.findall(r'SUB (\d+) (\[.*?\]\s*\])', t)
    g = defaultdict(list)
    for size, m in rows:
        g[canon(eval(m.replace(' ', '')))].append(int(size))
    return len(rows), g


def sample_classes(n=300000, seed=11178):
    D = 6
    J = np.zeros((D, D), np.int64)
    for k in range(3):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
    rng = np.random.default_rng(seed)
    cnt = Counter()
    done = 0
    while done < n:
        m = 50000
        M = np.broadcast_to(np.eye(D, dtype=np.int64), (m, D, D)).copy()
        for _ in range(80):
            v = rng.integers(0, 3, (m, D))
            c = rng.integers(1, 3, (m, 1))
            w = np.einsum('mi,ij,mjk->mk', v, J, M) % 3
            M = (M + c[:, :, None] * v[:, :, None] * w[:, None, :]) % 3
        B = M.reshape(-1, 3, 2, 3, 2).transpose(0, 1, 3, 2, 4)
        det = (B[..., 0, 0] * B[..., 1, 1] - B[..., 0, 1] * B[..., 1, 0]) % 3
        zero = (B == 0).all(axis=(-1, -2))
        code = np.where(zero, 0, np.where(det == 0, 1, np.where(det == 1, 2, 3)))
        best = None
        for r in PERMS:
            for c in PERMS:
                k = code[:, list(r)][:, :, list(c)].reshape(-1, 9)
                key = (k * (4 ** np.arange(8, -1, -1))).sum(1)
                best = key if best is None else np.minimum(best, key)
        for key, v in Counter(best.tolist()).items():
            cnt[tuple((key // 4 ** (8 - i)) % 4 for i in range(9))] += v
        done += m
    return cnt, done


def summarize():
    n_orb, g = gap_classes()
    cnt, n = sample_classes()
    table = []
    maxdev = 0.0
    for k, sizes in sorted(g.items(), key=lambda kv: sum(kv[1])):
        exp = sum(sizes) / 110565
        obs = cnt.get(k, 0) / n
        sd = np.sqrt(exp * (1 - exp) / n)
        maxdev = max(maxdev, abs(obs - exp) / max(sd, 1e-12))
        table.append(dict(signature=np.array(k).reshape(3, 3).tolist(), orbitals=sorted(sizes), total=sum(sizes),
                          sampled_fraction=obs, expected_fraction=exp))
    perfect = [r for r in table if all(x in (2, 3) for row in r['signature'] for x in row)]
    res = dict(pass_id=11178, orbitals=n_orb, signature_classes=len(g), total=sum(sum(v) for v in g.values()),
               merged=[r['orbitals'] for r in table if len(r['orbitals']) > 1], perfect=perfect,
               sample_size=n, max_sigma_deviation=maxdev, unseen_in_sample=len(set(g) - set(cnt)),
               extra_in_sample=len(set(cnt) - set(g)), table=table)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    r = summarize()
    print(json.dumps({k: v for k, v in r.items() if k != 'table'}, indent=1))
