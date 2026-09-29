"""Pass 11158 scan: perfect three-qutrit Clifford gates in Sp(6,3) by uniform sampling (200000 random products of 80 transvections); frozen as data/w33_pass11158_sp63_sample.json"""
import numpy as np, itertools, json, sys
rng = np.random.default_rng(11158)
n = 3; D = 2 * n
J = np.zeros((D, D), int)
for k in range(n): J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1, -1
def batch_random_symplectic(m, steps=80):
    M = np.tile(np.eye(D, dtype=int), (m, 1, 1))
    for _ in range(steps):
        v = rng.integers(0, 3, size=(m, D)); c = rng.integers(1, 3, size=m)
        # transvection x -> x + c * w(x, v) v applied on the left: M <- (I + c v (J v)^T ...) M ; w(x,v) = x^T J v
        Jv = (v @ J.T) % 3                      # J v per sample, shape (m, D): (J v)_i = sum_j J_ij v_j
        wv = np.einsum('mi,mij->mj', Jv, M) % 3   # (J v)^T M  -> w(M e_j, v) up to sign convention
        M = (M + c[:, None, None] * v[:, :, None] * wv[:, None, :]) % 3
    return M
def det2(A):  # (..., 2, 2)
    return (A[..., 0, 0] * A[..., 1, 1] - A[..., 0, 1] * A[..., 1, 0]) % 3
def det4(A):
    # batched determinant mod 3 via float det rounding (entries small)
    return np.round(np.linalg.det(A.astype(float))).astype(int) % 3
def blocks(M, R, C):
    rows = [2 * r + t for r in R for t in range(2)]; cols = [2 * c + t for c in C for t in range(2)]
    return M[:, rows][:, :, cols]
m = 200000
M = batch_random_symplectic(m)
ok_sym = np.all(((np.transpose(M, (0, 2, 1)) @ J @ M) % 3) == (J % 3), axis=(1, 2))
print('symplectic fraction', ok_sym.mean(), flush=True)
dets = np.stack([np.stack([det2(blocks(M, [i], [j])[:, :, :]) for j in range(n)], axis=1) for i in range(n)], axis=1)  # (m, i, j)
col_law = np.all(dets.sum(axis=1) % 3 == 1, axis=1); row_law = np.all(dets.sum(axis=2) % 3 == 1, axis=1)
single_ok = np.all(dets != 0, axis=(1, 2))
minors = []
for R in itertools.combinations(range(n), 2):
    for C in itertools.combinations(range(n), 2):
        minors.append(det4(blocks(M, list(R), list(C))) != 0)
double_ok = np.all(np.stack(minors, axis=1), axis=1)
perfect = single_ok & double_ok
pats = {}
for idx in np.nonzero(perfect)[0][:5000]:
    key = str(sorted(tuple(sorted(dets[idx][:, j].tolist())) for j in range(n)))
    pats[key] = pats.get(key, 0) + 1
out = dict(samples=m, symplectic_fraction=float(ok_sym.mean()), column_law=float(col_law.mean()), row_law=float(row_law.mean()),
           single_blocks_invertible=float(single_ok.mean()), perfect_fraction=float(perfect.mean()), perfect_count=int(perfect.sum()),
           perfect_column_det_patterns=pats, example=M[np.nonzero(perfect)[0][0]].tolist() if perfect.any() else None)
json.dump(out, open('sp63.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'example'}, indent=1))
