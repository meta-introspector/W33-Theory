"""Pass 11163: exact count of perfect three-qutrit Clifford gates. S <-> ordered symplectic basis (s1,s2 | s3,s4 | s5,s6);
the determinant of party block i of column pair (a,b) equals the local form w_i(a,b); the third column pair is fixed up to
SL(2,3) (24 choices) by the first two, and its block determinants are 1 - w_i(a,b) - w_i(c,d) (row law)."""
import numpy as np, itertools, json
V = np.array(list(itertools.product(range(3), repeat=6)), dtype=np.int64)       # 729 x 6
def local(u, v, i):
    return (u[..., 2 * i] * v[..., 2 * i + 1] - u[..., 2 * i + 1] * v[..., 2 * i]) % 3
Fi = np.stack([local(V[:, None, :], V[None, :, :], i) for i in range(3)], axis=-1)   # 729 x 729 x 3
W = Fi.sum(-1) % 3
good1 = (W == 1) & np.all(Fi != 0, axis=-1)
pairs1 = np.argwhere(good1)
print('symplectic pairs', int((W == 1).sum()), 'good first pairs', len(pairs1), flush=True)
total = 0
for (a, b) in pairs1:
    m = np.nonzero((W[:, a] == 0) & (W[:, b] == 0))[0]        # P1-perp (81 vectors)
    sub = Fi[np.ix_(m, m)]                                      # 81 x 81 x 3
    ok = (sub.sum(-1) % 3 == 1) & np.all(sub != 0, axis=-1)
    third = (1 - Fi[a, b][None, None, :] - sub) % 3
    ok &= np.all(third != 0, axis=-1)
    total += int(ok.sum())
count = 24 * total
SP63 = 3 ** 9 * 8 * 80 * 728
res = dict(symplectic_pairs=int((W == 1).sum()), good_first_pairs=len(pairs1), good_first_two=total, perfect_count=count,
           sp63=SP63, fraction=count / SP63)
json.dump(res, open('count63.json', 'w'), indent=1)
print(res)
