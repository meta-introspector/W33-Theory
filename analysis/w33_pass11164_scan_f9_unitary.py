"""Pass 11164: unitary matrices over F9 with all entries nonzero = F9-linear perfect gates (tetracode analogue)"""
import itertools, json
import numpy as np
EL = [(a, b) for a in range(3) for b in range(3)]                  # a + b i, i^2 = -1
def mul(x, y): return ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)
def add(x, y): return ((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)
def conj(x): return (x[0] % 3, (-x[1]) % 3)
def herm(u, v):
    s = (0, 0)
    for a, b in zip(u, v): s = add(s, mul(conj(a), b))
    return s
ZERO, ONE = (0, 0), (1, 0)
def unitary_group(n):
    vecs = list(itertools.product(EL, repeat=n))
    units = [v for v in vecs if herm(v, v) == ONE]
    out = []
    def extend(cols):
        if len(cols) == n:
            out.append(tuple(cols)); return
        for v in units:
            if all(herm(c, v) == ZERO for c in cols):
                extend(cols + [v])
    extend([])
    return out
def block(x):  # multiplication by x on F9 = F3^2 in basis (1, i): columns x*1, x*i
    a, b = x
    return np.array([[a, -b], [b, a]]) % 3
def embed(K):  # K given as columns
    n = len(K)
    S = np.zeros((2 * n, 2 * n), int)
    for j, col in enumerate(K):
        for i, x in enumerate(col):
            S[2 * i:2 * i + 2, 2 * j:2 * j + 2] = block(x)
    return S % 3
res = {}
for n in (2, 3):
    G = unitary_group(n)
    Jn = np.zeros((2 * n, 2 * n), int)
    for k in range(n): Jn[2 * k, 2 * k + 1], Jn[2 * k + 1, 2 * k] = 1, -1
    sym_ok = all(np.array_equal((embed(K).T @ Jn @ embed(K)) % 3, Jn % 3) for K in G[:200])
    nz = [K for K in G if all(x != ZERO for col in K for x in col)]
    res[n] = dict(order=len(G), embedding_symplectic=sym_ok, all_entries_nonzero=len(nz), example=[[list(x) for x in col] for col in nz[0]] if nz else None)
    print(n, res[n]['order'], sym_ok, len(nz), flush=True)
json.dump(res, open('f9.json', 'w'), indent=1)
