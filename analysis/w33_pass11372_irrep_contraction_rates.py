"""Pass 11372: contraction rates of the Clifford+T walk irrep by irrep, and their supremum.

The one-qutrit walk U_k = C_1 T C_2 T ... C_k T (C_i uniform in the Clifford group) has Fourier transform
mu^(lambda) = P_lambda lambda(T) on each irrep lambda of PU(3) (P_lambda = Clifford average).  Its spectrum lies in
that of P lambda(T) P, and rho_lambda = largest modulus < 1 governs how fast E f(U_k) converges for f in lambda.
Passes 11354/11359 computed rho on TENSOR spaces (mixtures of irreps).  Here rho is resolved irrep by irrep:

  * Sym^p V (x) Sym^q Vbar = (p,q) + (p-1,q-1) + ... + (p-m,q-m),  m = min(p,q)  (Dynkin labels),
    so the moduli that are new at (p,q) belong to the irrep (p,q);  p = q mod 3 (trivial centre);
  * matrix-free: Sym^p(g) by expanding products of linear forms (orthonormal monomial basis); the Pauli average is exact
    combinatorics (Z-weight 0, X-orbit sums); the remaining average is over 24 coset representatives (Pass 11312).

THEOREM (by citation).  sup_lambda rho_lambda < 1.  Benoist-de Saxce (Invent. Math. 2016): a probability measure with
algebraic support generating a dense subgroup of a compact simple Lie group has a spectral gap.  Applied to
nu = mu*mu*mu^vee*mu^vee (support C T C' T T^-1 C''^-1 T^-1 C'''^-1, containing Cl and T Cl T^-1), whose support
generates a dense subgroup -- checked here: the Lie algebra spanned by logarithms of infinite-order elements of
<Cl, T Cl T^-1> is all of su(3) -- one gets ||nu^(lambda)|| <= 1 - delta for every nontrivial lambda, hence
rho_lambda <= (1 - delta)^(1/4) uniformly.  The value of the supremum is not determined by that argument.

Findings: rho is not monotone in the degree; the largest value found is rho_(18,0) = largest root of
864 x^3 - 984 x^2 + 230 x + 9 = 0.7809... (p + q <= 24, and (p,0) to p = 48 by coherent states, --family).
rho_(8,8) = 0.7228 is close to the empirical ~0.72-per-gate decay of the reversible fraction (Pass 11312), but the sup
exceeds it, and the reversible set is Haar-null (its indicator is not an L^2 observable), so no Fourier rate theorem
ties the two: the match is recorded as unexplained.  Exact moment law (--moments): E|tr U_k|^6 = 6 + (3/8)^k,
E|tr U_k|^8 = 23 + 16 (3/8)^k, checked by enumeration.
"""

from __future__ import annotations

import itertools
import json
import sys
from math import factorial
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11213_cubic_t_violation as P1  # noqa: E402
import w33_pass11355_substrate_jarlskog as SJ  # noqa: E402

OUT = ROOT / "data" / "w33_pass11372_irrep_contraction_rates.json"
ZETA = np.exp(2j * np.pi / 9)


def monomials(p):
    return [a for a in itertools.product(range(p + 1), repeat=3) if sum(a) == p]


def sym_power(A, p, mons):
    """matrix of g = A on Sym^p V in the orthonormal basis u_a = sqrt(p!/a!) x^a (x^a -> prod_i (sum_j A_ji x_j)^a_i)"""
    idx = {m: i for i, m in enumerate(mons)}
    norm = np.array([np.sqrt(factorial(p) / np.prod([factorial(x) for x in m])) for m in mons])
    out = np.zeros((len(mons), len(mons)), complex)
    for col, a in enumerate(mons):
        poly = {(0, 0, 0): 1.0 + 0j}
        for i in range(3):
            for _ in range(a[i]):
                new = {}
                for m, c in poly.items():
                    for j in range(3):
                        if A[j, i] != 0:
                            mm = list(m)
                            mm[j] += 1
                            mm = tuple(mm)
                            new[mm] = new.get(mm, 0) + c * A[j, i]
                poly = new
        for m, c in poly.items():
            out[idx[m], col] += c
    # x^a is the symmetrisation average of e^(x)a, of norm sqrt(a!/p!) = 1/norm_a; orthonormal u_a = norm_a x^a
    return (out / norm[:, None]) * norm[None, :]


_SYM = {}
_MONS = {}


def sym_cached(i, g, p, conj):
    key = (i, p, conj)
    if key not in _SYM:
        if p not in _MONS:
            _MONS[p] = monomials(p)
        _SYM[key] = sym_power(g.conj() if conj else g, p, _MONS[p])
    return _SYM[key]


class Space:
    def __init__(self, p, q, reps):
        self.p, self.q = p, q
        self.mp, self.mq = monomials(p), monomials(q)
        self.np_, self.nq = len(self.mp), len(self.mq)
        # Pauli-invariant subspace: Z-weight 0, then orbit sums under the cyclic shift X (x_i -> x_{i+1})
        w = lambda m: sum(i * m[i] for i in range(3))
        cells = [(i, j) for i, a in enumerate(self.mp) for j, b in enumerate(self.mq) if (w(a) - w(b)) % 3 == 0]
        ip = {m: i for i, m in enumerate(self.mp)}
        iq = {m: i for i, m in enumerate(self.mq)}
        shift = lambda m: (m[2], m[0], m[1])
        seen = set()
        basis = []
        for i, j in cells:
            if (i, j) in seen:
                continue
            orb = []
            a, b = self.mp[i], self.mq[j]
            for _ in range(3):
                orb.append((ip[a], iq[b]))
                a, b = shift(a), shift(b)
            orb = sorted(set(orb))
            seen.update(orb)
            v = np.zeros((self.np_, self.nq), complex)
            for x, y in orb:
                v[x, y] = 1
            basis.append(v / np.linalg.norm(v))
        # the X phase convention: shift is a permutation of the ON basis (all norms equal along an orbit)
        self.B0 = np.array(basis)                                   # (m, np, nq)
        self.reps = [(sym_cached(i, g, p, False), sym_cached(i, g, q, True)) for i, g in enumerate(reps)]
        self.Tp = np.array([ZETA ** (sum((i ** 3) * m[i] for i in range(3)) % 9) for m in self.mp])
        self.Tq = np.array([ZETA ** (-(sum((i ** 3) * m[i] for i in range(3))) % 9) for m in self.mq])

    def clifford_invariants(self):
        """orthonormal basis (as (k, np, nq) tensors) of the Clifford-invariant subspace"""
        m = len(self.B0)
        Pm = np.zeros((m, m), complex)
        for A, Bc in self.reps:
            img = np.einsum('ij,mjk,lk->mil', A, self.B0, Bc)       # (A (x) B) applied to every basis vector
            Pm += np.einsum('nil,mil->nm', self.B0.conj(), img)
        Pm /= len(self.reps)
        w, V = np.linalg.eigh((Pm + Pm.conj().T) / 2)
        keep = V[:, w > 0.5]
        return np.einsum('mk,mil->kil', keep, self.B0)

    def moduli(self):
        Binv = self.clifford_invariants()
        if len(Binv) == 0:
            return np.array([])
        TB = Binv * self.Tp[None, :, None] * self.Tq[None, None, :]
        Mt = np.einsum('kil,jil->kj', Binv.conj(), TB)
        return np.sort(np.abs(np.linalg.eigvals(Mt)))[::-1]


def coset_reps():
    """one Clifford per symplectic class (representatives of Cl / Pauli): classified by C X C^dag, C Z C^dag"""
    w = np.exp(2j * np.pi / 3)
    X = np.roll(np.eye(3), 1, axis=0)
    Z = np.diag([1, w, w * w])
    W = {(a, b): np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b) for a in range(3) for b in range(3)}

    def label(P):
        for k, V in W.items():
            if abs(abs(np.trace(V.conj().T @ P)) - 3) < 1e-9:
                return k
        raise AssertionError

    out = {}
    for C in SJ.CF:
        key = (label(C @ X @ C.conj().T), label(C @ Z @ C.conj().T))
        out.setdefault(key, C / np.linalg.det(C) ** (1 / 3))       # det 1: phases are cube roots, trivial for p = q mod 3
    return list(out.values())


def new_moduli(cur, prev, tol=1e-8):
    cur = list(cur)
    for x in prev:
        k = int(np.argmin([abs(c - x) for c in cur]))
        assert abs(cur[k] - x) < 1e-6, "lower irreps must reappear"
        cur.pop(k)
    return cur


def lie_algebra_dimension():
    T = P1.T1
    CF = SJ.CF
    rng = np.random.default_rng(0)
    gens = []
    for _ in range(40):
        g = T @ CF[rng.integers(216)] @ T.conj().T @ CF[rng.integers(216)]
        w, V = np.linalg.eig(g / np.linalg.det(g) ** (1 / 3))
        ph = np.angle(w)
        H = (V @ np.diag(ph) @ np.linalg.inv(V))                    # log(g)/i, Hermitian up to rounding
        H = (H + H.conj().T) / 2
        H -= np.trace(H) / 3 * np.eye(3)
        gens.append(H)
    span = []
    def add(X):
        v = np.concatenate([X.real.ravel(), X.imag.ravel()])
        M = np.array(span + [v])
        if np.linalg.matrix_rank(M, tol=1e-8) > len(span):
            span.append(v)
            return True
        return False
    for H in gens:
        add(H)
    return len(span)


def run(maxdeg=24):
    res = dict(pass_id=11372)
    res["lie_algebra_dim_of_closure_of_<Cl,T Cl T^-1>"] = lie_algebra_dimension()
    reps = coset_reps()
    assert len(reps) == 24
    spec = {}
    table = {}
    for s in range(0, maxdeg + 1):
        for p in range(s + 1):
            q = s - p
            if (p - q) % 3:
                continue
            sp = Space(p, q, reps)
            mods = sp.moduli()
            spec[(p, q)] = list(mods)
            prev = spec.get((p - 1, q - 1), []) if min(p, q) > 0 else []
            nm = new_moduli(mods, prev)
            nm = sorted([x for x in nm if x < 1 - 1e-9], reverse=True)
            dim = (p + 1) * (q + 1) * (p + q + 2) // 2
            table[f"({p},{q})"] = dict(dim=dim, clifford_invariants=len(new_moduli(mods, prev)) if mods.size else 0,
                                       rho=float(nm[0]) if nm else 0.0,
                                       top=[float(f"{x:.10f}") for x in nm[:4]])
            print(p, q, table[f"({p},{q})"], flush=True)
    res["irreps"] = table
    rhos = [(v["rho"], k) for k, v in table.items()]
    rhos.sort(reverse=True)
    res["largest_rho"] = rhos[:8]
    by_deg = {}
    for k, v in table.items():
        p, q = map(int, k.strip("()").split(","))
        by_deg[p + q] = max(by_deg.get(p + q, 0.0), v["rho"])
    res["max_rho_by_p_plus_q"] = dict(sorted(by_deg.items()))
    return res


def trace_moment_law(ts=(3, 4), kmax=3):
    """EXACT LAW.  The steps C_i T are independent, so E pi(U_k) = (P pi(T))^k on any tensor representation and
        E|tr U_k|^{2t} = tr((P_t pi_t(T))^k) = sum of nu^k over the spectrum of P_t pi_t(T)
    (cyclicity: = tr((B^dag pi(T) B)^k), B an orthonormal basis of Clifford invariants in V^(x)t (x) Vbar^(x)t).
    The spectrum is computed (Pass 11359's matrix-free basis) and the law is checked by enumerating every
    coset-reduced word (Pass 11312) for k = 1..kmax."""
    import w33_pass11312_depth_law as DL
    import w33_pass11359_higher_degree_rates as H
    rng = np.random.default_rng(1)
    Cs = [H.det1(C) for C in P1.clifford1()]
    td = np.diag(H.det1(P1.T1))
    out = {}
    spec = {}
    for t in ts:
        B, _ = H.invariant_basis(t, t, Cs, rng)
        ph = np.ones((3,) * (2 * t), complex)
        for leg in range(2 * t):
            v = td if leg < t else td.conj()
            sh = [1] * (2 * t)
            sh[leg] = 3
            ph = ph * v.reshape(sh)
        ev = np.linalg.eigvals(B.conj().T @ (ph.reshape(-1, 1) * B))
        ev = ev[np.abs(ev) > 1e-10]
        spec[t] = ev
        vals = {}
        for z in ev:
            key = f"{z.real:.10f}" if abs(z.imag) < 1e-9 else f"{z.real:.10f}{z.imag:+.10f}i"
            vals[key] = vals.get(key, 0) + 1
        out[f"t={t}"] = dict(clifford_invariants=int(B.shape[1]), nonzero_spectrum_multiplicities=vals)
    CF = np.array(P1.clifford1())
    Rr = np.array(DL.coset_reps(list(CF)))
    cur = np.einsum('cij,jk->cik', CF, P1.T1)
    checks = []
    for k in range(1, kmax + 1):
        tr = np.abs(np.einsum('cii->c', cur))
        for t in ts:
            enum, law = float(np.mean(tr ** (2 * t))), float(np.sum(spec[t] ** k).real)
            checks.append(dict(t=t, k=k, enumerated=enum, law=law, dev=abs(enum - law)))
        if k < kmax:
            cur = np.einsum('rij,jk,ckl->rcil', Rr, P1.T1, cur).reshape(-1, 3, 3)
    out["checks"] = checks
    out["max_dev"] = max(c["dev"] for c in checks)
    out["closed_forms"] = {"t=3": "E|tr U_k|^6 = 6 + (3/8)^k", "t=4": "E|tr U_k|^8 = 23 + 16 (3/8)^k"}
    return out


def coherent_spectrum(p, rng):
    """moduli of P pi(T) on Sym^p V by coherent states: w_j = v_j^(x)p span Sym^p, <w_j, w_l> = (v_j^dag v_l)^p and
    <w_j, P pi(T) w_l> = avg_C (v_j^dag C T v_l)^p, so the operator has matrix G^-1 K in the basis w_j"""
    CF = np.array([C / np.linalg.det(C) ** (1 / 3) for C in SJ.CF])
    CT = np.einsum('cij,jk->cik', CF, P1.T1 / np.linalg.det(P1.T1) ** (1 / 3))
    N = (p + 1) * (p + 2) // 2
    v = rng.normal(size=(N, 3)) + 1j * rng.normal(size=(N, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    G = (v.conj() @ v.T) ** p
    K = np.mean(np.einsum('ji,cli->cjl', v.conj(), np.einsum('cij,lj->cli', CT, v)) ** p, axis=0)
    m = np.sort(np.abs(np.linalg.eigvals(np.linalg.solve(G, K))))[::-1]
    return [float(x) for x in m[m > 1e-8]], float(np.linalg.cond(G))


def coherent_pq(p, q, seed=3):
    """INDEPENDENT cross-check (no code shared with Space): moduli of P pi(T) on Sym^p V (x) Sym^q Vbar from coherent
    states w = v^(x)p (x) conj(u)^(x)q, <w_j, w_l> = (v_j^dag v_l)^p conj(u_j^dag u_l)^q"""
    rng = np.random.default_rng(seed)
    CT = np.einsum('cij,jk->cik', np.array(P1.clifford1()), P1.T1)
    N = ((p + 1) * (p + 2) // 2) * ((q + 1) * (q + 2) // 2)
    v = rng.normal(size=(N, 3)) + 1j * rng.normal(size=(N, 3))
    u = rng.normal(size=(N, 3)) + 1j * rng.normal(size=(N, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    u /= np.linalg.norm(u, axis=1)[:, None]
    G = (v.conj() @ v.T) ** p * np.conj(u.conj() @ u.T) ** q
    K = np.zeros((N, N), complex)
    for M in CT:
        K += (v.conj() @ (M @ v.T)) ** p * np.conj(u.conj() @ (M @ u.T)) ** q
    K /= len(CT)
    m = np.sort(np.abs(np.linalg.eigvals(np.linalg.solve(G, K))))[::-1]
    return [float(x) for x in m[m > 1e-7][:6]]


def charpoly(p, q, reps):
    sp = Space(p, q, reps)
    B = sp.clifford_invariants()
    Mt = np.einsum('kil,jil->kj', B.conj(), B * sp.Tp[None, :, None] * sp.Tq[None, None, :])
    c = np.poly(np.linalg.eigvals(Mt))
    assert np.abs(c.imag).max() < 1e-9
    return [float(x) for x in c.real]


def design_rates(tmax=10):
    """the exact convergence rate of E pi_t(U_k) = (P pi_t(T))^k to the Haar projector, i.e. of the walk to a unitary
    t-design: r_t = max rho over the irreps in V^(x)t (x) Vbar^(x)t (mixed weights, Pass 11370); needs the irrep table"""
    import w33_pass11370_twisted_indicator_degree_law as TW
    table = json.load(open(OUT))["irreps"]
    out = {}
    for t in range(1, tmax + 1):
        irr = {(w[0] - w[1], w[1] - w[2]) for _, _, _, w in TW.mixed_weights(t, 3)}
        irr = {x for x in irr if (x[0] - x[1]) % 3 == 0}
        assert all(f"({p},{q})" in table for p, q in irr), t
        best = max(irr, key=lambda x: table[f"({x[0]},{x[1]})"]["rho"])
        out[str(t)] = dict(rate=table[f"({best[0]},{best[1]})"]["rho"], attained_at=f"({best[0]},{best[1]})")
    return out


def family(pmax=48):
    """the (p,0) family to p = pmax (coherent states) and exact characteristic polynomials of the extreme blocks"""
    rng = np.random.default_rng(0)
    fam = {}
    for p in range(3, pmax + 1, 3):
        m, cond = coherent_spectrum(p, rng)
        fam[f"({p},0)"] = dict(moduli=[round(x, 10) for x in m], gram_condition=cond)
        print(p, fam[f"({p},0)"], flush=True)
    reps = coset_reps()
    cp = {}
    c = charpoly(18, 0, reps)                                        # Sym^18 V is the irrep (18,0) itself
    ints = [x * 864 for x in c]
    assert all(abs(x - round(x)) < 1e-6 for x in ints), ints
    cp["(18,0)"] = dict(integer_polynomial=[int(round(x)) for x in ints], largest_root=float(max(abs(np.roots(c)))))
    # (7,7): Sym^7 (x) Sym^7 also holds lower irreps; check that the two new eigenvalues solve 64x^2 - 16x - 13
    sp = Space(7, 7, reps)
    B = sp.clifford_invariants()
    ev = np.linalg.eigvals(np.einsum('kil,jil->kj', B.conj(), B * sp.Tp[None, :, None] * sp.Tq[None, None, :]))
    roots = sorted(float(z.real) for z in ev if abs(64 * z * z - 16 * z - 13) < 1e-8)
    cp["(7,7)"] = dict(integer_polynomial=[64, -16, -13], eigenvalues_solving_it=roots)
    assert len(roots) == 2
    # (8,8): the eigenvalues new at (8,8) (those of Sym^8 (x) Sym^8 minus those of Sym^7 (x) Sym^7)
    def block(p, q):
        sp = Space(p, q, reps)
        B = sp.clifford_invariants()
        return list(np.linalg.eigvals(np.einsum('kil,jil->kj', B.conj(), B * sp.Tp[None, :, None] * sp.Tq[None, None, :])))
    e8, e7 = block(8, 8), block(7, 7)
    for z in e7:
        k = int(np.argmin([abs(z - w) for w in e8]))
        assert abs(e8[k] - z) < 1e-7
        e8.pop(k)
    c = np.poly(e8).real * 10368
    assert np.abs(c - np.round(c)).max() < 1e-6
    cubic = np.polydiv(np.round(c), [4, -1])
    assert np.abs(cubic[1]).max() < 1e-6
    cp["(8,8)"] = dict(new_eigenvalue_polynomial=[int(x) for x in np.round(c)],
                       factorisation="(4x - 1)(2592x^3 - 1332x^2 - 585x + 140)",
                       cubic=[int(x) for x in np.round(cubic[0])], largest_root=float(max(abs(z) for z in e8)))
    cross = {f"({p},{q})": coherent_pq(p, q) for p, q in ((3, 3), (5, 5), (8, 8))}
    assert abs(cross["(8,8)"][1] - cp["(8,8)"]["largest_root"]) < 1e-6
    return dict(p0_family=fam, exact_characteristic_polynomials=cp, design_rates=design_rates(),
                independent_coherent_state_crosscheck=cross)


def main():
    if "--family" in sys.argv:
        res = json.load(open(OUT))
        res.update(family())
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    if "--moments" in sys.argv:
        res = json.load(open(OUT))
        res["trace_moment_law"] = trace_moment_law()
        print(json.dumps(res["trace_moment_law"], indent=1), flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "irreps"}, indent=1, default=str))


if __name__ == "__main__":
    main()
