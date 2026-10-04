"""Pass 11369: completeness of the substrate Jarlskog witness on all of PU(3), not only on Clifford+T words.

Pass 11355: U is reversible iff U^T ~ C U C^dag (C Clifford, up to phase), and J_2t = 1/2 ||A_t(U) - A_t(U^T)||^2.

THEOREM A (the hierarchy is complete, with a finite degree).  Write P(U) = U (x) Ubar.  The projective Clifford group
(order 216) acts LINEARLY on P by conjugation, and P(U^T) = P(U)^T.  Polynomial invariants of a finite group separate
its orbits, and a separating set exists in degree <= |G| (Derksen-Kemper, Computational Invariant Theory, separating invariants;
Kemper 2009).  The linear invariant <Phi|P|Phi> = tr(U U^dag)/3 = 1 lifts every invariant of lower degree to degree t,
so A_t(U) contains all invariants of bidegree <= (t,t).  Hence
        U reversible  <=>  J_2t(U) = 0   for every t >= t*,   with 3 <= t* <= 216,
and completeness is monotone in t.  The question is whether t* = 3.

RESULT: t* > 3.  J_6 is NOT complete on PU(3) (it is complete on the Clifford+T words of Pass 11355):
  (G) global: gradient descent on J_6 from Haar-random starts reaches zeros; a zero at distance > 1e-3 from the
      reversible set where the NEXT witness J_8 is > 1e-6 is a CERTIFIED non-reversible zero (every J_2t vanishes on
      the reversible set).  One such point is refined to ||A_3(U) - A_3(U^T)|| ~ 1e-14 by Gauss-Newton on the full
      729 x 729 residual (scratch computation recorded in the note): dF has rank 4, a smooth 4-dimensional family.
      The same search on J_8 finds no certified spurious zero -- conjecture t* = 4.
  (L) local: the reversible set is a union of strata R_C = {U : U^T ~ C U C^dag}, one per twisted class C ~ conj(D) C D^-1;
      the J_6 Hessian rank is compared with codim(R_C) at a random point of each stratum.
  (C) control: the same pipeline on d = 5 with J_4, which Pass 11357 found INCOMPLETE on words, must find certified
      non-reversible zeros -- otherwise the search would be blind.
"""

from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11357_jarlskog_degree_by_dimension as P7  # noqa: E402

OUT = ROOT / "data" / "w33_pass11369_j6_completeness.json"
_S = {}


def hermitian_basis(d):
    out = []
    for i in range(d):
        for j in range(i + 1, d):
            A = np.zeros((d, d), complex)
            A[i, j] = A[j, i] = 1
            out.append(A)
            B = np.zeros((d, d), complex)
            B[i, j], B[j, i] = -1j, 1j
            out.append(B)
    for k in range(1, d):
        D = np.zeros((d, d), complex)
        D[:k, :k] = np.eye(k)
        D[k, k] = -k
        out.append(D / np.sqrt(k * (k + 1) / 2))
    return np.array(out)                                         # d^2 - 1 traceless generators


def expi(H):
    w, V = np.linalg.eigh(H)
    return (V * np.exp(1j * w)) @ V.conj().T


def haar(d, rng):
    Z = (rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def j_and_dirs(U, t, Cl, dUs):
    """J_2t(U) and its derivatives along the tangent vectors dUs (exact in the traces)"""
    Ud = U.conj().T
    X = np.einsum('cij,jk,clk->cil', Cl, U, Cl.conj())          # C U C^dag
    Y = np.einsum('cij,jk,clk->cil', Cl, U.T, Cl.conj())        # C U^T C^dag
    a = np.einsum('ij,cji->c', Ud, X)
    b = np.einsum('ij,cji->c', Ud, Y)
    J = float(np.mean(np.abs(a) ** (2 * t)) - np.mean(np.abs(b) ** (2 * t)))
    grad = np.zeros(len(dUs))
    for k, dU in enumerate(dUs):
        dUd = dU.conj().T
        dX = np.einsum('cij,jk,clk->cil', Cl, dU, Cl.conj())
        dY = np.einsum('cij,jk,clk->cil', Cl, dU.T, Cl.conj())
        da = np.einsum('ij,cji->c', dUd, X) + np.einsum('ij,cji->c', Ud, dX)
        db = np.einsum('ij,cji->c', dUd, Y) + np.einsum('ij,cji->c', Ud, dY)
        ga = t * np.abs(a) ** (2 * t - 2) * 2 * np.real(np.conj(a) * da)
        gb = t * np.abs(b) ** (2 * t - 2) * 2 * np.real(np.conj(b) * db)
        grad[k] = float(np.mean(ga) - np.mean(gb))
    return J, grad


def j_and_grad(U, t, Cl, G):
    """J_2t(U) and its gradient along U -> U exp(i sum th_a G_a) at th = 0"""
    return j_and_dirs(U, t, Cl, [U @ (1j * g) for g in G])


def rev_distance(U, Cl):
    """min over Clifford C and phase lambda of ||U^T - lambda C U C^dag||_F"""
    X = np.einsum('cij,jk,clk->cil', Cl, U, Cl.conj())
    ov = np.abs(np.einsum('cij,ij->c', X.conj(), U.T))
    d = U.shape[0]
    return float(np.sqrt(max(0.0, 2 * d - 2 * ov.max())))


def descend(U0, t, Cl, G, rounds=30):
    """BFGS in the chart th -> U_base exp(i th.G) (chain rule through a numerically differentiated 3x3 exponential),
    re-based until the step vanishes"""
    U = U0
    h = 1e-6
    for _ in range(rounds):
        U_base = U

        def chart(th):
            return U_base @ expi(np.einsum('a,aij->ij', th, G))

        def fg(th):
            dUs = []
            for k in range(len(G)):
                e = np.zeros(len(G))
                e[k] = h
                dUs.append((chart(th + e) - chart(th - e)) / (2 * h))
            return j_and_dirs(chart(th), t, Cl, dUs)

        r = minimize(fg, np.zeros(len(G)), jac=True, method='BFGS', options=dict(gtol=1e-15, maxiter=300))
        U = chart(r.x)
        if np.linalg.norm(r.x) < 1e-10:
            break
    return U, j_and_grad(U, t, Cl, G)[0]


def _init(d, t):
    _S["d"], _S["t"] = d, t
    _S["Cl"] = P7.clifford_group(d)
    _S["G"] = hermitian_basis(d)


def _job(seed):
    rng = np.random.default_rng(seed)
    d, t, Cl, G = _S["d"], _S["t"], _S["Cl"], _S["G"]
    U, J = descend(haar(d, rng), t, Cl, G)
    return dict(seed=seed, J=J, dist=rev_distance(U, Cl), J_next=P7.J(U, t + 1, Cl))


def global_search(d, t, n, seed0):
    with Pool(10, initializer=_init, initargs=(d, t)) as pool:
        rows = pool.map(_job, [seed0 + i for i in range(n)], chunksize=4)
    zeros = [r for r in rows if r["J"] < 1e-12]
    pos = [r for r in rows if r["J"] >= 1e-12]
    spurious = [r for r in zeros if r["dist"] > 1e-3]
    uncertified = [r for r in spurious if r["J_next"] <= 1e-6]
    return dict(runs=n, zeros=len(zeros), max_distance_at_zeros=max((r["dist"] for r in zeros), default=None),
                spurious_zeros=len(spurious), spurious_not_certified_by_next_witness=len(uncertified),
                min_next_witness_at_spurious=min((r["J_next"] for r in spurious), default=None),
                spurious_examples=sorted(spurious, key=lambda r: -r["dist"])[:5],
                positive_minima=len(pos),
                positive_minima_values=sorted({round(r["J"], 6) for r in pos})[:20])


# ------------------------------------------------------------------ local analysis of the strata
def index_of(M, Cl):
    M = P7.canon(M)
    return int(np.argmin(np.abs(Cl - M[None]).reshape(len(Cl), -1).max(axis=1)))


def twisted_classes(Cl):
    n = len(Cl)
    seen = -np.ones(n, dtype=int)
    reps = []
    for i in range(n):
        if seen[i] >= 0:
            continue
        orb = {index_of(D.conj() @ Cl[i] @ D.conj().T, Cl) for D in Cl}
        for j in orb:
            seen[j] = len(reps)
        reps.append((i, len(orb)))
    return reps


def point_on_stratum(C, G, rng, tries=30):
    """solve U^T = lambda C U C^dag by Gauss-Newton on U(d) (lambda fitted), or None"""
    d = C.shape[0]
    Gfull = list(G) + [np.eye(d, dtype=complex)]
    for _ in range(tries):
        U = haar(d, rng)
        for _ in range(200):
            X = C @ U @ C.conj().T
            lam = np.vdot(X, U.T)
            lam /= abs(lam)
            res = (U.T - lam * X).ravel()
            if np.linalg.norm(res) < 1e-13:
                return U
            Jm = []
            for g in Gfull:
                dU = U @ (1j * g)
                Jm.append((dU.T - lam * C @ dU @ C.conj().T).ravel())
            Jm.append((-1j * lam * X).ravel())
            Jm = np.array(Jm).T
            A = np.concatenate([Jm.real, Jm.imag])
            b = -np.concatenate([res.real, res.imag])
            step = np.linalg.lstsq(A, b, rcond=None)[0]
            U = U @ expi(np.einsum('a,aij->ij', step[:len(Gfull)], np.array(Gfull)))
        X = C @ U @ C.conj().T
        if np.linalg.norm(U.T - np.vdot(X, U.T) / abs(np.vdot(X, U.T)) * X) < 1e-10:
            return U
    return None


def stratum_dimension(U, C, G):
    d = U.shape[0]
    X = C @ U @ C.conj().T
    lam = np.vdot(X, U.T)
    lam /= abs(lam)
    Jm = [(U @ (1j * g)).T.ravel() - lam * (C @ U @ (1j * g) @ C.conj().T).ravel() for g in G]
    Jm.append((-1j * lam * X).ravel())
    Jm = np.array(Jm).T
    A = np.concatenate([Jm.real, Jm.imag])
    s = np.linalg.svd(A, compute_uv=False)
    rank = int((s > 1e-8 * s[0]).sum())
    return (len(G) + 1) - rank                                    # kernel dim in (PU(d) tangent, lambda)


def hessian(U, t, Cl, G, h=1e-5):
    n = len(G)
    H = np.zeros((n, n))
    for k in range(n):
        Up = U @ expi(h * G[k])
        Um = U @ expi(-h * G[k])
        H[k] = (j_and_grad(Up, t, Cl, G)[1] - j_and_grad(Um, t, Cl, G)[1]) / (2 * h)
    return (H + H.T) / 2


def local_analysis(d, t, seed=1):
    rng = np.random.default_rng(seed)
    Cl = P7.clifford_group(d)
    G = hermitian_basis(d)
    out = []
    for i, size in twisted_classes(Cl):
        U = point_on_stratum(Cl[i], G, rng)
        if U is None:
            out.append(dict(twisted_class_size=size, nonempty=False))
            continue
        dim = stratum_dimension(U, Cl[i], G)
        ev = np.sort(np.abs(np.linalg.eigvalsh(hessian(U, t, Cl, G))))[::-1]
        codim = len(G) - dim
        out.append(dict(twisted_class_size=size, nonempty=True, stratum_dim=dim, codim=codim,
                        J=j_and_grad(U, t, Cl, G)[0], hessian_eigs=[float(f"{x:.3e}") for x in ev],
                        hessian_rank=int((ev > 1e-6 * max(ev[0], 1e-300)).sum()) if ev[0] > 1e-9 else 0,
                        rank_equals_codim=bool(ev[0] > 1e-9 and (ev > 1e-6 * ev[0]).sum() == codim)))
    return out


def run():
    res = dict(pass_id=11369)
    res["local_d3_J6"] = local_analysis(3, 3)
    print(json.dumps(res["local_d3_J6"], indent=1), flush=True)
    res["global_d3_J6"] = global_search(3, 3, 2000, 113690000)
    print(res["global_d3_J6"], flush=True)
    res["global_d3_J8"] = global_search(3, 4, 1000, 113692000)
    print(res["global_d3_J8"], flush=True)
    res["local_d3_J8"] = local_analysis(3, 4)
    print(json.dumps(res["local_d3_J8"], indent=1), flush=True)
    res["control_d5_J4"] = global_search(5, 2, 400, 113695000)
    print(res["control_d5_J4"], flush=True)
    return res


def refine_spurious(seed=2, iters=3):
    """the spurious zero reached from Haar seed `seed`, refined by Gauss-Newton on the full residual
    F(U) = A_3(U) - A_3(U^T) (729 x 729); reports ||F|| per step, the singular values of dF, J_8, J_10, distance"""
    Cl = P7.clifford_group(3)
    G = hermitian_basis(3)

    def pi3(V):
        a = np.kron(np.kron(V, V), V)
        return np.kron(a, a.conj())

    def F(U):
        A = sum(pi3(C @ U @ C.conj().T) for C in Cl) / len(Cl)
        B = sum(pi3(C @ U.T @ C.conj().T) for C in Cl) / len(Cl)
        return (A - B).ravel()

    U, _ = descend(haar(3, np.random.default_rng(seed)), 3, Cl, G)
    norms = []
    for _ in range(iters):
        f0 = F(U)
        h = 1e-6
        Jm = np.array([(F(U @ expi(h * g)) - F(U @ expi(-h * g))) / (2 * h) for g in G]).T
        A = np.concatenate([Jm.real, Jm.imag])
        sv = np.linalg.svd(A, compute_uv=False)
        norms.append(float(np.linalg.norm(f0)))
        step = np.linalg.lstsq(A, -np.concatenate([f0.real, f0.imag]), rcond=1e-9)[0]
        U = U @ expi(np.einsum('a,aij->ij', step, G))
    norms.append(float(np.linalg.norm(F(U))))
    return dict(seed=seed, residual_norms=norms, dF_singular_values_relative=[float(x) for x in sv / sv[0]],
                dF_rank=int((sv > 1e-6 * sv[0]).sum()), J6=P7.J(U, 3, Cl), J8=P7.J(U, 4, Cl), J10=P7.J(U, 5, Cl),
                distance_to_reversible=rev_distance(U, Cl))


def _deep_job(args):
    k, seed, n = args
    import w33_pass11213_cubic_t_violation as P1
    import w33_pass11252_exact_reversibility as R
    R.WEYL[1] = R.Weyl(1)
    CF = np.array(P1.clifford1())
    rng = np.random.default_rng(seed)
    out = dict(words=0, violating=0, j6_zero_but_violating=0, j6_positive_but_reversible=0, min_j6_violator=np.inf,
               undecided=0, numerically_ambiguous=0, argmin_word=None)
    for _ in range(n):
        U = np.eye(3, dtype=complex)
        word = [int(x) for x in rng.integers(216, size=k)]
        for c in word:
            U = CF[c] @ P1.T1 @ U
        try:
            v = R.decide(U, 1, rng)[0]
        except AssertionError:                        # floating magnitudes too close to class: reported, not guessed
            out["numerically_ambiguous"] += 1
            continue
        if v is None:
            out["undecided"] += 1
            continue
        j = P7.J(U, 3, CF)
        out["words"] += 1
        if v is False:
            out["violating"] += 1
            if j < out["min_j6_violator"]:
                out["min_j6_violator"], out["argmin_word"] = j, word
            out["j6_zero_but_violating"] += j < 1e-9
        else:
            out["j6_positive_but_reversible"] += j > 1e-9
    return out


def deep_words(ks=(5, 6, 8, 10, 12, 16, 20), per_k=200000, nproc=3):
    """does J_6 stay complete on DEEP Clifford+T words, where the spurious family is not excluded?"""
    res = {}
    with Pool(nproc) as pool:
        for k in ks:
            parts = pool.map(_deep_job, [(k, 1136900 + 1000 * k + i, per_k // 20) for i in range(20)])
            agg = {key: sum(p[key] for p in parts) for key in parts[0] if key not in ("min_j6_violator", "argmin_word")}
            best = min(parts, key=lambda p: p["min_j6_violator"])
            agg["min_j6_violator"] = float(best["min_j6_violator"])
            agg["argmin_word_clifford_indices"] = best["argmin_word"]
            res[str(k)] = agg
            print(k, agg, flush=True)
    return res


def main():
    if "--deep" in sys.argv:
        res = json.load(open(OUT))
        res["deep_words"] = deep_words()
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    if "--refine" in sys.argv:
        res = json.load(open(OUT))
        res["refined_spurious_zero"] = refine_spurious()
        print(res["refined_spurious_zero"], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
