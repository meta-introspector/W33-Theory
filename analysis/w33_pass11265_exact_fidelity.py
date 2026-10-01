"""Pass 11265: the exact best time-reversal fidelity from Weyl coefficients, and the spectrum of the optimal residual.

FORMULA (validated against the matrix twirl).  With U = 3^{-n} sum_p c(p) W(p) in symmetric Weyl operators
(Pass 11252), for an anti-symplectic L the best reversal in the Weyl coset of V_L has
    F(L) = max_b | sum_p conj(c(p)) c(Lp) omega^{<b, Lp>} | / D^2,      D = 3^n,
and F_T(U) = max_L F(L).  The inner maximum is a Fourier transform over F_3^{2n} (all b at once).

NUMERICAL TREE EXHAUSTION (branch and bound over L).  L is built basis vector by basis vector.  For a partial map on the
span S of the first k basis vectors, every completion satisfies
    D^2 F(L) <= max_{chi character of L(S)} | sum_{p in S} g(p) chi(Lp) | + max_pairing sum_{p not in S} |c(p)||c(Lp)|,
with g(p) = conj(c(p)) c(Lp): the span part is maximised over its 3^k characters exactly, and the rest by the
rearrangement inequality (L maps the complement of S bijectively onto the complement of L(S)).  Best-first search: the
largest open bound is an upper bound in exact arithmetic.  This implementation uses floating complex arithmetic and
prunes at 1e-12: it returns a numerically exhausted tree or a floating bracket if a node budget is hit, not a directed-
rounding interval certificate.

RESIDUAL SPECTRUM.  For two qutrits, every optimal reversal (up to 600 per level) of every level of the Pass 11235
stream is examined: does SOME optimal residual R = V U^* V^dag U have eigenphases in mu_9 / mu_27 / mu_81, and for the
F_min^2 words, is it the sumset {-1,0,1}+{-1,0,1} found in Pass 11251?
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11252_exact_reversibility as R  # noqa: E402

OUT = ROOT / "data" / "w33_pass11265_exact_fidelity.json"


def full_F(wl, c, L):
    D = wl.D
    Lp = (wl.labels @ L.T) % 3
    g = np.conj(c) * c[Lp @ wl.pow3]
    forms = (wl.labels @ wl.Om @ Lp.T) % 3
    return float(np.abs((R.OM ** forms) @ g).max() / D ** 2)


def best_fidelity(wl, U, incumbent=0.0, node_budget=2_000_000):
    """Numerical max_L F(L) by branch and bound; returns (value, tree exhausted?, nodes)."""
    n2 = 2 * wl.n
    D = wl.D
    c = wl.coeffs(U)
    a = np.abs(c)
    labels = wl.labels
    # basis: greedily the largest |c| labels, linearly independent
    order = np.argsort(-a)
    basis, span = [], {0}
    for i in order:
        if i not in span:
            v = labels[i]
            basis.append(v)
            span |= {wl.index(labels[s] + lam * v) for s in span for lam in (1, 2)}
            if len(basis) == n2:
                break
    B = np.array(basis)
    Binv = R._inv_mod3(B.T)
    lams = [np.array(list(itertools.product(range(3), repeat=k)), dtype=int).reshape(3 ** k, k) for k in range(n2 + 1)]
    OmB = B @ wl.Om
    state = dict(best=incumbent, nodes=0, exact=True, L=None)
    a_sorted_all = np.sort(a)[::-1]

    def bound(k, Im):
        coords = lams[k]
        src = ((coords @ B[:k]) % 3) @ wl.pow3 if k else np.array([0])
        dst = ((coords @ Im) % 3) @ wl.pow3 if k else np.array([0])
        g = np.conj(c[src]) * c[dst]
        # characters of L(S) ~ F_3^k: chi_t(L p) = omega^{t . lambda(p)}
        if k:
            ph = R.OM ** ((lams[k] @ coords.T) % 3)          # (3^k chars, 3^k points)
            span_part = np.abs(ph @ g).max()
        else:
            span_part = abs(g[0])
        ms = np.ones(len(a), bool)
        ms[src] = False
        md = np.ones(len(a), bool)
        md[dst] = False
        rest = float(np.sort(a[ms])[::-1] @ np.sort(a[md])[::-1])
        return (span_part + rest) / D ** 2

    import heapq

    def children(images):
        k = len(images)
        Im = np.array(images) if k else np.zeros((0, n2), int)
        bk = B[k]
        cand = labels
        if k:
            need = (-(OmB[k] @ B[:k].T)) % 3
            have = (cand @ wl.Om @ Im.T) % 3
            cand = cand[(have == need).all(axis=1)]
            img_span = set((((lams[k] @ Im) % 3) @ wl.pow3).tolist())
            cand = [v for v in cand if wl.index(v) not in img_span]
        else:
            cand = list(cand[1:])
        return [images + [v] for v in cand]

    heap = [(-bound(0, np.zeros((0, n2), int)), 0, [])]
    counter = 1
    while heap:
        negub, _, images = heapq.heappop(heap)
        if -negub <= state["best"] + 1e-12:
            heap = []
            break
        state["nodes"] += 1
        if state["nodes"] > node_budget:
            state["exact"] = False
            heapq.heappush(heap, (negub, counter, images))
            break
        for ch in children(images):
            k = len(ch)
            Im = np.array(ch)
            if k == n2:
                L = (Im.T @ Binv) % 3
                f = full_F(wl, c, L)
                if f > state["best"]:
                    state["best"], state["L"] = f, L
            else:
                ub = bound(k, Im)
                if ub > state["best"] + 1e-12:
                    heapq.heappush(heap, (-ub, counter, ch))
                    counter += 1
    state["upper"] = max(state["best"], -heap[0][0]) if heap else state["best"]
    return state["best"], state["exact"], state["nodes"], state["L"], state["upper"]


def validate_formula(rng, trials=40):
    wl = R.WEYL[2]
    worst = 0.0
    for _ in range(trials):
        U = R.random_word(2, int(rng.integers(1, 4)), rng)
        c = wl.coeffs(U)
        VM = R.random_clifford(2, rng)
        M = np.zeros((4, 4), int)
        for j in range(4):
            e = np.zeros(4, int)
            e[j] = 1
            img = VM @ wl.W[wl.index(e)] @ VM.conj().T
            M[:, j] = wl.labels[int(np.argmax(np.abs(np.einsum('pij,ij->p', wl.W.conj(), img))))]
        L = (-M @ wl.J) % 3
        direct = max(abs(np.trace(wl.W[i] @ VM @ U.conj() @ VM.conj().T @ wl.W[i].conj().T @ U)) for i in range(81)) / 9
        worst = max(worst, abs(direct - full_F(wl, c, L)))
    return worst


SUMSET = sorted((i + j) % 9 for i in (-1, 0, 1) for j in (-1, 0, 1))


def _quantum_order(ev, Ns=(9, 27, 81)):
    ph = np.angle(ev / ev[0])
    for N in Ns:
        k = ph / (2 * np.pi / N)
        if np.allclose(k, np.round(k), atol=1e-6):
            return N
    return None


def _is_sumset(ev):
    """eigenvalues = global phase x zeta^{i+j}, i, j in {-1, 0, 1}?"""
    for e0 in ev:
        k = np.angle(ev / e0) / (2 * np.pi / 9)
        if not np.allclose(k, np.round(k), atol=1e-6):
            return False
        ks = np.round(k).astype(int) % 9
        for s0 in range(9):
            if sorted((ks + s0) % 9) == SUMSET:
                return True
    return False


def residual_spectra(reps, max_maximisers=600):
    """for every level of the Pass 11235 stream: does SOME optimal reversal leave a residual V U* V^dag U with
    eigenphases in mu_9 / mu_27 / mu_81?  For F_min^2 levels: is that residual spectrum the sumset {-1,0,1}+{-1,0,1}?"""
    import w33_pass11253_cubic_field_levels as L53
    out = {}
    for d, fac in L53.two_qutrit_stream():
        U = L53.prod(fac)
        vals = P2_vals(U, reps)
        f = vals.max()
        key = round(float(f), 9)
        if f > 1 - 1e-9 or key in out:
            continue
        arg = np.argwhere(vals > f - 1e-9)[:max_maximisers]
        best, sumset = None, False
        for r, a_ in arg:
            V = P2.PA[a_] @ reps[r]
            ev = np.linalg.eigvals(V @ U.conj() @ V.conj().T @ U)
            q = _quantum_order(ev)
            if q is not None:
                best = q if best is None else min(best, q)
                sumset = sumset or _is_sumset(ev)
        out[key] = dict(maximisers_checked=int(len(arg)), quantised_residual_exists=best, sumset_residual=sumset)
    return out


def fmin2_sumset_check(reps):
    """every F_min^2 word of Pass 11238's target-sector set: some optimal residual is the sumset"""
    T, I3 = P2.T, P2.I3
    words = {"(I(x)T)SUM(I(x)T^2)": np.kron(I3, T) @ P2.SUM @ np.kron(I3, T @ T),
             "(I(x)T^2)SUM(I(x)T)": np.kron(I3, T @ T) @ P2.SUM @ np.kron(I3, T),
             "(I(x)T)SUM^2(I(x)T^2)": np.kron(I3, T) @ P2.SUM @ P2.SUM @ np.kron(I3, T @ T)}
    out = {}
    for k, U in words.items():
        vals = P2_vals(U, reps)
        f = vals.max()
        arg = np.argwhere(vals > f - 1e-9)[:600]
        ss = any(_is_sumset(np.linalg.eigvals(P2.PA[a_] @ reps[r] @ U.conj() @ (P2.PA[a_] @ reps[r]).conj().T @ U))
                 for r, a_ in arg)
        out[k] = dict(F_T=float(f), is_Fmin2=bool(abs(f - FMIN2) < 1e-9), sumset_residual=bool(ss))
    return out


FMIN2 = ((1 + 2 * np.cos(2 * np.pi / 9)) / 3) ** 2


def P2_vals(U, reps):
    import w33_pass11251_exact_reversal as E
    return E.float_maximisers(U, reps)


def run():
    rng = np.random.default_rng(11265)
    for n in (2, 3):
        R.WEYL[n] = R.Weyl(n)
    res = dict(pass_id=11265)
    res["formula_max_deviation_vs_twirl"] = validate_formula(rng)
    print("formula", res["formula_max_deviation_vs_twirl"], flush=True)
    reps = np.load(P2.CACHE)
    # branch and bound vs exhaustive on two qutrits
    T, I3 = P2.T, P2.I3
    named2 = {
        "(T(x)T)SUM": np.kron(T, T) @ P2.SUM,
        "(I(x)T)SUM(I(x)T^2)": np.kron(I3, T) @ P2.SUM @ np.kron(I3, T @ T),
    }
    val = []
    for k in range(12):
        named2[f"random d={1 + k % 4} #{k}"] = R.random_word(2, 1 + k % 4, rng)
    for name, U in named2.items():
        ex = float(P2_vals(U, reps).max())
        t0 = time.time()
        bb, exact, nodes, _, _ = best_fidelity(R.WEYL[2], U)
        val.append(dict(word=name, exhaustive=ex, branch_and_bound=bb, search_exhausted=exact,
                        numerical_pruning_tolerance=1e-12, nodes=nodes,
                        seconds=round(time.time() - t0, 2), agree=bool(abs(ex - bb) < 1e-9)))
        print(val[-1], flush=True)
    res["two_qutrit_bb_validation"] = val
    res["bb_agrees_all"] = all(v["agree"] for v in val)
    # three qutrits
    S12, S23 = R.sum_gate(3, 0, 1), R.sum_gate(3, 1, 2)
    k3 = lambda x, y, z: np.kron(np.kron(x, y), z)  # noqa: E731
    named3 = {
        "(T(x)T)SUM (x) I": (np.kron(np.kron(T, T) @ P2.SUM, I3), 0.8440296287459851),
        "(T(x)T(x)T) SUM12": (k3(T, T, T) @ S12, 0.8440296287459851),
        "(T(x)T(x)T) SUM23 SUM12": (k3(T, T, T) @ S23 @ S12, 0.23746200473369541),   # 11239's lower bound
    }
    out3 = {}
    for name, (U, inc) in named3.items():
        t0 = time.time()
        bb, exact, nodes, _, ub = best_fidelity(R.WEYL[3], U, incumbent=inc - 1e-9, node_budget=20_000)
        out3[name] = dict(F_T_lower=bb, F_T_upper=ub, search_exhausted=exact,
                          numerical_pruning_tolerance=1e-12, nodes=nodes, seconds=round(time.time() - t0, 1))
        print(name, out3[name], flush=True)
    res["three_qutrit"] = out3
    res["residual_spectra_two_qutrit"] = residual_spectra(reps)
    rs = res["residual_spectra_two_qutrit"]
    res["levels_with_a_quantised_optimal_residual"] = sorted(k for k, v in rs.items() if v["quantised_residual_exists"])
    res["levels_total"] = len(rs)
    res["fmin2_words"] = fmin2_sumset_check(reps)
    print(res["levels_with_a_quantised_optimal_residual"], res["fmin2_words"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "residual_spectra_two_qutrit"}, indent=1, default=str))


if __name__ == "__main__":
    main()
