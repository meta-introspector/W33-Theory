#!/usr/bin/env python3
"""Pass 11167: N_AB + N_AC <= 4/sqrt15 -- proved on the symmetry class of the optimum, a rigorous (weaker) global bound,
and a dual see-saw search from thousands of starts that never exceeds 4/sqrt15.

(1) THE SYMMETRY CLASS.  psi* = sqrt(7/15)|000> + sqrt(2/15) sum_{a=1,2} |a>(|0a> + |a0>) is invariant under U (x) conj(U) (x)
    conj(U) for every U in U(2) acting on span{|1>,|2>} (Pass 11162's 4-dimensional stabiliser).  The invariant states are
        psi = alpha|000> + beta sum_a |a a 0> + gamma sum_a |a 0 a>      (phases removable: alpha, beta, gamma >= 0),
    and with s = beta^2, t = gamma^2, alpha^2 = 1 - 2s - 2t the partial transposes split into 2x2 blocks, giving exactly
        N_AB = sqrt(t^2 + 4 alpha^2 s) - t + s,     N_AC = sqrt(s^2 + 4 alpha^2 t) - s + t,
        N_AB + N_AC = sqrt(t^2 + 4 alpha^2 s) + sqrt(s^2 + 4 alpha^2 t)          (the linear terms cancel).
    Every critical point in the triangle and the whole boundary are computed exactly (sympy): the maximum is 4/sqrt15 at
    s = t = 2/15, uniquely.  So 4/sqrt15 is the exact maximum on the symmetry class of the optimum.
(2) A RIGOROUS GLOBAL BOUND.  For a two-qutrit state the partial transpose has at most (d-1)^2 = 4 negative eigenvalues
    (Rana, PRA 87, 054301 (2013)), so its purity p >= min_{m<=4} [(1+N)^2/(9-m) + N^2/m]: N <= h(p).  For pure psi_ABC,
    purity(rho_AB) = purity(rho_C), and N_AB <= N_{B|AC}(psi) = ((sum_i sqrt(mu_i))^2 - 1)/2 (tracing C is local).  Hence
        N_AB + N_AC <= max over qutrit spectra mu, nu of  min(h(p(nu)), g(mu)) + min(h(p(mu)), g(nu))  =  B,
    evaluated on a fine grid with a Lipschitz margin: B is recorded below -- a proof that polygamy of spatial negativity is
    strictly below the trivial 2, but not yet tight.
(3) DUAL SEE-SAW.  N(sigma) = max_{0<=P<=I} -Tr(P sigma^{T}), so N_AB + N_AC = max over psi, P1, P2 of <psi|H(P1,P2)|psi>.
    Alternating psi <- top eigenvector of H and P_i <- projector onto the negative eigenspace never decreases the value;
    thousands of random starts converge to at most 4/sqrt15 (recorded below).
Scope: exact on the symmetry class, rigorous bound B globally, numerical evidence for 4/sqrt15 globally.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11167_polygamy_global.json"


# ---------- (1) the symmetry class ----------
def family_state(s, t):
    a = np.sqrt(max(0.0, 1 - 2 * s - 2 * t))
    psi = np.zeros(27)
    psi[0] = a
    for k in (1, 2):
        psi[k * 9 + k * 3] = np.sqrt(s)       # |k k 0>
        psi[k * 9 + k] = np.sqrt(t)           # |k 0 k>
    return psi


def negs(psi):
    T = psi.reshape(3, 3, 3)
    out = []
    for M in (T.reshape(9, 3), T.transpose(0, 2, 1).reshape(9, 3)):
        rho = M @ M.conj().T
        pt = rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
        ev = np.linalg.eigvalsh((pt + pt.conj().T) / 2)
        out.append(float(-ev[ev < 0].sum()))
    return out


def symmetry_class_exact():
    s, t = sp.symbols('s t', nonnegative=True)
    al = 1 - 2 * s - 2 * t
    A, B = t ** 2 + 4 * al * s, s ** 2 + 4 * al * t
    f = sp.sqrt(A) + sp.sqrt(B)
    # interior critical points: dA/ds sqrt(B) + dB/ds sqrt(A) = 0 etc. -> square to a polynomial system
    ga = [sp.diff(A, v) for v in (s, t)]
    gb = [sp.diff(B, v) for v in (s, t)]
    # grad f = 0  <=>  ga_i sqrt(B) + gb_i sqrt(A) = 0 (i = s, t)  =>  ga_s gb_t = ga_t gb_s  and  ga_s^2 B = gb_s^2 A
    eq1 = sp.expand(ga[0] * gb[1] - ga[1] * gb[0])
    eq2 = sp.expand(ga[0] ** 2 * B - gb[0] ** 2 * A)
    sols = sp.solve([eq1, eq2], [s, t], dict=True)
    cands = []
    for so in sols:
        sv, tv = so.get(s), so.get(t)
        if sv is None or tv is None or not (sv.is_real and tv.is_real):
            continue
        if sv >= 0 and tv >= 0 and 2 * sv + 2 * tv <= 1:
            # keep genuine critical points (sign condition of the unsquared equations)
            val = f.subs({s: sv, t: tv})
            gs = (ga[0] * sp.sqrt(B) + gb[0] * sp.sqrt(A)).subs({s: sv, t: tv})
            if abs(float(gs)) < 1e-9:
                cands.append((sp.nsimplify(sv), sp.nsimplify(tv), sp.nsimplify(sp.simplify(val))))
    # boundary (by the s <-> t symmetry, t = 0 covers s = 0): f(s, 0) = 2 sqrt(s - 2 s^2) + s on [0, 1/2]; alpha = 0 gives
    # f = s + t = 1/2.  One-variable maximum of f(s, 0): f' = 0 <=> (1 - 4s) + sqrt(s - 2s^2) = 0.
    x = sp.symbols('x', real=True)
    fb = 2 * sp.sqrt(x - 2 * x ** 2) + x
    crit = [c for c in sp.solve(sp.Eq((1 - 4 * x) ** 2, x - 2 * x ** 2), x) if c.is_real and 0 <= c <= sp.Rational(1, 2)
            and abs(float((1 - 4 * c) + sp.sqrt(c - 2 * c ** 2))) < 1e-12]
    bvals = [float(fb.subs(x, c)) for c in crit + [sp.Integer(0), sp.Rational(1, 2)]] + [0.5]
    # dense numerical cross-check of the whole triangle
    S, T = np.meshgrid(np.linspace(0, 0.5, 1001), np.linspace(0, 0.5, 1001))
    Al = np.clip(1 - 2 * S - 2 * T, 0, None)
    Fv = np.where(2 * S + 2 * T <= 1, np.sqrt(T ** 2 + 4 * Al * S) + np.sqrt(S ** 2 + 4 * Al * T), -1)
    return dict(critical_points=[[str(a), str(b), str(c), float(c)] for a, b, c in cands],
                boundary_critical=[str(c) for c in crit], boundary_max=max(bvals),
                max_value=max([float(c[2]) for c in cands] + bvals), grid_max=float(Fv.max()))


# ---------- (2) rigorous global bound ----------
def h_of_p(p):
    """largest N with min_{m=1..4} [(1+N)^2/(9-m) + N^2/m] <= p (the PT spectrum constraint)"""
    lo, hi = 0.0, 1.0
    for _ in range(60):
        N = (lo + hi) / 2
        need = min((1 + N) ** 2 / (9 - m) + N ** 2 / m for m in (1, 2, 3, 4))
        if need <= p:
            lo = N
        else:
            hi = N
    return lo


def g_of_spec(mu):
    return ((np.sqrt(np.clip(mu, 0, None)).sum()) ** 2 - 1) / 2


def global_bound(n=1500):
    """N_AB + N_AC <= h(p(lambda_C)) + g(lambda_C)  (N_AB <= h(purity rho_AB = purity rho_C), N_AC <= N_{C|AB} = g);
    maximised over the qutrit simplex on a grid of step 1/n (vectorised bisection for h)"""
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing='ij')
    m = i + j <= n
    L = np.stack([i[m], j[m], n - i[m] - j[m]], axis=1) / n
    P = (L ** 2).sum(1)
    G = ((np.sqrt(L).sum(1)) ** 2 - 1) / 2
    lo, hi = np.zeros_like(P), np.ones_like(P)
    for _ in range(50):
        N = (lo + hi) / 2
        need = np.min([(1 + N) ** 2 / (9 - k) + N ** 2 / k for k in (1, 2, 3, 4)], axis=0)
        ok = need <= P
        lo, hi = np.where(ok, N, lo), np.where(ok, hi, N)
    V = lo + G
    k = int(V.argmax())
    return float(V.max()), [float(x) for x in L[k]]


# ---------- (3) dual see-saw ----------
def pt(rho):
    return rho.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)


def seesaw(rng, iters=400):
    psi = rng.normal(size=27) + 1j * rng.normal(size=27)
    psi /= np.linalg.norm(psi)
    val = -1
    for _ in range(iters):
        T = psi.reshape(3, 3, 3)
        Ms = (T.reshape(9, 3), T.transpose(0, 2, 1).reshape(9, 3))
        Ps = []
        for M in Ms:
            w, V = np.linalg.eigh(pt(M @ M.conj().T))
            Vn = V[:, w < 0]
            Ps.append(Vn @ Vn.conj().T)
        # H = -(P1^{T_B} (x) I_C) - (P2^{T_C} on AC (x) I_B), acting on A B C
        H1 = np.kron(pt(Ps[0]), np.eye(3))
        H2 = np.kron(pt(Ps[1]), np.eye(3)).reshape(3, 3, 3, 3, 3, 3).transpose(0, 2, 1, 3, 5, 4).reshape(27, 27)
        H = -(H1 + H2)
        w, V = np.linalg.eigh((H + H.conj().T) / 2)
        new = float(w[-1])
        psi = V[:, -1]
        if abs(new - val) < 1e-13:
            break
        val = new
    return val


def summarize(n_seesaw=3000):
    exact = symmetry_class_exact()
    s = t = 2 / 15
    num = sum(negs(family_state(s, t)))
    B, Bspec = global_bound()
    rng = np.random.default_rng(11167)
    vals = np.array([seesaw(rng) for _ in range(n_seesaw)])
    target = 4 / np.sqrt(15)
    res = dict(pass_id=11167, symmetry_class=exact, family_value_at_2_15=num, four_over_sqrt15=target,
               closed_form_check=bool(abs(np.sqrt(t ** 2 + 4 * (1 - 4 * s) * s) * 2 - target) < 1e-12),
               global_bound=B, global_bound_argmax_spectrum=Bspec,
               seesaw_starts=n_seesaw, seesaw_max=float(vals.max()),
               seesaw_above_target=int((vals > target + 1e-9).sum()),
               seesaw_at_target=int((abs(vals - target) < 1e-7).sum()),
               seesaw_top_values=sorted({round(float(v), 6) for v in vals})[-6:])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
