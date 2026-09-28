#!/usr/bin/env python3
"""One-loop engine for the SO(16)xSO(16) string on T6/(Z2W x Z3) with order-3 Wilson lines (Pass 11106).

Z(tau) = (1/6) sum_{g,h in G} Z(g,h),  G = Z2W x Z3 = Z6 (twists v0 = (0,1,1,1), v = (0,1/3,1/3,-2/3); shifts V0, V;
Wilson lines W1 = W2, W3 = W4, W5 = W6 in the E8xE8 basis of the orbifolder model files):
  * the 4 beta-type sectors (g, h in {1, beta}, beta = the Witten element, trivial on T6) = (1/3) Z[SO16^2 on T6]:
    Hamiltonian lattice sums over the Narain lattice Gamma_{6,22} with Wilson lines.  Narain vectors P(pi, m, n):
    l = pi + A n,  w_+- = m - A^T pi - 1/2 A^T A n + B n +- G n,  p_{L,R}^2 = 1/2 w^T G^-1 w  (l^2 + pL^2 - pR^2 = pi^2 + 2mn).
    The Wilson-line torus sums are decomposed into the cosets pi.W mod 1 (E8 theta functions with characteristics);
    the Wilson-line-free torus factorises as Gamma_{2,2}(T*, rho).  Vacuum phases: Z[0,1](-1/tau) = Z[1,0](tau) and
    Z[1,0](tau+1) = -H[1,1](tau) (fixed numerically, then asserted).
  * the other 32 sectors are one SL(2,Z) orbit per fixed point: the seed B_f(0,1) (untwisted sector with insertion
    beta theta_f, local shift V0 + V + sum_t c_t W_t) is Gamma_1(6)-invariant (checked), and the 24 coset images
    B_f(0,1)(gamma tau) give every twisted sector; the 8 pure-Z3 pairs vanish identically (supersymmetric).
Integral over F: strip 1 <= tau2 <= 3 (tau1 midpoint rule, exact for integer spins; Gauss-Legendre in tau2) + cap
(|tau| >= 1) + tail from a fit tau2 <Z>_tau1 = c0 + c1 exp(-a tau2).  Validated on the 10D SO(16)xSO(16) integral of
Pass 11099 (I = -725.97).  Lambda_4 = -(1/2) M^4 I  (same normalisation as Pass 11099).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import gcd

import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit


NTH = 60


def qpow(tau, x):
    return np.exp(2j * np.pi * tau * x)


def theta(a, b, tau, N=None):
    """theta[a,b](tau) = sum_n q^{(n+a)^2/2} e^{2 pi i (n+a) b}; a, b arrays (broadcast) ; tau scalar"""
    if N is None:
        N = int(np.ceil(np.sqrt(40.0 / (np.pi * tau.imag)))) + 3
    n = np.arange(-N, N + 1).reshape((-1,) + (1,) * np.ndim(a))
    x = n + np.asarray(a, float)
    return np.sum(np.exp(1j * np.pi * tau * x * x + 2j * np.pi * x * np.asarray(b, float)), axis=0)


def eta(tau):
    q = np.exp(2j * np.pi * tau)
    K = int(40.0 / (2 * np.pi * tau.imag)) + 10
    k = np.arange(1, K + 1)
    return np.exp(2j * np.pi * tau / 24) * np.prod(1 - q ** k)




def fermions(h, g, tau, eps=+1):
    """R(h,g) = sum_{r in V u S} sign(r) q^{(r+h)^2/2} e^{2 pi i (r+h).g}; V: Z^4 odd sum (+), S: (Z+1/2)^4 (-) with parity
    projection (1 + eps (-1)^{sum r})/2"""
    h, g = np.asarray(h, float), np.asarray(g, float)
    ph = np.exp(-1j * np.pi * np.sum(h))
    V = 0.5 * (np.prod(theta(h, g, tau)) - ph * np.prod(theta(h, g + 0.5, tau)))
    S = 0.5 * (np.prod(theta(h + 0.5, g, tau)) + eps * ph * np.prod(theta(h + 0.5, g + 0.5, tau)))
    return V - S


def e8theta(w, z, tau):
    """sum_{pi in E8} q^{(pi+w)^2/2} e^{2 pi i pi.z}"""
    w, z = np.asarray(w, float), np.asarray(z, float)
    ph = np.exp(-1j * np.pi * np.sum(w))
    t = (np.prod(theta(w, z, tau)) + ph * np.prod(theta(w, z + 0.5, tau)) + np.prod(theta(w + 0.5, z, tau))
         + ph * np.prod(theta(w + 0.5, z + 0.5, tau)))
    return 0.5 * np.exp(-2j * np.pi * np.dot(w, z)) * t


def gauge16(w, z, tau):
    return e8theta(w[:8], z[:8], tau) * e8theta(w[8:], z[8:], tau)


def parse_model(block):
    rows = [[float(Fr(x.strip())) for x in r.split(',')] for r in block]
    V0, V = np.array(rows[0]), np.array(rows[1])
    W = [np.array(r) for r in rows[2:8]]
    return V0, V, W


V0V = np.array([0, 1, 1, 1], float)
VV = np.array([0, 1 / 3, 1 / 3, -2 / 3])


def seed_block(tau, U, eps):
    """B_f(0,1): untwisted sector with insertion gamma_f = beta theta_f (twist u = v0 + v, shift U)"""
    u = V0V + VV
    t2 = tau.imag
    et = eta(tau)
    bos = 1.0
    for i in (1, 2, 3):
        bos *= abs(et / theta(0.5, 0.5 + u[i], tau)) ** 2
    F = np.conj(fermions(np.zeros(4), -u, tau, eps)) / np.conj(et) ** 4       # right movers: qbar
    G = gauge16(np.zeros(16), U, tau) / et ** 16
    return bos * F * G / (t2 * abs(et) ** 4)


def cosets(N=6):
    """gamma in SL(2,Z) with (0,1) gamma running over all (a,b) mod N with gcd(a,b,N) = 1"""
    out = {}
    for a in range(N):
        for b in range(N):
            if gcd(gcd(a, b), N) != 1:
                continue
            # find c = a + N x, d = b + N y coprime
            done = False
            for x in range(0, 6):
                for y in range(0, 6):
                    c, d = a + N * x, b + N * y
                    if gcd(c, d) == 1:
                        # alpha d - beta c = 1
                        g, s, t = ext_gcd(d, c)          # s d + t c = 1
                        alpha, beta = s, -t
                        assert alpha * d - beta * c == 1
                        out[(a, b)] = (alpha, beta, c, d)
                        done = True
                        break
                if done:
                    break
            assert done
    return out


def ext_gcd(a, b):
    if b == 0:
        return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, s, t = ext_gcd(b, a % b)
    return g, t, s - (a // b) * t


def act(gm, tau):
    a, b, c, d = gm
    return (a * tau + b) / (c * tau + d)


def fixed_point_shifts(V0, V, W):
    """local shifts U_f = V0 + V + sum_t c_t W_(2t-1), c in Z3^3"""
    out = []
    for c1 in range(3):
        for c2 in range(3):
            for c3 in range(3):
                out.append(V0 + V + c1 * W[0] + c2 * W[2] + c3 * W[4])
    return out


def Z_twisted(tau, shifts, cos, eps):
    tot = 0
    for gm in cos.values():
        tp = act(gm, tau)
        for U in shifts:
            tot += seed_block(tp, U, eps)
    return tot / 6.0


HEX = np.array([[1.0, -0.5], [-0.5, 1.0]])
EPSM = np.array([[0.0, 1.0], [-1.0, 0.0]])


def GB(T):
    """T = b + i sqrt(det G), hexagonal shape (U = rho)"""
    g = T.imag / np.sqrt(0.75)
    return g * HEX, T.real * EPSM


def torus_sum(tau, T, o, nw, M=None):
    """sum_{m in Z^2} q^{pL^2/2} qbar^{pR^2/2}, offsets o (2,), windings nw (2,)"""
    G, B = GB(T)
    Gi = np.linalg.inv(G)
    if M is None:
        M = int(np.ceil(np.sqrt(30.0 / (np.pi * tau.imag) * max(1.0, G[0, 0])))) + 3
    m = np.arange(-M, M + 1)
    m1, m2 = np.meshgrid(m, m, indexing='ij')
    base = np.stack([m1 - o[0], m2 - o[1]], -1) + B @ nw
    wp = base + G @ nw
    wm = base - G @ nw
    pl = 0.5 * np.einsum('...i,ij,...j->...', wp, Gi, wp)
    pr = 0.5 * np.einsum('...i,ij,...j->...', wm, Gi, wm)
    return np.sum(np.exp(1j * np.pi * tau * pl - 1j * np.pi * np.conj(tau) * pr))


def gamma22(tau, T, K=None):
    """Gamma_{2,2}(T, rho) without Wilson lines"""
    G, B = GB(T)
    if K is None:
        K = int(np.ceil(np.sqrt(30.0 / (np.pi * tau.imag) / min(1.0, G[0, 0] * 0.75)))) + 2
    tot = 0
    for n1 in range(-K, K + 1):
        for n2 in range(-K, K + 1):
            tot += torus_sum(tau, T, np.zeros(2), np.array([n1, n2], float))
    return tot


class WLPart:
    """Lambda'[s,t](tau) for the two Wilson-line tori (indices tw = [t1, t2] into W = [W1, W3, W5])"""

    def __init__(self, V0, Wt, Tw, K=4):
        self.V0, self.Wt, self.Tw, self.K = V0, Wt, Tw, K     # Wt: the two nonzero Wilson lines (E8xE8 vectors)
        self.WW = np.array([[a @ b for b in Wt] for a in Wt])
        self.V0W = np.array([V0 @ a for a in Wt])

    def lam(self, tau, s, t, phase_vac=1.0):
        K, V0, Wt = self.K, self.V0, self.Wt
        # E8^2 thetas: G(N mod 3, j) with z = t V0 + j.W
        Gc = {}
        for N1 in range(3):
            for N2 in range(3):
                w = s * V0 + N1 * Wt[0] + N2 * Wt[1]
                for j1 in range(3):
                    for j2 in range(3):
                        z = t * V0 + j1 * Wt[0] + j2 * Wt[1]
                        Gc[(N1, N2, j1, j2)] = gauge16(w, z, tau)
        # per torus: S_t(N_t, o) = sum_{n1 + n2 = N_t} torus_sum
        cache = {}

        def S(ti, Nt, o):
            key = (ti, Nt, round(o % 1.0, 10))
            if key not in cache:
                oo = o % 1.0
                acc = 0
                for n1 in range(-K, K + 1):
                    n2 = Nt - n1
                    if abs(n2) > K:
                        continue
                    acc += torus_sum(tau, self.Tw[ti], np.array([oo, oo]), np.array([n1, n2], float))
                cache[key] = acc
            return cache[key]
        tot = 0
        for N1 in range(-2 * K, 2 * K + 1):
            for N2 in range(-2 * K, 2 * K + 1):
                lamv = (N1 - N1 % 3) * Wt[0] + (N2 - N2 % 3) * Wt[1]
                for k1 in range(3):
                    for k2 in range(3):
                        g = 0
                        for j1 in range(3):
                            for j2 in range(3):
                                z = t * V0 + j1 * Wt[0] + j2 * Wt[1]
                                g += np.exp(-2j * np.pi * (j1 * k1 + j2 * k2) / 3) * np.exp(-2j * np.pi * lamv @ z) \
                                    * Gc[(N1 % 3, N2 % 3, j1, j2)]
                        g *= np.exp(2j * np.pi * t * s * (V0 @ V0)) / 9.0
                        if abs(g) < 1e-300:
                            continue
                        o = [k1 / 3 + s * self.V0W[0] + 0.5 * (N1 * self.WW[0, 0] + N2 * self.WW[0, 1]),
                             k2 / 3 + s * self.V0W[1] + 0.5 * (N1 * self.WW[1, 0] + N2 * self.WW[1, 1])]
                        tot += g * S(0, N1, o[0]) * S(1, N2, o[1])
        return tot * phase_vac


def Zbeta_block(tau, wl, s, t, Tstar, phase_vac=1.0):
    """E8xE8 on T6 with Wilson lines, sector beta^s, insertion beta^t (twist v0, shift V0)"""
    et = eta(tau)
    F = np.conj(fermions(s * V0V, -t * V0V, tau)) / np.conj(et) ** 4
    L = wl.lam(tau, s, t, phase_vac) * gamma22(tau, Tstar)
    return F * L / (et ** 22 * np.conj(et) ** 6) / (tau.imag * abs(et) ** 4)


TOP = 3.0
NT1 = 24


def grid():
    pts = []
    # strip
    x, w = np.polynomial.legendre.leggauss(24)
    t2 = 1 + (TOP - 1) * (x + 1) / 2
    w2 = w * (TOP - 1) / 2
    t1 = (np.arange(NT1) + 0.5) / NT1 - 0.5
    for a, wa in zip(t2, w2):
        for b in t1:
            pts.append((b + 1j * a, wa / NT1 / a ** 2, 'strip', a))
    # cap: tau1 in [-1/2, 1/2], tau2 in [sqrt(1 - tau1^2), 1]
    x1, w1 = np.polynomial.legendre.leggauss(24)
    xx, ww = np.polynomial.legendre.leggauss(12)
    for a, wa in zip(x1 / 2, w1 / 2):
        lo = np.sqrt(1 - a * a)
        for b, wb in zip(lo + (1 - lo) * (xx + 1) / 2, ww * (1 - lo) / 2):
            pts.append((a + 1j * b, wa * wb / b ** 2, 'cap', b))
    return pts


def integrate(values, pts):
    """values: Z at pts (complex).  Returns strip, cap, tail, total and the tau2 * <Z> profile"""
    strip = sum(v.real * p[1] for v, p in zip(values, pts) if p[2] == 'strip')
    cap = sum(v.real * p[1] for v, p in zip(values, pts) if p[2] == 'cap')
    prof = {}
    for v, p in zip(values, pts):
        if p[2] == 'strip':
            prof.setdefault(round(p[3], 12), []).append(v.real)
    t2s = sorted(prof)
    g = np.array([np.mean(prof[t]) * t for t in t2s])       # tau2 <Z>
    # tail: int_TOP^inf dtau2 tau2^-3 c(tau2); c ~ last value (+ exponential correction estimate)
    c_top = g[-1]
    tail = c_top / (2 * TOP ** 2)
    return dict(strip=strip, cap=cap, tail=tail, total=strip + cap + tail, c_top=c_top,
                profile=[(float(t), float(x)) for t, x in zip(t2s, g)])

def fitted(pts, vals, top=3.0):
    strip = sum(v * p[2] for v, p in zip(vals, pts) if p[3] == 'strip')
    cap = sum(v * p[2] for v, p in zip(vals, pts) if p[3] == 'cap')
    prof = {}
    for v, p in zip(vals, pts):
        if p[3] == 'strip':
            prof.setdefault(round(p[4], 12), []).append(v)
    t = np.array(sorted(prof)); g = np.array([np.mean(prof[x]) * x for x in t])
    sel = t > 2.0
    f = lambda x, c0, c1, a: c0 + c1 * np.exp(-a * x)
    try:
        (c0, c1, a), _ = curve_fit(f, t[sel], g[sel], p0=(g[-1], g[sel][0] - g[-1], 2.0), maxfev=20000)
        tail = quad(lambda x: f(x, c0, c1, a) / x ** 3, top, np.inf)[0]
    except Exception:
        c0, c1, a = g[-1], 0, 0; tail = g[-1] / (2 * top ** 2)
    return dict(strip=strip, cap=cap, tail=tail, I=strip + cap + tail, c0=c0, decay=a / (4 * np.pi), c_top=g[-1])



def Z_beta(tau, wl, Tstar):
    """(1/6) [Z(1,beta) + Z(beta,1) + Z(beta,beta)] = (1/3) Z[SO16^2 on T6]; Z(1,1) = 0 (E8xE8 is supersymmetric)"""
    return (Zbeta_block(tau, wl, 0, 1, Tstar) + Zbeta_block(tau, wl, 1, 0, Tstar) - Zbeta_block(tau, wl, 1, 1, Tstar)) / 6
