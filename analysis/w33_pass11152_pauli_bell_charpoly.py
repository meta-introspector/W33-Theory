#!/usr/bin/env python3
"""Pass 11152: the best spooky correlation available with W(3,3)'s own measurements is an exact algebraic number -- the
Pauli-only CGLMP Bell operator has characteristic polynomial (x^3 - 6x - 2)(x^3 - 3x + 1)^2.

Pass 11149 found that with Pauli measurements only, the best two-qutrit state reaches 2.6016791.  The optimal Pauli
strategy is the simplest one: Alice and Bob both measure in the Z and X eigenbases (natural labelling).  Here the Bell
operator of that strategy is built exactly: the spectral projectors of Z and X are (I + w^-k P + w^-2k P^2)/3 with entries
in Q(w), so every coefficient of its characteristic polynomial is a rational number with denominator dividing 9^9; a
60-digit Faddeev-LeVerrier computation therefore identifies the coefficients exactly (deviation < 1e-58):
    char(x) = x^9 - 12x^7 + 45x^5 - 6x^4 - 57x^3 + 18x^2 + 6x - 2 = (x^3 - 6x - 2)(x^3 - 3x + 1)^2.
Hence:
  * the Pauli-only maximum is the largest root of x^3 - 6x - 2: 2 sqrt2 cos((1/3) arccos(1/(2 sqrt2))) = 2.6016791 --
    exceeding the local bound 2 with nothing but W(3,3) measurements, provided the STATE carries magic (the maximising
    state is not a stabilizer state: mana 0.392, Pass 11149);
  * the rest of the spectrum is 2cos(2 pi k/9), k = 1, 2, 4 (the roots of x^3 - 3x + 1, each twice): ninth roots of unity
    in a qutrit (order-3) problem.  Why ninth roots appear is NOT explained here (open).
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_pass11152_pauli_bell_charpoly.json"


def bell_charpoly(dps=60):
    mp.mp.dps = dps
    w = mp.mpc(-0.5, mp.sqrt(3) / 2)
    mat = lambda f: mp.matrix([[f(i, j) for j in range(3)] for i in range(3)])
    X = mat(lambda i, j: 1 if i == (j + 1) % 3 else 0)
    Z = mat(lambda i, j: w ** i if i == j else 0)
    proj = lambda P: [(mp.eye(3) + P * (w ** (-k)) + P * P * (w ** (-2 * k))) / 3 for k in range(3)]
    PZ, PX = proj(Z), proj(X)
    kron = lambda A, B: mp.matrix([[A[i // 3, j // 3] * B[i % 3, j % 3] for j in range(9)] for i in range(9)])

    def term(Am, Bm, plus, minus):
        t = mp.zeros(9, 9)
        for a in range(3):
            for b in range(3):
                s = (1 if plus(a, b) else 0) - (1 if minus(a, b) else 0)
                if s:
                    t += s * kron(Am[a], Bm[b])
        return t
    B = (term(PZ, PZ, lambda a, b: a == b, lambda a, b: a == (b - 1) % 3)
         + term(PX, PZ, lambda a, b: b == (a + 1) % 3, lambda a, b: b == a)
         + term(PX, PX, lambda a, b: a == b, lambda a, b: a == (b - 1) % 3)
         + term(PZ, PX, lambda a, b: b == a, lambda a, b: b == (a - 1) % 3))
    n, M, c, I = 9, mp.zeros(9, 9), [mp.mpf(1)], mp.eye(9)
    for k in range(1, n + 1):
        M = I if k == 1 else B * M + c[-1] * I
        c.append(-sum((B * M)[i, i] for i in range(n)) / k)
    return c, B


def summarize():
    c, B = bell_charpoly()
    x = sp.symbols('x')
    target = [int(v) for v in sp.Poly(sp.expand((x ** 3 - 6 * x - 2) * (x ** 3 - 3 * x + 1) ** 2), x).all_coeffs()]
    dev = float(max(abs(mp.re(a) - b) for a, b in zip(c, target)))
    imag = float(max(abs(mp.im(a)) for a in c))
    top = float(max(mp.re(e) for e in mp.eig(B)[0]))
    res = dict(pass_id=11152, charpoly_coeffs=target, deviation=dev, max_imag=imag, top_eigenvalue=top,
               cubic_root=float(max(sp.Poly(x ** 3 - 6 * x - 2).nroots())),
               trig=float(2 * mp.sqrt(2) * mp.cos(mp.acos(1 / (2 * mp.sqrt(2))) / 3)),
               ninth_roots=[float(2 * mp.cos(2 * mp.pi * k / 9)) for k in (1, 2, 4)])
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
