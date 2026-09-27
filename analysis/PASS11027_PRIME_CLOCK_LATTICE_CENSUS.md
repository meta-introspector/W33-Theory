# Pass 11027 — prime clock lattices through q=13

Producer: `analysis/w33_pass11027_prime_clock_lattice_census.py`
Certificate: `data/w33_pass11027_prime_clock_lattice_census.json`
Regression: `tests/test_w33_pass11027_prime_clock_lattice_census.py`

Pass 10970 proved the all-prime root firewall for the projective-line
self-dual MDS clock code. This pass pushes the exact Euclidean-glue census
through two new rungs.

The exhaustive minima are

[
egin{array}{c|ccccc}
q&3&5&7&11&13\ hline
min|v_{
m glue}|^2&2&4&6&10&12
end{array}
]

so every tested rung satisfies (min|v_{
m glue}|^2=q-1).
The q=13 computation checks all (13^7-1=62,748,516) nonzero codewords.
That last equality is deliberately **not** promoted to an all-primes theorem.

The genuine all-prime theorem remains the MDS inequality. With

[
d=
rac{q+3}{2},
qquad
|v_{
m glue}|^2ge d
rac{q-1}{q},
]

comparison with root norm two reduces to

[

rac{(q+3)(q-1)}{2q}>2
iff
(q-3)(q+1)>0.
]

Thus every odd prime (qge5) is strictly above root norm, while (q=3)
hits it exactly.
Consequently the first root shell is

[
q=3:quad E_8,;240	ext{ roots},
]

because tetracode glue adds 216 roots to the 24 roots of (A_2^4). For every
odd prime (qge5),

[
oxed{Phi(L_q)=A_{q-1}^{,q+1}},
]

with exactly

[
q(q^2-1)
]

roots and root/rank ratio (q).

So the exceptional discontinuity is now very concrete: q=3 has
(240/8=30) roots per rank, whereas the higher-prime tower stays on the
base-root law. The deep q=13 scan is replayable with `--deep-q13`; default
regression uses its frozen exact result so ordinary CI does not enumerate
62.7 million words on every run.
