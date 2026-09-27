# Pass 11027 — prime clock lattices through q = 13

Producer: `analysis/w33_pass11027_prime_clock_lattice_census.py`
Certificate: `data/w33_pass11027_prime_clock_lattice_census.json`
Regression: `tests/test_w33_pass11027_prime_clock_lattice_census.py`

Pass 10970 proved the all-prime root firewall for the projective-line
self-dual MDS clock code. This pass extends the exact Euclidean-glue census.

Exact exhaustive minima:

| q | 3 | 5 | 7 | 11 | 13 |
|---|---:|---:|---:|---:|---:|
| minimum glue norm | 2 | 4 | 6 | 10 | 12 |
| minimum codewords | 8 | 72 | 240 | 552 | 756 |

Every tested rung satisfies minimum glue norm = q − 1. For q = 13 the deep
scan checks all 13⁷ − 1 = 62,748,516 nonzero codewords. The equality q − 1
is certified only through q = 13 and is not promoted to an all-primes theorem.

The genuine all-prime theorem remains the MDS inequality. With

d = (q + 3)/2

and minimum nonzero glue norm at least d(q − 1)/q, comparison with root norm
two reduces to

(q − 3)(q + 1) > 0.
Therefore every odd prime q ≥ 5 is strictly above root norm, while q = 3
hits root norm exactly. Hence q = 3 is the unique odd-prime rung capable of
creating new roots through nontrivial glue.

Root shells:

- q = 3: E₈ with 240 roots; tetracode glue adds 216 roots to A₂⁴.
- q ≥ 5: no new glue roots, so the root system remains A(q−1)^(q+1).
- The base root count is q(q² − 1), with root/rank ratio q.

This isolates a real exceptional discontinuity: q = 3 has 240/8 = 30 roots
per rank, while the higher-prime tower remains on the base-root law.

The q = 13 exhaustive scan is replayable with `--deep-q13`. Ordinary CI uses
the frozen exact result rather than enumerating 62.7 million words every run.
No continuum or physical-dimensional interpretation is inferred.
