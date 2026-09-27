# Pass 11032 — the minimum-glue law survives exactly through q = 17

Producer: `analysis/w33_pass11032_prime_clock_q17_reduction.py`
Certificate: `data/w33_pass11032_prime_clock_q17_reduction.json`
Regression: `tests/test_w33_pass11032_prime_clock_q17_reduction.py`

The previous exact minima were 2,4,6,10,12 for q = 3,5,7,11,13. The observed law is minimum glue norm = q−1.

At q = 17 the full code has 17⁹−1 = 118,587,876,496 nonzero words. We avoid that enumeration using the projective symmetry of the extended Reed–Solomon clock code.

Its homogeneous degree is m=(q−1)/2. Projective normalization multiplies values by μᵐ, which by Euler’s criterion is ±1. The A(q−1) discriminant norm is invariant under sign, so PGL₂(q) acts by norm-preserving signed permutations.

Any word with norm below q−1 must, by evenness, have norm at most q−3. Since each nonzero symbol contributes at least (q−1)/q, such a word has at most q−3 nonzero coordinates and therefore at least four zeros. Three-transitivity moves three zeros to 0, 1, and infinity.

The counterexample search reduces to f(x)=x(x−1)g(x) with deg g≤(q−7)/2: only q^((q−5)/2) words, a q³ reduction.

For q=17 the complete reduced search checks 17⁶−1 = 24,137,568 words. Its exact minimum is 16 = q−1.

So the law is now exact through q=17. This is not yet an all-primes proof; the next exact target is q=19.
