# Pass 11046 — the existing 81D matter carrier becomes the regular K81 module

Producer: `analysis/w33_pass11046_full_k81_dressed_equivariant_compiler.py`
Certificate: `data/w33_pass11046_full_k81_dressed_equivariant_compiler.json`

Keep the carrier

`C9_multiplicity tensor C3_internal tensor C3_external`.

Replace the trivial H27 action on the multiplicity factor by the latent `A9` action, while leaving the old internal Schrödinger action and external `C3` register in place.

Then

`(A9 tensor V_omega) tensor Reg(C3) = Reg(H27) tensor Reg(C3) = Reg(K)`.

So the scheduler and dressed execution modules are now isomorphic on the same 81-dimensional vector space. No ancilla enlargement is needed.

The explicit full compiler is `T81 = T_H tensor F3`, and it intertwines the H27 `X,Z,C` generators plus the external shift at both split primes.

Most importantly, the old compiler's three 27D blocks now have a structural meaning:

`S1 = 3chi tensor V` (compatible),
`S2 = V tensor V -> 3Vbar`,
`L = Vbar tensor V -> sum_9 chi`.

Thus the previous 27 + 27 + 27 permutation was the Clebsch–Gordan shadow of the minimal commutant dressing.
