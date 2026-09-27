# Pass 11041 — the phase-weld 12 is a dressed pairwise sector

Producer: `analysis/w33_pass11041_phase_weld_pairwise12_dressing.py`
Certificate: `data/w33_pass11041_phase_weld_pairwise12_dressing.json`

The diagonal phase weld has a 12-dimensional intersection with the 27-dimensional operator-compatible S1 sector. The uploaded note suggested identifying this directly with the triqutrit pairwise-correlation 12.

In the frozen S1 coefficient factorization `(t,r,i) in F3^3`, choosing reference value 0 gives the exact correlation-order split
`27 = 1 + 6 + 12 + 8`.

At split primes 103 and 109, for both diagonal slopes `c+p` and `c-p`, the weld intersection K12 has projection ranks
`C0:1, C1:6, C2:12, C3:8`,
but intersects every pure sector trivially:
`K12 ∩ Ck = 0` for k=0,1,2,3.

Since `dim K12 = dim C2 = 12` and the C2 projection has rank 12, K12 is the graph of a unique map from the full pairwise C2 sector into C0⊕C1⊕C3.

So the pairwise 12 is fully activated, but never in isolation: the diagonal weld dresses it simultaneously with scalar, one-body, and irreducibly tripartite components.
