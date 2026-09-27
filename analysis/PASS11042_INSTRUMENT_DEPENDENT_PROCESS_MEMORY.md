# Pass 11042 — memory visibility depends on the intervention family

Producer: `analysis/w33_pass11042_instrument_dependent_process_memory.py`
Certificate: `data/w33_pass11042_instrument_dependent_process_memory.json`

Two explicit qutrit one-slot process combs are constructed. Both are positive and satisfy the causal trace rule
`Tr_F Upsilon = I_O tensor I_I/3`.

The classical hidden-memory comb has zero normalized process negativity. The coherent memory comb has negativity 1.

For the classical control, a hidden trit H is copied to past P and future F, while the exposed middle system carries `|H>`.

A computational-basis instrument reveals H. Conditioning on its outcome screens past from future exactly:
`I(P:F|M_Z)=0`.

A Fourier-basis instrument produces an outcome independent of H. Past and future remain perfectly correlated:
`I(P:F|M_X)=log2(3)`.

Both instruments have normalized uniform outcomes and the same forward causal order.

This is an exact finite example of instrument-dependent Markov order: memory can disappear under one intervention family and remain visible under another without any retrocausal influence.
