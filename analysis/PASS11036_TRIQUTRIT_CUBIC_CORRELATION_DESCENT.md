# Pass 11036 — exact CCZ temporal echo and a correction to the “missing 12” story

Producer: `analysis/w33_pass11036_triqutrit_cubic_correlation_descent.py`
Certificate: `data/w33_pass11036_triqutrit_cubic_correlation_descent.json`
Regression: `tests/test_w33_pass11036_triqutrit_cubic_correlation_descent.py`

For the qutrit phase polynomial f(a,b,c)=abc, successive forward finite differences give exactly

abc → bc → c → 1.

Exponentiation is the gate ladder

CCZ → controlled two-qutrit phase → local Z-like phase → central ωI.

Independently, choosing one reference ray in each qutrit gives the correlation-order dimension split

27 = 1 + 6 + 12 + 8.

The 12 is three pairwise blocks of dimension four.

The important correction is linear-algebraic. Model the tripartite sector as U₁⊗U₂⊗U₃ with each Uᵢ two-dimensional. Each fixed slant contraction is surjective onto its own four-dimensional pair block. But the three fixed contractions combined map an 8D source into the 12D direct sum with rank only 7 and a one-dimensional common kernel.

Therefore “8→12→6→1” cannot literally be a surjective linear chain. It is a branching correlation hierarchy. A diagonal weld may activate all three pairwise blocks, but CCZ alone cannot prove that it restores all 12 independent pairwise degrees of freedom.

This is a useful firewall: the uploaded PDFs’ stronger missing-12 claim now has a precise executable condition for the still-open E6 intertwiner.
