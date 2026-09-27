# Pass 11034 — exhaustive GL₂(3) search for an independent proton-protection symmetry

Producer: `analysis/w33_pass11034_independent_flavour_symmetry_search.py`
Certificate: `data/w33_pass11034_independent_flavour_symmetry_search.json`
Regression: `tests/test_w33_pass11034_independent_flavour_symmetry_search.py`

We enumerate every semisimple complex GL₂(3) representation of total dimension three: 12 family-representation types. Higgs fields may carry either one-dimensional character 1 or det. The FI-forced singlet is required to be trivial so its VEV preserves the candidate symmetry.

We then require nondegenerate, full-rank invariant Yukawa pairings for QUHᵤ, QDH_d, and LEH_d and test UDD, QLD, and LLE.

There are 576 full-rank Yukawa-compatible assignments.

If every family is one of the two irreducible 3D triplets, there are 16 assignments and every one allows all three RPV tensors with multiplicity one. Nonabelian triplet flavour therefore does not solve the problem.

Across the entire 576-assignment search, exactly four assignments forbid all three RPV tensors. Every survivor uses only 3·1 and 3·det family representations. Hence every survivor factors through the generation-blind determinant quotient

det : GL₂(3) → F₃* ≅ Z₂.

Importantly det(−I)=+1, so this Z₂ is independent of the clock’s central −I parity. It is an algebraic escape route, but it is an extra parity rather than a flavour texture.

The next decisive test is objectwise: determine whether one of the four determinant-charge patterns is realized by an exact geometric/space-group automorphism of the committed orbifold vacua and remains unbroken by every required singlet VEV.
