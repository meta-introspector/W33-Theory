# Pass 11190 — every perfect five-qutrit Clifford gate is C10 or cross+perm: a Klein-quadric proof via the signs of AME(10,3) states

Producer: `analysis/w33_pass11190_ame10_sign_structure.py`
Frozen: `data/w33_pass11190_ame10_sign_structure.json`
Regression: `tests/test_w33_pass11190_ame10_sign_structure.py`

**The conjecture (Pass 11186).** Perfect five-qutrit Clifford gates realise only two of the five allowed orientation
patterns, 'C10' and 'cross+perm'. The other three, 'C4+C6', 'double cross' and 'all', never occur. **This pass proves
it**, for every such gate, not only the F₉-linear ones or the 71 found by tabu search.

**Setup.** A perfect gate's Choi state is an AME(10,3) stabilizer state: a Lagrangian L ⊂ W₁ ⊕ … ⊕ W₁₀ with
L ∩ W_T = 0 for every 5-set T.

**1. Signs.** For a 4-set K, V_K = L ∩ W_{K^c} is 2-dimensional, and each projection to a party p ∉ K is an isomorphism.
Take s_p = ω_p(a_p, b_p) = ±1 for a basis (a, b). This splits the 6-set K^c into sign classes. Isotropy gives Σ s_p ≡ 0
(mod 3), so K^c is either **uniform** or split **3 + 3**.

**2. Dictionary.** From the Choi Lagrangian {(u, S T u)}:

    det S[o, i] = −1   ⇔   i and o lie in the same sign class of the 6-set {i} ∪ O

This was checked on 12 600 blocks. It means the pattern of a gate is read off six of the 210 sign partitions.

**3. The Klein quadric.** For a 3-set A, V_A = L ∩ W_{A^c} is 4-dimensional. The seven forms β_r = π_r^*ω_r are points
of the Klein quadric Q⁺(5,3) in Λ²V_A^*:
* their radicals V_{A∪r} are pairwise transverse (AME), so β_r ∧ β_{r'} ≠ 0;
* they sum to zero (isotropy).

With G_{rr'} = β_r ∧ β_{r'}:

    σ_{A^c∖r'}(r, r'') = G_{rr'} G_{r''r'},     every row of G sums to 0 (mod 3).

So the **sign graph Γ_A** (the −1 entries of G, up to complement) has all degrees in {0, 3, 6}. There are 1052 labelled
graphs in four classes. Realisability as a Gram matrix in the hyperbolic Klein space rules out two of them:

| class | Gram rank | disc | realisable |
|---|---|---|---|
| 7K₁ | 6 | 1 (elliptic) | no |
| K₃,₃ + K₁ | 6 | 1 (elliptic) | no |
| K₄ + 3K₁ | 6 | 2 (hyperbolic) | yes |
| prism + K₁ | 4 | — | yes |

**4. Duality.** Take disjoint 4-sets K, K′ leaving {i, j}. Isotropy between V_K and V_{K′} forces the two transfer maps
W_i → W_j to be negative contragredients. Hence σ_{K^c}(i, j) = σ_{K′^c}(i, j).

**5. SAT.** Each of the 210 six-sets carries one of 11 sign structures. The constraints are a Klein-realisable Γ on every
7-set and duality on 1575 pairs: 56 910 variables and 488 580 clauses. Both CaDiCaL and Glucose4 give:
* **Q1**: a K₄ + 3K₁ sign graph on {0,…,6} is **UNSAT**. The constraints are S₁₀-invariant, so every Γ_A is prism + K₁.
  Its isolated vertex is the only uniform 4-set through A. Hence **the uniform 4-sets of every AME(10,3) stabilizer state
  form a Steiner system S(3,4,10)**, the inversive plane of order 3.
* **Q2**: any pattern other than C10 or cross+perm at inputs {0,…,4} is **UNSAT**. The controls are SAT: C10 and
  cross+perm are each realisable. **Theorem: every perfect five-qutrit Clifford gate has pattern C10 or cross+perm.**

**Corollary: the counts are forced.** cross+perm ⇔ some I∖{i} is uniform ⇔ I contains a block of the S(3,4,10). A 5-set
holds at most one block, since two blocks share at most two points. So every AME(10,3) state gives cross+perm for exactly
the 180 input sets containing a block, and C10 for the other 72. The 2 : 5 ratio is forced in both sources: the F₉ census
(6 291 456 : 15 728 640) and the 71 tabu graphs (5112 : 12 780).

**Positive control.** All 71 AME(10,3) graphs of Pass 11186 satisfy every constraint: zero-sum, prism + K₁ on all 120
seven-sets, duality on all 1575 pairs, and S(3,4,10).

**Scope and prior art.**
* The proof is computer-assisted. Every constraint is derived above from the axioms (Lagrangian + AME), so the SAT
  instances are relaxations of the true problem, and UNSAT is conclusive. The solver verdicts agree across two solvers;
  no DRAT certificate is stored.
* Klein correspondence and Witt-type classification are standard.
* The link between AME states and combinatorial designs is known in the orthogonal-array sense (Goyeneche et al.,
  arXiv:1506.08857; Helwig et al., arXiv:1306.2879).
* The sign structure, the forced S(3,4,10) and the pattern theorem were not found in the corpus (`RESULTS_INDEX.md`, grep
  for S(3,4,10)) or in a web search. The rediscovery guard's 'clifford+steiner' candidates are Steiner trihedral pairs
  and S(5,6,12) (`analysis/PASS20260916_execute_all5_intertwiner_magic_closure.md`, `W33_FOR_EVERYONE.tex`), which are
  unrelated objects.
