# Pass 11249: every AME(10,3) state with a twisted 10-cycle symmetry is F₉-linear

Producer: `analysis/w33_pass11249_ame10_cyclic_search.py`
Certificate: `data/w33_pass11249_ame10_cyclic_search.json`
Regression: `tests/test_w33_pass11249_ame10_cyclic_search.py`

**Question.** Is every AME(10,3) stabilizer state F₉-linear (Glynn's arc)? Passes 11224 and 11229 answered by random
search. Here the search is exhaustive inside a structured family. The family is every ten-qutrit stabilizer state
invariant under a twisted cyclic shift σ: the 10 qutrits are permuted cyclically, with a local twist T ∈ SL(2,3)
applied at the wrap-around.

**Why this family.** Any automorphism that permutes the qutrits as a 10-cycle with local Cliffords attached is
LC-conjugate to such a σ: move every local part onto one coordinate, and the conjugacy class of T is what remains.
Glynn's state has 10-cycles in its permutation group PGL(2,9), so it must appear. It does, which is the positive
control.

**Method (exact over F₃).**
* For the semisimple twists T ∈ {I, −I, J} (J² = −1), σ is semisimple.
* F₃²⁰ splits into the components ker f(σ), one for each irreducible factor f of x^N − 1. The symplectic form pairs f
  with its reciprocal f*.
* Every invariant Lagrangian is a direct sum of K-subspaces, and all of them are enumerated.
* Each one is tested for:
  * minimum distance 6, over all 3¹⁰ codewords;
  * F₉-linearity: the algebra of local endomorphisms preserving the code is solved by linear algebra and searched for
    an element with A_i² = −1 on every qutrit.

| twist T | order of σ | invariant Lagrangians | AME(10,3) | F₉-linear | not F₉-linear |
|---|---|---|---|---|---|
| I | 10 | 1600 | **48** | **48** | **0** |
| −I | 20 | 336 | 0 | – | 0 |
| J | 40 | 8 | 0 | – | 0 |

**Result.**
* Every AME(10,3) stabilizer state with a 10-cycle symmetry whose twist has order 1, 2 or 4 is F₉-linear. By Pass
  11224 it is therefore Glynn's state.
* No other twist class admits an AME state at all.

**Not covered.**
* Twists of order 3 or 6, where σ is not semisimple in characteristic 3.
* States whose symmetry group has no 10-cycle.

So full uniqueness remains **open**, but a non-F₉ state, if it exists, cannot carry this symmetry. This is evidence
of a different kind from the random searches of Passes 11186, 11224 and 11229.
