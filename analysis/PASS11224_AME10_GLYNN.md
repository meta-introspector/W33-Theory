# Pass 11224: the ten-qutrit AME state is Glynn's arc, unique up to local Cliffords

Producer: `analysis/w33_pass11224_ame10_f9_structure.py`
Certificate: `data/w33_pass11224_ame10_f9_structure.json`
Regression: `tests/test_w33_pass11224_ame10_f9_structure.py`

**Question.** How many AME(10,3) stabilizer states are there up to local Cliffords (LC) and relabelling of the qutrits?
Equivalently, how many perfect five-qutrit Clifford gates are there? This was the fifth next step of the session.

**Background (master's passes, cited rather than repeated).**
* Pass 11175 counted the F₉-linear perfect five-qutrit gates exhaustively, in one fixed F₉ structure: 2 642 411 520.
* Pass 11186 found 71 AME(10,3) graph states by tabu search.
* Pass 11190 proved that the sign structure is unique, an S(3,4,10) system with symmetry PGL(2,9).
* Pass 11214 identified that system with the Baer sublines of PG(1,9).

## Method

**F₉-linearity up to LC.**
* A local F₉ structure on a qutrit is an order-8 element M_v of GL(2,3), acting as multiplication by a primitive
  ξ ∈ F₉. There are 12 such elements: 6 of trace 1 and 6 of trace 2.
* Every five qutrits form an information set, so C = {(x, Ex)}. C is invariant under ⊕M_v exactly when E·M_S·E⁻¹ is
  block diagonal with order-8 blocks.
* This is a batch test over 12⁵ choices. Allowing all 12 elements per qutrit tests F₉-linearity up to LC.

**Reed–Solomon or Glynn.**
* An F₉-linear AME(10,3) code is a [10,5,6]₉ MDS code, i.e. a 10-arc in PG(4,9).
* Glynn (Discrete Math. 59 (1986) 43–51) showed there are two such arcs: the normal rational curve (Reed–Solomon / GRS
  codes) and his non-classical arc.
* They are told apart by the dimension of the Schur square C⋆C: 2k − 1 = 9 for GRS, and 10 for Glynn's code. Both
  values are checked on explicit codes.

## Results

1. **All 71 census states are F₉-linear up to LC.** Each has exactly 4 F₉ structures: ±ξ and the Frobenius conjugate,
   all of one trace.
2. **None is Reed–Solomon.** Each has F₉-dimension 5 and Schur square of dimension 10. **All 71 are Glynn's code.**
3. **No Reed–Solomon AME(10,3) exists.** No GRS [10,5]₉ code is Hermitian self-dual for any coordinate weighting:
   0 of the 1024 sign patterns. Glynn's code is Hermitian self-dual for exactly the uniform weights (±1).
4. **Automorphisms.**
   * Glynn's arc has 360 projective and 360 Frobenius-twisted stabilising maps in PΓL(5,9).
   * In all 720, the induced coordinate scalings have constant norm, so each lifts to 4 local-Clifford automorphisms.
     Hence **|Aut_LC| = 2880**.
   * On the 10 qutrits these act as a sharply 3-transitive group of order 720, with element orders 1, 2, 3, 4, 5, 8 and
     10: **PGL(2,9)**. The linear part is PSL(2,9) ≅ A₆. This is exactly the sign-structure symmetry of Pass 11190, and
     PGL(2,9) rather than PΓL(2,9), as Pass 11214 found.
5. **One orbit, and the count.**

       24¹⁰ · 10! / 2880 = 79 888 260 016 373 760 = 6¹⁰ · 2 642 411 520 / 2.

   * The right-hand side is master's exhaustive F₉ count, times the 6¹⁰ local F₉ structures of one trace, divided by
     the 2 structures each code has in that class.
   * The Glynn state's orbit under LC and relabelling therefore exhausts every F₉-linear AME(10,3) state. **The F₉-type
     AME(10,3) stabilizer state is unique.**
   * There are exactly 79 888 260 016 373 760 perfect five-qutrit Clifford gates of F₉ type.

## Open

Whether an AME(10,3) stabilizer state that is not F₉-linear exists. None appeared in 300 tabu restarts. If none
exists, the counts above are the complete answer.

## Prior art

* **Pass 11170** (`analysis/PASS11170_PERFECT_GATE_LADDER.md`) records that AME(10,3), i.e. [[10,0,6]]₃, comes from
  Grassl–Gulliver's circulant Hermitian self-dual [10,5,6]₉ code (Grassl–Rötteler arXiv:1502.05267; Danielsen
  arXiv:1106.2428). Its explicit circulant perfect gates are instances of it.
* Result 3 above implies that **Grassl–Gulliver's code is necessarily Glynn's code**: every Hermitian self-dual
  [10,5,6]₉ code is a non-GRS MDS code, hence Glynn's by Glynn's 1986 classification.
  * This identification may already be in the coding-theory literature. That could not be checked from here (no
    network access), and it is not claimed as new outside this corpus.
* Glynn's arc, MDS codes of length q+1 and Schur-product tests are classical.
* New in the corpus:
  * the identification with Glynn's arc ("Glynn" had no hits before this pass);
  * the non-existence of Hermitian self-dual GRS codes of this length (an exhaustive check here);
  * |Aut_LC| = 2880 with permutation image PGL(2,9);
  * uniqueness of the F₉-type state up to LC, by orbit–stabiliser against master's exhaustive count;
  * the exact count of F₉-type perfect five-qutrit Clifford gates.
* The rediscovery guard's "Baer" compounds point to unrelated files: the Baer sublines here are master's Pass 11214.
