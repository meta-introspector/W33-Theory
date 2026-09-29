# Pass 11164 — the three-qutrit tetracode gate lives over F₉: K = P² + i·𝟙

Producer: `analysis/w33_pass11164_three_qutrit_tetracode.py`
Scan: `analysis/w33_pass11164_scan_f9_unitary.py`, frozen as `data/w33_pass11164_f9_unitary.json`
Regression: `tests/test_w33_pass11164_three_qutrit_tetracode.py`

**No F₃ analogue.** The tetracode gate K⊗I is built from F₃ scalars. A three-qutrit analogue would need a 3×3 matrix
over F₃ whose entries and 2×2 minors are all nonzero.
* None exists: normalise each row's first entry to 1, and the three rows would need pairwise distinct second ratios
  in {±1}.
* Equivalently, there is no [6,3,4] MDS code over F₃.

**The F₉ analogue.** Work over F₉ = F₃[i], with conjugation x ↦ x³ and ⟨u,v⟩ = Σ ū_i v_i.
* The imaginary part of ⟨·,·⟩ is the three-qutrit symplectic form.
* So a unitary matrix over F₉ embeds, via x = a+bi ↦ [[a,−b],[b,a]], as a symplectic 6×6 matrix with F₉-linear blocks.
* The block determinant is the norm, so the gate is perfect iff every entry is nonzero.

**Enumeration.**

| | order | all entries nonzero (perfect, F₉-linear) |
|---|---|---|
| U(2, F₉) | 96 | 64 |
| U(3, F₉) | 24 192 (= \|U(3,3)\|) | **12 288** |

* Every column has norms (1, 1, 2). This is forced: three nonzero norms must sum to 1. It is the (1, 1, −1)
  determinant pattern of Pass 11158.
* 24 of the perfect gates are circulant.

**Closed form.**

    K = P² + i·𝟙,   K†K = (P − i𝟙)(P² + i𝟙) = I + 𝟙² = I,

because 𝟙² = 3·𝟙 = 0 in characteristic 3. The all-ones matrix is nilpotent exactly because **the number of parties
equals the characteristic**. The entries are i and 1+i, all nonzero, so the gate is perfect. Pass 11165 checks its
information pattern.

**Codes (cross-reference).** The graph {(u, Ku)} of the two-qutrit gate K = [[1,1],[1,−1]] is (a, b, a+b, a−b). This
is the tetracode of Pass 10946's clock code (a, b, a+b, b−a) with the last coordinate's sign flipped, which is exactly
the orientation sign discussed there. It is also the q = 3 rung [4,2,3]₃ of Pass 10970's projective-line tower. The
three-qutrit K = P² + i𝟙 has all nine 2×2 minors nonzero (det K = 1), so its graph is a **[6,3,4]₉ MDS code**, and
that is why F₉ is needed. The [6,3,4]₅ of Pass 10970 is over F₅ and carries no qutrit symplectic structure.

**Scope.** Embedding U(n,q²) in Sp(2n,q) is standard, and so is building AME states from MDS codes over extension
fields. New here: the three-qutrit perfect gates inside W(3,3)'s Clifford group, the F₉-linear count, and the closed
form K.
