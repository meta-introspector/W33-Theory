# Pass 11035 — cubic cocycle contracts to qutrit Weyl algebra and doubles to W33

Producer: `analysis/w33_pass11035_cubic_cocycle_weyl_w33_descent.py`
Certificate: `data/w33_pass11035_cubic_cocycle_weyl_w33_descent.json`
Regression: `tests/test_w33_pass11035_cubic_cocycle_weyl_w33_descent.py`

Motivated by the uploaded temporal-cohomology and anomaly PDFs, choose the alternating representative

ν(g,h,k) = ω^det(g,h,k),   g,h,k ∈ F₃³.

The normalized 3-cocycle identity is checked exhaustively on all 27⁴ quadruples.

Contracting any one coordinate axis leaves a bilinear 2-cocycle on the complementary F₃² plane. Its alternating commutator form has trivial radical, so its twisted group algebra is the full qutrit Weyl algebra M₃(C). All three coordinate contractions are nondegenerate.

Now double future and past Weyl planes with opposite orientation. Projectivizing F₃⁴ gives exactly 40 points. The maximal isotropic 2-planes give exactly 40 lines, four points per line and four lines per point: W(3,3).

The diagonal {(u,u)} is a Lagrangian Bell line. Exactly 27 Lagrangian lines are transverse to it. This makes the PDF’s “27 boundary trivializations” count exact inside W33.

The q=3 tetracode has eight nonzero weight-three words. They split exactly as four omitted coordinates times two opposite codewords, giving a precise 4×2 omission/chirality index set.

Boundary: the 4×2 count is an exact combinatorial bijection, not yet the missing intertwiner proving that tetracode sign equals the committed E8 tensor-sector cocycle chirality.
