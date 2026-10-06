# Pass 11512 — the cubed-phase test is NOT exact on PU(3): an explicit non-reversible unitary passes it (it stays exact on every Clifford+T word computed)

Producer: `analysis/w33_pass11512_cubed_phase_exhaustive.py` (main run; stage `strong`)
Certificate: `data/w33_pass11512_cubed_phase_exhaustive.json`
Regression: `tests/test_w33_pass11511_11515.py`

**Question (Passes 11500 and 11504).**
* On every Clifford+T word, the cubed-phase test decides reversibility: for some anti-symplectic L, |c(Lp)| = |c(p)| and
  (c(Lp)/c(p))³ is constant.
* The test forgets whether the residual phase f in c(Lp) = μ·ω^{f(p)}·c(p) is affine.
* A counterexample on PU(3) would be a unitary with a non-affine f that is not reversible.

## Exhaustive reduction

* **Symmetries.** Two reduce all 24 × 3⁹ pairs (L, f) without loss.
  * f matters modulo affine functions, because an affine part is absorbed by μ and the Weyl translation.
  * Clifford conjugation sends (L, f) to (SLS⁻¹, f∘S⁻¹) for S ∈ SL(2, 3), since c(p) ↦ c(S⁻¹p).
* **What remains:** **816 orbit representatives** with non-affine f, carrying 5442 eigenspaces of the monomial matrices
  P_{(L,f)}. Their dimensions are 4291 × 1, 664 × 2, 281 × 3, 153 × 4, 48 × 5 and 5 × 6.
* **The candidates.** Each unitary candidate is U = Σ_p c(p) W(p) with c in such an eigenspace.

## Two passes

1. **BFGS, 12 starts per eigenspace.**
   * Only 4 orbits contain a unitary, and all 4 are reversible.
   * Each was hit by about one start in twelve. A single eigenspace could therefore be missed with probability about
     (11/12)¹² ≈ 0.35, which motivated the second pass.
2. **Strong pass (stage `strong`).**
   * The 4291 one-dimensional eigenspaces are decided **exactly**: U is fixed up to scale, so it is unitary iff UU† is a
     multiple of 1.
   * The 1151 higher-dimensional ones were searched by Levenberg–Marquardt from 200 starts each, with hit rates recorded.

## Result

**The strong pass overturns the first.**

| | BFGS, 12 starts | strong pass (exact for d = 1; LM × 200 for d ≥ 2) |
|---|---|---|
| unitary-containing eigenspaces | 4 | **1212** |
| … of dimension 1, 2, 3, 4, 5, 6 | — | 414, 355, 237, 153, 48, 5 |
| … reversible | 4 | 1210 |
| … **not reversible** | 0 | **2** |
| minimum LM hit rate (d ≥ 2) | — | 0.495 |

* BFGS missed about 99.7% of the unitary-containing eigenspaces. Levenberg–Marquardt hits each of them from at least half of
  its starts.
* **The two non-reversible unitaries** both lie in 6-dimensional eigenspaces of L₀ with non-affine
  f = (0,0,1,0,0,1,2,2,0) and f = (0,0,2,0,0,2,1,1,0) (labels (a,b) in lexicographic order).
* **One is written into the certificate** (stage `counter`, key `counterexample`), with every check:

| check | value |
|---|---|
| unitarity error | 7·10⁻¹⁶ |
| cubed-phase test | **passes** |
| max Clifford overlap (Theorem 1) | **2.985 < 3** |
| Pass 11252's exact Weyl criterion (independent decider) | **not reversible** |

> **The cubed-phase test is a necessary condition for reversibility on PU(3), but not a sufficient one.** The affine
> condition that cubing discards matters for some unitaries.

**Reading.**
* On one-qutrit Clifford+T words the test is exact through depth 9 (Pass 11500). The counterexamples are therefore not
  words, or at least not words of depth ≤ 9.
* The words' coefficients lie in a cyclotomic ring. That is the natural suspect for why the discarded condition never
  bites on words, but this is not proved.
* Pass 11504's reading ("4 unitaries, all reversible: evidence for exactness") is **withdrawn**. Its optimizer could not
  see the unitaries.
* **Lesson.** A negative search is only as good as its solver's hit rate. Measure the hit rate before reading an empty
  result.

**Prior art.** The criterion is Pass 11252's, and the word-level observation is Pass 11500's. The symmetry reduction and
both searches are new here.
