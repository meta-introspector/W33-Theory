# Pass 11504 — is the cubed-phase test exact beyond words? A search of PU(3) finds no counterexample

Producer: `analysis/w33_pass11504_cubed_phase_test.py`
Certificate: `data/w33_pass11504_cubed_phase_test.json`
Regression: `tests/test_w33_pass11498_11505.py`

## Question

* **What Pass 11500 found.** On every one-qutrit Clifford+T word to depth 9, the following condition coincides with
  reversibility: for some anti-symplectic L of F₃², |c(Lp)| = |c(p)| and (c(Lp)/c(p))³ is constant.
* **What is exact.** Pass 11252's exact criterion needs the stronger c(Lp) = μ ω^{f(p)} c(p) with f **affine**. Cubing
  forgets whether f is affine.
* **What a counterexample would be.** A unitary obeying the relation with a **non-affine** f that is *not* reversible.

## Method

* The relation is **linear** in the coefficient vector c. For fixed (L, f), c must lie in the μ-eigenspace of the 9×9
  monomial matrix P_{(L,f)}, defined by (P c)(Lp) = ω^{f(p)} c(p).
* So every candidate is U = Σ_p c(p) W(p) with c in such an eigenspace.
* **Search.** 4000 random draws of (L, f) were made, and the non-affine ones kept.
  * In every eigenspace, ‖UU† − 1‖² was minimised by BFGS from 3 starts.
  * Each unitary found was tested with Theorem 1's exact overlap decider.

## Result

| quantity | value |
|---|---|
| non-affine (L, f) drawn | 3 999 |
| eigenspaces searched | 26 958 |
| unitaries found in them | **4** |
| … reversible (overlap = 3) | **4** |
| … not reversible (counterexample) | **0** |

## Reading

* **Unitarity is almost incompatible with a non-affine phase relation.** Only a handful of the ~27 000 eigenspaces
  contain a unitary at all.
* **Every unitary that does exist is reversible**, by some other anti-symplectic L with an affine relation.
* No counterexample was found. This is evidence, not a proof, that the cubed test is exact on all of PU(3).
  * BFGS from 3 starts can miss unitaries in an eigenspace.
  * The draws cover a fraction of the 24 × (3⁹ − 27) pairs (L, f).
* **Why exactness would matter.** The test involves only the cubes c(p)³ and the moduli |c(p)|². If it is exact on PU(3),
  then time reversal of a qutrit tick is decided by **cubic data of its Weyl coefficients**, with no phase-linearity
  check.

**Prior art.** The criterion is Pass 11252's, and the observation on words is Pass 11500's. No external source is known.

## Correction (Pass 11512)

* **The search here was blind, and its reading is withdrawn.** BFGS from 3 starts found unitaries in 4 eigenspaces.
  Levenberg–Marquardt, with the one-dimensional eigenspaces decided exactly, finds **1212** unitary-containing eigenspaces
  among the 816 symmetry classes.
* **Two of them hold non-reversible unitaries that pass the cubed-phase test.** One is verified by Theorem 1 (overlap
  2.985) and by Pass 11252's criterion.
* So the test is **not** exact on PU(3). The word-level coincidence of Pass 11500 stands.
