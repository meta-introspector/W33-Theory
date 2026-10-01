# Pass 11252 — time reversibility is decidable exactly; three qutrits decided

Producer: `analysis/w33_pass11252_exact_reversibility.py`
Certificate: `data/w33_pass11252_exact_reversibility.json`
Regression: `tests/test_w33_pass11250_11254.py`

**Problem.** Pass 11239 could only bound three-qutrit reversal fidelities from below.
* The anti-unitary three-qutrit Clifford group has ~3.3×10¹² elements.
* Generic search over it is blind: it reached 0.41 on an exactly reversible tick.
* Magnitude-only certificates miss the minimal violator (Pass 11228).

## The criterion (exact)

Expand U = 3⁻ⁿ Σ_p c(p) W(p) in symmetric Weyl operators W(x,z) = ω^{2xz} X^x Z^z. These satisfy:
* W(p)W(q) = ω^{⟨p,q⟩} W(p+q), W(p)† = W(−p), and W(x,z)* = W(x,−z);
* the Weil representation acts **without phases**: V_M W(p) V_M† = W(Mp).

Matching Weyl coefficients in V U* V† = λU† then gives:

> **U has a substrate time reversal ⟺ c(Lp) = μ·ω^{⟨b,Lp⟩}·c(p) for all p,
> for some anti-symplectic linear L, some b ∈ F₃²ⁿ and some |μ| = 1.**

So reversibility is a property of the Weyl coefficients alone.

## Algorithm

* Backtrack over the images of a basis.
* Prune every partial map on three conditions:
  * the anti-symplectic form;
  * |c(Lp)| = |c(p)| on the whole partial span;
  * solvability over F₃ of ψ(p)/ψ(p₀) = ω^{f(p−p₀)} with f linear, where ψ = c∘L / c.
* A map that reaches a leaf satisfies the criterion on all of F₃²ⁿ, so the first leaf proves reversibility.
* Exhausting the search proves violation. Magnitude classes use gap clustering, never rounding, so equal values are
  never split.

**Independent check.** For every map found, V_M is also built by the Weil twirl, and |tr(V U* V† U)| is evaluated over
all Weyl cosets. The two methods are asserted to agree.

## Results

**Validation against the exhaustive 51 840 × 81 search (two qutrits).**
* **204/204 agree:** 200 random words with 0–3 cubic gates, plus 4 named ticks.
* The named ticks are (T⊗T)SUM and (I⊗T)SUM(I⊗T²) (violating), and (I⊗T)SUM and (I⊗T)SUM(I⊗T) (reversible).

**Three qutrits, exact verdicts.**

| tick | verdict | Pass 11239 (generic search) |
|---|---|---|
| [(T⊗T)·SUM] ⊗ I | **violating** | ≥ 0.844 |
| (T⊗T⊗T)·SUM₁₂ | **violating** | ≥ 0.844 |
| (I⊗I⊗T)·SUM₂₃ | reversible | 1 (product reversal) |
| (T⊗T⊗T)·SUM₂₃·SUM₁₂ | **violating** | 0.237 |
| (T⊗I⊗T)·SUM₂₃·SUM₁₂ | **reversible** | 0.200: the generic search missed it |
| (I⊗I⊗T)·SUM₂₃·SUM₁₂ | reversible | — |
| (T⊗I⊗I)·SUM₁₂·SUM₂₃ | reversible | — |

* **Positive controls:** 6/6 Clifford conjugates C U C† of a reversible tick are found reversible.
* **Census** (random three-qutrit Clifford words with k cubic gates, 100 each): **10%, 94%, 100%, 100%** violate for
  k = 1, 2, 3, 4. None is undecided. The two-qutrit fractions were 8/65/93/98% (Pass 11235), so T-violation becomes
  generic faster with a third qutrit.
* Every violation at k ≥ 3 is certified at the root: not even a partial map survives.

**What this replaces.** The question "is this three-qutrit tick T-violating?" is now answered exactly, in seconds. The
best fidelity F_T is still a maximisation, which this criterion does not compute.
