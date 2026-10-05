# Pass 11499 — the "good" rules: a frame-level reversibility theorem for M z₁ = −z₁ and M²z₁ = −z₁, valid for every n

Producer: `analysis/w33_pass11499_good_rules.py` (main run; stage `union`)
Certificate: `data/w33_pass11499_good_rules.json`
Regression: `tests/test_w33_pass11498_11505.py`

**Background.** Pass 11373 found four universal one-gate rules:
* the fixed axis (Mz₁ = z₁) and M²z₁ = z₁ are bad;
* Mz₁ = −z₁ and M²z₁ = −z₁ are good, meaning every frame is reversible.

Pass 11498 proved the bad half (the magic-axis law) for every n. This pass treats the good half, in the same variables:
Pass 11350's (S)/(F) at k = 0, with A = QJ.

## Theorem (any n)

Let Mz₁ = −z₁ or M²z₁ = −z₁, and let Q be a symplectic solution of (S) at k = 0, with A = QJ.

> Every frame a with ω((A⁻¹ − I)x, a) = 0 for all x ∈ ker(M − I) is **reversible**.
> In particular, if A fixes ker(M − I) pointwise, **every frame is reversible**.

**Proof.** At k = 0 we have t₀ = 0, and r_Q = 0 because Qz₁ = z₁ (Pass 11498, step 3). So (F) reads
(I − M⁻¹)e = −(M⁻¹ + A)a, with e_{x₁} = 0.
1. **Solvability.** Im(I − M⁻¹) = Im(M − I) = ker(M − I)^⊥. For x ∈ ker(M − I),
   ω(x, (M⁻¹ + A)a) = ω(x, a) − ω(A⁻¹x, a) = −ω((A⁻¹ − I)x, a), which vanishes by hypothesis. So a solution e exists.
2. **The constraint e_{x₁} = ω(e, z₁) = 0 holds automatically.**
   * Choose y with (I − M)y = z₁: y = −z₁ if Mz₁ = −z₁, and y = −(z₁ + Mz₁) if M²z₁ = −z₁.
   * Then ω(e, z₁) = ω((I − M⁻¹)e, y) = −ω((M⁻¹ + A)a, y).
   * Use Az₁ = −z₁ and AM = M⁻¹A (from (S): A = MAM), so AMz₁ = Mz₁ when M²z₁ = −z₁.
   * The two terms then cancel, so ω(e, z₁) = 0, and every solution already satisfies the constraint. ∎

The pointwise-fixing hypothesis is **linear in Q**: QJx = x on a basis of ker(M − I). So its existence is decided by
Gaussian elimination plus a symplecticity filter.

## Coverage

**n = 2 (all 2592 good classes): every one is covered.**

| route | M z₁ = −z₁ | M²z₁ = −z₁ |
|---|---|---|
| a single A fixing ker(M − I) | 620 | 1890 |
| union of the frame-level criteria over all k = 0 solutions (dim ker(M − I) = 2) | 27 | 54 |
| M = −I: ker(M − I) = 0, Q = I solves (S) | 1 | — |
| **total** | **648/648** | **1944/1944** |

* **Soundness.** In every case the certified frames are reversible in the decider's verdicts.
* **Exact cell shares.** Pr(Mz₁ = −z₁) = 1/(3^{2n} − 1) and Pr(M²z₁ = −z₁) = 3/(3^{2n} − 1). These give exactly the 648
  and 1944 classes at n = 2.

**n = 3 (Pass 11373's orbits in the good cells).**
* The single-map hypothesis holds on **96.8%** of the Mz₁ = −z₁ mass and **95.8%** of the M²z₁ = −z₁ mass.
* Every orbit where it holds has all frames reversible: 77 and 69 orbits, all decided.
* The remainder is 2.8% and 3.2% where the hypothesis fails, plus 0.5% and 1.0% beyond the enumeration cap.

**Together with Pass 11498, all four of Pass 11373's rules now have proofs.**
* The two bad rules are proved for every n.
* The two good rules are proved at n = 2 on every class, and at n = 3 on about 96% of the mass.
* What is still missing for the good rules at all n is a proof that the frame-level criteria always cover F₃^{2n}. The
  n = 2 failures of the single-map form all have dim ker(M − I) = 2, and unions cover them there.

**Prior art.** The rules are Pass 11373's. The proof uses Pass 11350's normal form and Pass 11498's step 3. No external
source is known.
