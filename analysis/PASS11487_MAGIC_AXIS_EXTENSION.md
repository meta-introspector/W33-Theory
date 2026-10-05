# Pass 11487 — the magic-axis law when the magic axis lies inside Im(M − I): a second theorem, sound on every case

Producer: `analysis/w33_pass11487_magic_axis_extension.py` (stages `n2`, `n3`, `n4`)
Certificate: `data/w33_pass11487_magic_axis_extension.json`
Regression: `tests/test_w33_pass11486_11492.py`

## Setting

* **The law.** A one-gate tick is U = W(a) V_M T₁, with z₁ the magic axis. The magic-axis law F′ (Pass 11420) says: if
  M²z₁ = z₁, then every frame a with ω(a, v) ≠ 0 violates, where v = z₁ + Mz₁.
* **What Pass 11457 proved.** It proved F′ when R_M = Im(M − I) is nondegenerate. In general it proved F′ for frames
  a ⊥ Rad, where Rad = R_M ∩ R_M^⊥, provided z₁ ∉ R_M.
* **The gap.** When z₁ ∈ R_M, v lies in Rad itself. So "a ⊥ Rad" excludes every F′ frame, and Pass 11457 says nothing.

## Lemma (level hyperplanes)

If z₁ ∈ R_M, all three Weyl shifts of T₁ land on the same coset. The coefficient moduli of G = V_M T₁ are then supported
on R_M and equal

|R_M|^(−1/2) · |Σ_j τ_j ω^{j s(p) + q j²}|,  with s affine on R_M.

* **q ≠ 0:** the three values 0.778, 0.508 and 1.462 (× |R_M|^(−1/2)) are distinct. The level sets are therefore the
  three cosets of a hyperplane K ⊂ R_M.
* **q = 0:** the modulus is constant, so the class is *magnitude-blind*.
* **Checked per class.** The support equals R_M, and the level sets are three cosets of one index-3 subspace.

## Theorem

In the hyperplane case, every frame a with ω(a, v) ≠ 0 and ω(a, Rad ∩ K) = 0 is violating.

**Proof.**
1. Suppose U is reversible. Since m_U(p) = m_G(p − a), Pass 11457's Lemma 2 gives an anti-symplectic L preserving m_U.
2. The three level sets a + X_k carry distinct values, so L fixes each of them.
3. Hence L(a + R_M) = a + R_M. So L R_M = R_M and t := La − a ∈ R_M.
4. Also t + Lr lies in the coset of r, so Lr − r ∈ K for all r ∈ R_M.
5. L preserves R_M, R_M^⊥ and Rad. Since v ∈ Rad, r₀ := Lv − v ∈ Rad ∩ K.
6. Then −ω(a, v) = ω(La, Lv) = ω(a + t, v + r₀) = ω(a, v) + ω(a, r₀) = ω(a, v).
7. So 2ω(a, v) = 0, a contradiction. ∎

**Corollary.** If Rad ∩ K = 0, the class satisfies F′ completely.

## Verification

**n = 2, exhaustive (all F′-setting classes with z₁ ∈ R_M).**

| | count |
|---|---|
| classes with z₁ ∈ R_M | 920 |
| level hyperplane (q ≠ 0) | 752 |
| — theorem sound (every certified frame violates in the exhaustive verdicts) | **752 / 752** |
| — F′ fully proved (Rad ∩ K = 0) | 704 |
| — F′ partly proved (dim Rad ∩ K = 1) | 48 |
| magnitude-blind (q = 0) | 168 (the open case, Pass 11486) |

* **Combined with Pass 11457:**
  * Pass 11457 fully proves F′ on 352 of its 376 classes.
  * This pass fully proves it on 704 of 920.
  * Together, F′ is **proved** on 1056 of the 1296 F′-setting classes of n = 2.

**n = 3** (orbits of Pass 11373 in the two bad cells).
* **Soundness.** The extension is sound on all 76 hyperplane orbits that the decider reaches. A further 40 are undecided
  there.
* **Coverage by orbit mass:**

| outcome | mass |
|---|---|
| fully proved by Pass 11457 | 21.2% |
| **fully proved by the extension** | **42.4%** |
| extension partial (Rad ∩ K ≠ 0) | 7.3% |
| Pass 11457 partial | 3.7% |
| magnitude-blind (q = 0) | 25.5% |

**n = 4** (the 317 Stab(z₁) classes the decider could not reach, Pass 11433). There are no verdicts here; the certificate
is the theorem itself.
* Pass 11457 already proved 44 of these classes (25.4% of their mass).
* **The extension proves 70 more (13.7%).**
* Left over:
  * 68 magnitude-blind classes (31.2%);
  * 90 extension-partial classes (19.8%);
  * 45 Pass 11457 partial classes (9.9%).
* Pass 11486 applies its phase theorems to these.

**What remains.** Magnitudes certify nothing on the q = 0 classes, and only part of each partial class. Pass 11486 closes
both at n = 2 with phases: every one of the 1296 classes gets a uniform proof.

**Prior art.** This extends Pass 11457's argument from the cosets R_M + jz₁ to the hyperplane cosets of K. Lemma 2 and the
support lemma are Pass 11457's. No external source is known for the magic-axis law itself.
