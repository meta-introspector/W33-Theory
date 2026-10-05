# Pass 11486 — the magic-axis law has a uniform proof on all 1296 classes of two qutrits: phases close what magnitudes cannot, because every reversing map fixes v

Producer: `analysis/w33_pass11486_blind_classes_phases.py` (stages `n2`, `n3`)
Certificate: `data/w33_pass11486_blind_classes_phases.json`
Regression: `tests/test_w33_pass11486_11492.py`

## The gap this closes

* Passes 11457 and 11487 prove the magic-axis law F′ from the **moduli** of the Weyl coefficients.
* When z₁ ∈ R_M = Im(M − I) and the twisted moduli |Σ_j τ_j ω^{js+qj²}| are constant (q = 0), magnitudes certify
  nothing. These are the 168 *magnitude-blind* classes of n = 2 (Pass 11487).
* This pass uses the **phases**.

## Lemma 1 (pure phases)

* For q = 0, Σ_j τ_j ω^{js} is a diagonal entry of T, namely ζ^{s³} for s ∈ {0, 1, −1} (ζ = e^{2πi/9}).
* So on R_M:

  c_G(r) = const · ω^{Q(r)} · ζ^{s(r)³},  with Q quadratic and s = ℓ + κ affine.

* **Checked on every blind class:** the exponent function e = arg c / (2π/9) admits e ≡ 3Q + s³ (mod 9).

## Lemma 2 (the exact criterion becomes F₃ algebra)

* Pass 11252's criterion is c_U(Lp) = μ ω^{⟨b,Lp⟩} c_U(p). Here U = W(a)G, L is anti-symplectic with L R_M = R_M, and
  t = La − a ∈ R_M.
* The criterion becomes e(Lr + t) − e(r) ∈ 3·Aff + const.
* **Mod 3** (s³ ≡ s), this gives **(i)** ℓ(Lr) = ℓ(r).
* **Mod 9**, with c₀ = ℓ(t) and (s + c₀)³ − s³ = 3c₀s² + 3c₀²s + c₀³, it gives
  **(ii)** B(Lr, Lr′) = B(r, r′) − 2c₀ ℓ(r) ℓ(r′), where B is the polar form of Q.
* **Checked at n = 2:** the (L, t) satisfying the exact criterion reproduce the exhaustive verdicts on all 168 blind classes.

## Theorem

If ℓ(v) = 0 and ℓ = μ·B(v, ·) on R_M with μ ≠ 0, then every frame a with ω(a, v) ≠ 0 and a ⊥ (Rad ∩ rad B|_{R_M})
is violating.

**Proof.**
1. Put r = v in (ii). Since ℓ(v) = 0, this gives B(Lv, Lr′) = B(v, r′).
2. By (i), B(v, Lr′) = μ⁻¹ℓ(Lr′) = μ⁻¹ℓ(r′) = B(v, r′).
3. Subtracting, B(Lv − v, ·) vanishes on L R_M = R_M. So r₀ := Lv − v ∈ rad B|_{R_M}.
4. Also r₀ ∈ Rad, because L preserves R_M and R_M^⊥.
5. Then −ω(a, v) = ω(La, Lv) = ω(a + t, v + r₀) = ω(a, v) + ω(a, r₀) = ω(a, v).
6. So 2ω(a, v) = 0, a contradiction. ∎

## A second phase theorem (q ≠ 0): Pass 11487's 48 partial classes

* For q ≠ 0, the distinct moduli force s(Lr + t) = s(r) exactly, so c₀ = 0.
* (ii) then reads B(Lr, Lr′) = B(r, r′), with no ℓ(v) term.
* The same proof gives: if ℓ = μB(v, ·) on R_M, every frame with ω(a, v) ≠ 0 and a ⊥ (Rad ∩ rad B|_{R_M}) violates.
* **Checked on all 48 partial classes:**
  * c_G = const·ω^Q·τ̃_q(s), with Q quadratic and s affine;
  * ℓ ∝ B(v, ·);
  * B|_{R_M} nondegenerate.
  * **The theorem is sound and fully proves F′ on 48/48.**

## A third theorem (transvections): Pass 11457's 24 partial classes

* **Statement.** Suppose z₁ ∉ R_M and R_M = span(r₁), so that M is a transvection x ↦ x + c·ω(x, r₁)·r₁. Then every frame
  with ω(a, v) ≠ 0 violates.
* **Proof.**
  1. M²z₁ = z₁ forces ω(z₁, r₁) = 0. Hence Mz₁ = z₁ and v = −z₁.
  2. A reversing L fixes the three cosets jz₁ + R_M, since the |τ_j| are distinct. So Lz₁ = z₁ + κr₁.
  3. Because ω(r₁, z₁) = 0, the coefficient phase on each coset is ω^{q₁λ²}, with q₁ = χ(r₁, r₁) the Wall form. The Wall
     form is nondegenerate, so q₁ ≠ 0.
  4. The exact criterion allows only affine phase changes. The j·λ cross term 2q₁κελ must therefore vanish, so κ = 0.
  5. Hence Lv = v, and the usual step gives 2ω(a, v) = 0. ∎
* **Checked on all 24:** each is a transvection with Mz₁ = z₁ and q₁ ≠ 0, and F′ holds (sound 24/24).

## Results

**n = 2 (exhaustive).** On all 168 magnitude-blind classes:
* Lemma 1's decomposition holds;
* ℓ(v) = 0 and ℓ ∝ B(v, ·);
* B|_{R_M} is nondegenerate (GF(3) rank via `DomainMatrix`);
* the theorem is sound and fully proves F′;
* Lemma 2's criterion reproduces the exhaustive verdicts;
* **every admissible (L, t) fixes v**.

**The magic-axis law at n = 2 now has a uniform proof on every class.**

| theorem | hypothesis family | classes |
|---|---|---|
| Pass 11457 (moduli) | z₁ ∉ R_M, Rad = 0 | 352 |
| Pass 11487 (moduli, level hyperplane) | z₁ ∈ R_M, q ≠ 0, Rad ∩ K = 0 | 704 |
| this pass, phases q ≠ 0 | z₁ ∈ R_M, q ≠ 0, Rad ∩ K ≠ 0 | 48 |
| this pass, phases q = 0 | z₁ ∈ R_M, q = 0 (magnitude-blind) | 168 |
| this pass, transvection | z₁ ∉ R_M, dim R_M = 1 | 24 |
| **total** | | **1296 = all F′-setting classes** |

**Scope.**
* Each class is covered by a theorem whose hypotheses are F₃ linear-algebra facts, checked class by class. Reversibility
  is never decided in the proof.
* Soundness was nevertheless cross-checked against the exhaustive verdicts on every class.
* This is a uniform *explanation* of F′ at n = 2, which was already known to hold by exhaustion. It is not yet a
  single n-independent proof.

**n = 3** (the blind orbits among Pass 11373's orbits in the two bad cells).
* All 33 blind orbits admit the decomposition e = 3Q + s³, and all have ℓ(v) = 0.
* ℓ ∝ B(v, ·) holds on 24 of them.
* Those 24 are **fully proved**: 49.1% of the blind mass.
* The theorem is sound on all 21 of them that the decider reaches.
* On the other 9 orbits (50.9% of the mass) the hypothesis fails.

**n = 4** (the 317 Stab(z₁) classes the decider could not reach in Pass 11433). There are no verdicts here; the theorems,
with hypotheses checked per class, are the certificate.

| outcome | classes | share of the unreached mass |
|---|---|---|
| fully proved by Pass 11457 | 44 | 25.4% |
| fully proved by Pass 11487 | 70 | 13.7% |
| **blind (q = 0): fully proved here** | **68** | **31.2%** |
| **partial hyperplane (q ≠ 0): fully proved here** | **90** | **19.8%** |
| transvection: fully proved here | 2 | ~0 |
| open: z₁ ∉ R_M, Rad ≠ 0, dim R_M > 1 | 43 | 9.9% |

* **Every blind and every partial-hyperplane class at n = 4 satisfies the hypotheses.** This contrasts with n = 3, where
  9 blind orbits fail them.
* **The fixed-axis law at n = 4 is now:**
  * 97.60% of Stab(z₁) decided (Pass 11456);
  * **2.16% proved** (up from 0.61%);
  * **0.24% open** (down from 1.79%).
* **The open residue** is the non-transvection case z₁ ∉ R_M with Rad ≠ 0.
  * The coset-phase conditions there are B∘L = B on R_M and 2B̃(Lu − u, Lr) = ℓ′(r) − ℓ′(Lr), with u = Mz₁ and
    ℓ′(r) = −½ω(r, u).
  * They are derived but not yet shown to force Lv = v.

* **All theorems together at n = 3** (`n3_all_theorems`). Of the F′-setting orbit mass, **79.7% is fully proved**:

| outcome | orbits | mass |
|---|---|---|
| Pass 11457 | 38 | 21.2% |
| Pass 11487 | 76 | 42.4% |
| phases, q = 0 | 24 | 12.5% |
| phases, q ≠ 0 | 28 | 3.6% |
| transvection | 2 | ~0 |
| open: blind, ℓ ∝ B(v, ·) fails | 9 | 13.0% |
| open: q ≠ 0, ℓ ∝ B(v, ·) fails | 12 | 3.7% |
| open: z₁ ∉ R_M, Rad ≠ 0, dim R_M > 1 | 18 | 3.7% |

  The theorems are sound on all 102 of these orbits that the decider reaches.
* The decomposition and ℓ(v) = 0 survive at n = 3.
* The proportionality ℓ ∝ B(v, ·) does **not** always hold, and B|_{R_M} can be degenerate. The theorem is stated with
  these hypotheses rather than for all classes.
* **Open:** a hypothesis-free proof of Lv = v from Lemma 2. The n = 2 data say Lv = v always, but the general mechanism
  is not identified.

**Prior art.**
* Lemma 2 rests on Pass 11252's exact criterion.
* No external source is known.
