# Pass 11372 — the Clifford+T walk irrep by irrep: exact moment laws, exact t-design rates, algebraic contraction rates, and no link to the 0.72

Producer: `analysis/w33_pass11372_irrep_contraction_rates.py` (main table; `--moments`; `--family`)
Certificate: `data/w33_pass11372_irrep_contraction_rates.json`
Regression: `tests/test_w33_pass11369_11373.py`

**Setting.** The one-qutrit walk is U_k = C_k T ⋯ C_1 T, with C_i uniform in the Clifford group.
* On an irrep λ of PU(3) its Fourier transform is P_λ λ(T), where P_λ is the Clifford average.
* ρ_λ is the largest modulus below 1 of P_λ λ(T) P_λ.
* Passes 11354/11359 computed ρ on tensor spaces, which mix irreps. Here ρ is resolved **irrep by irrep**, using
  Symᵖ V ⊗ Sym^q V̄ = (p,q) ⊕ (p−1,q−1) ⊕ ⋯: the moduli that are new at (p,q) belong to (p,q).
* It is matrix-free. The Pauli average is done exactly by combinatorics; then 24 coset representatives are averaged.
* Validation: it reproduces 3/8 on (3,3), 11/24 on (5,5) and 1/2 on (9,0), the irrep values behind Pass 11359.

## 1. Exact moment laws (theorem plus enumeration check)

* The steps are independent, so E π(U_k) = (Pπ(T))^k on any tensor representation.
* Taking traces: **E|tr U_k|^{2t} = Σ ν^k** over the nonzero spectrum of P_t π_t(T).

> **E|tr U_k|⁶ = 6 + (3/8)^k,  E|tr U_k|⁸ = 23 + 16·(3/8)^k.**

* The spectra are: t = 3, 1 (×6) and 3/8 (×1); t = 4, 1 (×23) and 3/8 (×16).
* 6 and 23 are the Haar values: the permutations of 3 and 4 letters with no increasing run longer than 3.
* 16 is the multiplicity of the irrep (3,3) in V^{⊗4} ⊗ V̄^{⊗4}.
* Checked against full enumeration of all coset-reduced words for k = 1, 2, 3: 6.375, 6.140625, 6.052734375 and 29,
  25.25, 23.84375, to 10⁻¹³.

## 2. Exact unitary t-design rates of the Clifford+T walk

* E π_t(U_k) converges to the Haar projector at rate **r_t = max ρ_λ over the irreps λ ⊂ V^{⊗t} ⊗ V̄^{⊗t}**. Which
  irreps occur comes from the mixed weights of Pass 11370.

| t | 1–2 | 3–4 | 5 | 6 | 7 | 8–10 |
|---|---|---|---|---|---|---|
| r_t | **0** (exact 2-design after one step) | **3/8** | **11/24** | **7/12** | **1/8 + √(7/32)** | **ρ₍₈,₈₎ = 0.72276…** |
| attained at | — | (3,3) | (5,5) | (8,2), (2,8) | (7,7) | (8,8) |

* The exact values:
  * ρ₍₇,₇₎ is a root of 64x² − 16x − 13.
  * **ρ₍₈,₈₎ is the largest root of 2592x³ − 1332x² − 585x + 140.** The new-eigenvalue polynomial of the (8,8) block
    is (4x − 1)(2592x³ − 1332x² − 585x + 140).
* **Independent cross-check.** Coherent states w = v^{⊗p} ⊗ ū^{⊗q} on Symᵖ ⊗ Sym^q V̄ share no code with the
  Pauli-orbit method. They give 3/8, 11/24 and **0.7227601** (with 0.59270717 and −0.3970744 also present) for
  (3,3), (5,5) and (8,8).
* At t = 12 the irrep (18,0) enters, with **ρ₍₁₈,₀₎ = the largest root of 864x³ − 984x² + 230x + 9 = 0.78092…**.
* All denominators divide a small multiple of 216 = |Cl|.

## 3. The table and the supremum

* All irreps with p + q ≤ 24 (109 irreps with p ≡ q mod 3), plus the (p,0) family up to p = 48 by coherent states.
  * Coherent states: Symᵖ is spanned by v^{⊗p}, with ⟨v^{⊗p}, u^{⊗p}⟩ = (v†u)^p, so P π(T) has matrix G⁻¹K.
  * Validated on (9,0), (12,0) and (18,0).
* **ρ is not monotone in the degree.**
  * Along (p,0): 0, 0, 1/2, 7/12, 0, 0.7809, 1/2, 37/108, 1/2, 0.5965, 89/216, 0.6673, 1/2, 0.6126, 0.6067, 0.7515
    (p = 3, …, 48).
  * Largest value found: **ρ₍₁₈,₀₎ = 0.78092…** The next largest are (14,8) and (8,14) at 0.7427, (8,8) at 0.7228, and (18,6) and (6,18) at 0.7007.
* **Theorem (by citation):** sup_λ ρ_λ < 1.
  * Benoist–de Saxcé (Invent. Math. 2016) give a spectral gap for algebraic measures generating a dense subgroup of a
    compact simple Lie group.
  * Apply it to μ²μ̌², whose support contains Cl and T Cl T⁻¹. The closure of ⟨Cl, T Cl T⁻¹⟩ has Lie algebra of
    dimension 8 (checked numerically), i.e. all of su(3).
  * Then ρ_λ ≤ (1 − δ)^{1/4} uniformly. The value of the supremum is open.

## 4. Correction: the 0.72 is not a representation-theoretic rate

* Passes 11354/11359 suggested that the reversible fraction's empirical decay, about 0.72 per gate (Pass 11312), is
  governed by high-degree representations.
* Two facts now undercut that:
  * the supremum of ρ is at least 0.781 > 0.72;
  * the reversible set is Haar-null, so its indicator is not an L² observable, and no Fourier rate theorem applies to
    it.
* ρ₍₈,₈₎ = 0.7228 sits strikingly close to the empirical rate. It is recorded as an **unexplained coincidence**, not a
  mechanism.

**Prior art.** Convergence of Clifford+T circuits to unitary designs is studied for multi-qubit systems (e.g.
Haferkamp et al. 2020, "Efficient unitary designs with a system-size independent number of non-Clifford gates").
We found no exact single-qutrit rates or moment laws. The search covered the corpus for "3/8", "11/24", "7/12",
"t-design rate" and "moment law", plus RESULTS_INDEX.
