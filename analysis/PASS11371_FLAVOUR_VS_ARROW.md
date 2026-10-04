# Pass 11371 — the arrow of time is the Altland–Zirnbauer class of the tick; the flavour Jarlskog is blind to it by theorem and exclusive with it at one gate

Producer: `analysis/w33_pass11371_flavour_vs_arrow.py`
Certificate: `data/w33_pass11371_flavour_vs_arrow.json`
Regression: `tests/test_w33_pass11369_11373.py`

**The map.** This answers the open point of the Codex flavour audit (decision 035b65c4: "an explicit full-rank flavour
map is still needed").
* A one-qutrit tick U is a 3×3 unitary. Read it as the mixing matrix between the clock basis, where the magic gate T
  is diagonal (the "mass basis"), and its image:
  H_u = diag(m_u), H_d = U diag(m_d) U†.
* Then (checked to 10⁻¹⁵, sign included):

> det[H_u, H_d] = 2i · J_CKM(U) · Π_{i<j} (m_{u,i} − m_{u,j})(m_{d,i} − m_{d,j}),  J_CKM(U) = Im(U₁₁U₂₂Ū₁₂Ū₂₁).

## Theorem 1 (exact, elementary): flavour CP cannot see the arrow

* **The parities.** J_CKM(Uᵀ) = J_CKM(U), while J_CKM(Ū) = J_CKM(U†) = −J_CKM(U).
* **Monomial invariance.** J_CKM(CUC†) = J_CKM(U) for every monomial C. Rows and columns move by the same
  permutation, and sgn² = 1.
* **Consequences.**
  * **Flavour CP is even under the substrate time reversal U ↦ Uᵀ** (Pass 11355). No function of J_CKM can witness
    the arrow.
  * If U has an antiunitary symmetry (Ū ~ CUC†) or a unitary reversal (U† ~ CUC†) realised by a **Borel**
    Clifford, then J_CKM(U) = 0. The Borel Cliffords are those preserving the clock basis: 54 of 216.

## Lemma 2 (exact): the reversal types form a Klein four-group lattice

* The four operations 1, T: U ↦ Uᵀ, K: U ↦ Ū and I: U ↦ U† form a Klein four-group.
* The set of operations X for which X(U) is Clifford-conjugate to U, up to phase, is a **subgroup**. This uses that
  Cl is closed under conjugation and that Cᵀ = C̄⁻¹.
* So any two of {reversible, antiunitary-symmetric, unitarily reversible} imply the third. Possible types: none, T,
  K, I, TKI.
* **Observed** on every word with up to 3 cubic gates: only **none, T, I, TKI** occur. A tick with an antiunitary
  symmetry but no reversal (type K) never appears. This is computed, not proved.
* Type K does occur in PU(3). A generic real rotation R has R̄ = R, but R⁻¹ is not Clifford-conjugate to R. So its
  absence is a property of the Clifford+T words, not of the group.

## The Klein four-group IS the Altland–Zirnbauer tenfold way (odd-dimensional part)

Write U = e^{−iHt}. The three channels are the three AZ symmetries of H:

| channel | operator | on U | on H | AZ symmetry |
|---|---|---|---|---|
| ᵀ (the arrow) | antiunitary Θ = C̄K | ΘUΘ⁻¹ ∝ U⁻¹ | ΘHΘ⁻¹ = +H | time reversal **T** |
| K | antiunitary C†K | ΘUΘ⁻¹ ∝ U | ΘHΘ⁻¹ = −H | particle–hole **C** |
| I | unitary S | SUS⁻¹ ∝ U⁻¹ | SHS⁻¹ = −H | chiral **S = T·C** |

* Lemma 2, "any two imply the third", is the AZ relation S = TC.
* In odd dimension a scalar square must be **+1**: det(V V̄) = |det V|² > 0, while (−1)³ < 0.
* So the five symmetry types are exactly the five odd-dimensional AZ classes:
  none = **A**, T = **AI**, I = **AIII**, K = **D**, TKI = **BDI**.
* **The substrate's arrow of time is the AZ symmetry class of its Floquet operator.** Reversible ticks are AI or BDI;
  violating ticks are A or AIII.

**Checked on every word with 1 and 2 cubic gates** (`--az`, certificate key `altland_zirnbauer`):

| gates | A | AI | AIII | D | BDI |
|---|---|---|---|---|---|
| 1 | 0 | 180 | **18** (all the violators) | 0 | 18 |
| 2 | 2052 | 2808 | 54 | **0** | 270 |

* Every Clifford-realised T or C whose square is a scalar has square **+1**: 6804 T and 810 C at two gates.
* Every reversible word has at least one reversal with T² = +1.
* Some realisations square to a non-scalar Clifford instead. These are ticks with an extra unitary symmetry, e.g. 162
  BDI words at two gates whose particle–hole operators all square non-trivially.
* **All chiral (AIII and BDI) ticks show exact chiral spectral pairing** e^{iθ} ↔ λe^{−iθ}, with 0 failures.
* **Consequences:**
  * every single-magic-gate violator is **chiral class AIII**: it breaks the arrow while keeping a chiral symmetry;
  * **but only for one qutrit.** For two qutrits (`--az2`; 3000 uniformly sampled one-gate ticks, exact verdicts
    from the Pass 11350 decider), **151 of 274 violators have no spectral reflection symmetry at all**. They are class
    A, chiral under no unitary whatsoever, since pairing is necessary for chirality. Violators are still about ten
    times more often spectrally paired than reversible ticks: 123/274 = 45%, against 127/2726 = 4.7%;
  * the level statistics of Pass 11353 (reversible COE, violating CUE) is the AI/A dichotomy;
  * class D (particle–hole without time reversal) never occurs among Clifford+T words. In PU(3) it does occur, e.g.
    real rotations.

**Prior art.** The tenfold way is Altland–Zirnbauer (1997). Floquet unitaries are classified by the same symmetry
classes (e.g. Roy–Harper 2017). In this repository:
* Pass 11353 already noted that Θ² = +1 is forced in odd dimension.
* Pass 11210 owns the phase-space Kramers statement: θ² = −1, meaning T² is the parity operator, needs an even
  number of qutrits and is exactly the 36 E₆ reflections for two qutrits.
* Pass 5632 used tenfold-way vocabulary ("class-D-like") for a different object, the stabiliser carrier.

The identification of the substrate arrow of a Floquet tick with its AZ class, and the census, are new here.

## Computed: the joint census (exact Clifford search, all coset-reduced words of Pass 11312)

**The Cliffords themselves.**
* The 54 Borel Cliffords have J_CKM = 0.
* **All 162 others are complex Hadamard matrices with the maximal Jarlskog |J_CKM| = 1/(6√3)**, the trimaximal value.

**Violating fraction by flavour-CP class.**

| cubic gates | words | violating when J_CKM = 0 | when \|J_CKM\| = J_max | when J_CKM is another value |
|---|---|---|---|---|
| 1 | 216 | 18/54 | **0/162** | — |
| 2 | 5184 | 1998/3240 = 0.62 | **108/1944 = 0.056** | — |
| 3 | 124 416 | 23 814/48 600 = 0.49 | **648/23 328 = 0.028** | 43 740/52 488 = 0.83 |

* **At one gate the arrow and flavour CP are mutually exclusive.**
  * All 18 arrow-breaking ticks have J_CKM = 0: their Clifford is a diagonal shear (Pass 11266), hence Borel.
  * All 162 ticks with maximal flavour CP are reversible.
* **Deeper words.** Maximal flavour CP strongly suppresses the arrow: 2.8% violate at three gates, against 49% with no
  flavour CP and 83% for the generic Jarlskog values that first appear at three gates.

**Reading.**
* The substrate's T-violation and flavour CP violation sit in **different channels** of the Klein four-group.
  * The arrow is the ᵀ-channel.
  * The CKM Jarlskog is odd in the K- and I-channels and even in ᵀ.
* So the substrate arrow is **not** the origin of CKM-type CP violation under this map, and CKM-type CP does not
  imply the arrow.
* The two are correlated on Clifford+T words: maximally mixing ticks are almost always reversible. That correlation
  is computed, not proved.
* **Scope.** This holds for one specific, natural flavour map (the clock basis as the mass basis). Pass 11360's
  J₆ values are a different, Clifford-invariant quantity.
* **Relation to the Codex flavour bridge.** Passes 11361 and 11374 (extended to full ray orbits and native E₆ fields
  in 11379 and 11384) build a different map: three supplied rays with
  supplied positive weights, and a CP-even selector that picks φ = π/2 or 3π/2. That map uses an ansatz; the map here
  uses the tick itself. The two are not compared here. Both use the same basis-independent Jarlskog commutator.

**Prior art.** J_max = 1/(6√3), attained by the trimaximal (Fourier) matrix, is classical (Jarlskog 1985). The
Klein-four structure of reversing symmetries is standard in the theory of reversing symmetry groups (Lamb–Roberts,
Physica D 112, 1998; Baake–Roberts). The corpus has no earlier statement that J_CKM is transpose-even, and no flavour
census of substrate ticks. The search covered "transpose-even", "J_CKM(U^T)", "reversing symmetry" and
"strongly real".
