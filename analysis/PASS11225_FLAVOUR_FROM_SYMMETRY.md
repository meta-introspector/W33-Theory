# Pass 11225: can the geometry's own symmetry fix the mixing angles? An exhaustive answer

Producer: `analysis/w33_pass11225_flavour_from_symmetry.py`
GAP export: `analysis/gap/w33_pass11225_three_dim_images.g` → `data/w33_pass11225_three_dim_images.txt.gz`
Certificate: `data/w33_pass11225_flavour_from_symmetry.json`
Regression: `tests/test_w33_pass11225_flavour_from_symmetry.py`

**Question.**
* The paper's scorecard lists "masses, mixing angles, coupling constants" as OPEN.
* The older manuscripts give the angles as SRG-parameter ratios: sin²θ₁₂ = 4/13, sin²θ₂₃ = 7/13, sin²θ₁₃ = 2/91. The
  paper rightly refuses to count those.
* The principled version of the question is the residual-symmetry method of flavour physics:
  * the three generations carry a 3-dimensional representation of a finite group;
  * the charged-lepton and neutrino mass terms keep subgroups G_e and G_ν;
  * then group theory alone fixes U_PMNS = U_e† U_ν.
* For Majorana neutrinos G_ν is a Klein group (full prediction) or a single Z₂ (one predicted column). For Dirac
  neutrinos it is any abelian group with nondegenerate spectrum. A Z₂ on the charged-lepton side predicts one row.

**Scope.**
* Every 3-dimensional irreducible representation of every subgroup, up to conjugacy, of:
  * W(E₆) = Aut(W33), with 350 subgroup classes;
  * Sp(4,3), with 162 subgroup classes.
* That is 606 representations and **34 distinct image groups in U(3)**. Among them are A₄, S₄, A₅, Δ(27), Δ(54),
  Σ(36×3), Σ(72×3), Δ(81) and Δ(162)-type groups, and the single-qutrit Clifford group Σ(216×3) (order 648, with its
  C₂ extension).
* Data, as quoted from NuFIT 6.0, normal ordering, 3σ: sin²θ₁₂ ∈ [0.275, 0.345], sin²θ₁₃ ∈ [0.02030, 0.02388],
  sin²θ₂₃ ∈ [0.430, 0.596], δ free. Generations may be assigned freely (row and column permutations).
* CKM: |V_us| ∈ [0.2230, 0.2271].

## Results

1. **No complete pattern survives.** No pair of residual symmetries in any subgroup fixes all three lepton angles inside
   the data. This holds for Majorana neutrinos (Klein groups) and for Dirac neutrinos (any nondegenerate abelian
   group). **Mixing angles cannot come from the geometry's finite symmetry alone.**
2. **No Cabibbo angle.** No pair of residual symmetries gives a quark mixing with |V_us| in range and the other mixing
   small.
3. **The legacy numerology is not a symmetry pattern.** No residual pattern has (sin²θ₁₂, sin²θ₁₃, sin²θ₂₃) =
   (4/13, 2/91, 7/13).
4. **Five single columns survive** (one predicted column plus a free angle and phase). One row survives as well: the A₅
   row (0.127, 0.436, 0.436).

| surviving column | group it comes from | as PMNS column | predicted sin²θ₁₂ | forced cos δ |
|---|---|---|---|---|
| **TM₁: (2/3, 1/6, 1/6)** | S₄, the quotient of the **line stabiliser 3³:S₄** | 1st | **0.317–0.320** | −0.33 … +0.45 |
| (2/3)·sin²(kπ/9), k = 4, 1, 2: (0.6466, 0.0780, 0.2755) | **order-9 elements of the single-qutrit Clifford group** | 1st | 0.338–0.340 | −1 … −0.95 (δ ≈ π) |
| golden: (0.6545, 0.0955, 0.25) | A₅ (from 2⁴:A₅) | 1st | 0.329–0.332 | ±(0.66 … 1) |
| golden: (0.2764, 0.3618, 0.3618) | A₅ | 2nd | 0.282–0.283 | broad |
| TM₂: (1/3, 1/3, 1/3) | A₄, the quotient of the point stabiliser | 2nd | 0.340–0.341 | broad |

* Only TM₁ lies within about 1σ of the global fit, sin²θ₁₂ = 0.307 (+0.012/−0.011). The others sit 2–3σ away and will
  be decided by a precise θ₁₂ measurement.
* TM₁ comes from the S₄ by which the stabiliser of a W33 **line** acts. A line is a maximal set of commuting
  two-qutrit Paulis, i.e. a measurement context.
* The Clifford column is new as far as the corpus goes. It is exactly ((1 − cos 2π/9)/3, (1 − cos 4π/9)/3,
  (1 − cos 8π/9)/3). Its sharp corollary is that δ is near π (CP nearly conserved).

**Correction during the pass.** GAP's Dixon representations are not unitary in general.
* A first run used them directly and reported six "viable" complete patterns and Cabibbo patterns for the qutrit
  Clifford group. They were not even doubly stochastic: rows summed to 1.5, an artefact of non-orthogonal eigenbases.
* Every image is now unitarised by the square root of its invariant Hermitian form. Unitarity of every matrix,
  orthonormality of every eigenbasis and unistochasticity of every pattern are asserted.
* Several first-run rows and columns, all from non-unitary images, disappeared with the fix.

## Reading for the scorecard

* **Mixing angles from the finite symmetry alone: NO-GO.** The geometry's symmetry hosts the classic flavour groups as
  stabiliser quotients, and with them the known trimaximal and golden-ratio partial patterns. It cannot fix the
  complete mixing matrix, and it has no Cabibbo angle.
* **What it does supply is one column.** The best of them, TM₁ from a measurement context, gives a falsifiable pair:
  sin²θ₁₂ ≈ 0.318 and cos δ ∈ [−0.33, 0.45].
* Masses and the remaining angle need dynamics: vacuum alignment, which subgroups are kept, and RG running. That remains
  OPEN.

## Prior art

* The residual-symmetry method and the TM₁, TM₂ and golden-ratio patterns are standard. Representative references:
  Lam (2008); Altarelli–Feruglio (A₄); Albright–Rodejohann (TM₁/TM₂); Everett–Stuart and Ding (A₅); Fonseca–Grimus
  (classification of residual mixing in finite subgroups of U(3)). These were not checked against network sources from
  here.
* A legacy script in the corpus, `scripts/w33_democratic_mixing.py`, proves exactly democratic mixing (all 1/3) for one
  Steinberg-sector construction. It is a different construction (81-dimensional projectors) and is not superseded.
* The rediscovery guard's Cabibbo/PMNS-with-Clifford/qutrit hits do not overlap with this pass:
  * `exploration/CKM_FROM_27_LINES.py` derives CKM from the 27-line intersection structure, not from residual
    symmetries;
  * `W33_FOR_EVERYONE.tex` and the 2026-09-23 frontier note mention CKM/PMNS phases only as open targets.
* New here:
  * the exhaustive application to every subgroup of Aut(W33) and Sp(4,3);
  * the NO-GO for complete patterns and for the Cabibbo angle;
  * the identification of TM₁ with the line stabiliser;
  * the (2/3)sin²(kπ/9) column of the qutrit Clifford group;
  * the check that the legacy numerology is not a symmetry pattern.
