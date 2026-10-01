# Pass 11230: what TM₁ predicts, and whether a flavon potential selects it

Producer: `analysis/w33_pass11230_tm1_alignment.py`
Certificate: `data/w33_pass11230_tm1_alignment.json`
Regression: `tests/test_w33_pass11230_tm1_alignment.py`

**Starting point.** Pass 11225 found that among all residual-symmetry patterns of Aut(W33) and Sp(4,3), the only
column within about 1σ of the data is TM₁, (2/3, 1/6, 1/6). It comes from the S₄ quotient of the line stabiliser
3³:S₄. Two questions follow:
* what TM₁ then predicts for δ and θ₂₃;
* whether the most general S₄-invariant flavon potential actually leaves the residual symmetries TM₁ needs.

## A. Sum rules: δ against θ₂₃

Fixing the first PMNS column fixes cos δ as a function of (θ₁₃, θ₂₃) through |U_μ1|². Over the NuFIT 6.0 3σ box
quoted in Pass 11225:

| sin²θ₂₃ | 0.43 | 0.47 | 0.50 | 0.53 | 0.561 | 0.596 |
|---|---|---|---|---|---|---|
| TM₁ cos δ | −0.31 | −0.13 | **0** | +0.13 | +0.27 | +0.43 |
| TM₁ δ (or 360° − δ) | 108° | 98° | **90°** | 82° | 74° | 65° |

At these values sin²θ₁₃ = 0.02195. Across the full θ₁₃ range, cos δ spans [−0.33, +0.45], and sin²θ₁₂ = 0.318.

* **TM₁: cos δ = 0 exactly at maximal θ₂₃, for every θ₁₃** (checked exactly). The sign of cos δ is the octant of θ₂₃.
  So TM₁ predicts near-maximal CP violation, |δ − 90°| < 26° (or the mirror value about 270°), tied to the octant.
  This is a sharp test for DUNE and Hyper-K.
* **Clifford column** ((2/3)sin²(kπ/9), from the order-9 elements of the qutrit Clifford group). It fits only for
  **sin²θ₂₃ ∈ [0.586, 0.596]**, the upper edge of the 3σ range, and there δ ≈ π. The other assignment of its entries
  fits nowhere.
* **Golden A₅ columns.** One fits for sin²θ₂₃ ≥ 0.516, the other for sin²θ₂₃ ≤ 0.480.

The TM₁ atmospheric sum rule is standard (Albright–Rodejohann 2009; King's mixing sum rules). New here are the two
facts that it is the line stabiliser's column (Pass 11225) and that the Clifford column is confined to the
upper-octant edge.

## B. Vacuum alignment

The TM₁ involutions are computed, not assumed. In the basis where the charged-lepton Z₃ is the cyclic permutation, the
three involutions whose isolated eigenvector has Fourier weights (2/3, 1/6, 1/6) are diag(1,−1,−1)·P₂₃ and its two
cyclic conjugates. A neutrino flavon whose stabiliser is exactly that Z₂ points along (0,1,−1) in the triplet **3**
(or (0,1,1) in **3′**).

Method:
* the most general S₄-invariant potential is sampled as the Reynolds average of random symmetric tensors, so every
  invariant is present with a generic coefficient;
* the mass terms are tachyonic, so the symmetry breaks;
* the potential is minimised globally (24 BFGS starts);
* the vacuum's stabilisers are read off.

Scenarios: one triplet, and two triplets (φ_e, φ_ν) in all four 3/3′ assignments. The two-triplet cases come with and
without a shaping Z₂×Z₂ that separates the flavons, at degree 4 and degree 6. Vacua in which a flavon vanishes are
resampled.

| scenario | outcome |
|---|---|
| one triplet, renormalisable (degree ≤ 4), 3 or 3′ (300 potentials) | minima only (1,0,0) [Klein] and (1,1,1) [S₃ or Z₃]; **the TM₁ direction never** |
| one triplet with sextic terms, 3′ (150) | (1,1,0) with stabiliser Z₂ in **37/150**; the TM₁ direction becomes reachable |
| two triplets, every allowed cross-coupling, all 12 scenarios (**1500 vacua**) | **TM₁ in 0/1500** |

Why the two-flavon case fails:
* the cross-couplings align the two vacua;
* whenever φ_e keeps a Z₃, φ_ν sits on the same (1,1,1) axis and keeps the same group, which makes the neutrinos
  degenerate (all 1500 vacua: "G_e has no Z₃" or "G_ν too large").

**Reading.**
* The geometry hosts TM₁ (Pass 11225), and TM₁ makes a sharp δ–θ₂₃ prediction.
* An S₄ flavon potential cannot by itself put the vacuum there. With one flavon the TM₁ direction needs
  non-renormalisable terms. With the two flavons that TM₁ requires, generic cross-couplings misalign it completely
  (0/1500).
* This reproduces, for the line stabiliser's S₄, the known vacuum-alignment problem of discrete flavour symmetry. The
  standard escapes are SUSY driving fields with F-term alignment, extra dimensions, and shaping symmetries that forbid
  the dangerous cross-couplings. These are additional structure, not supplied by W33 here.
* **Masses and the remaining mixing stay OPEN.** What this pass adds is one falsifiable statement (TM₁ ⇒ cos δ ∝ the
  octant, zero at maximal θ₂₃) and one no-go (no renormalisable S₄ potential selects it).

**Scope.** These are random samples, not a proof over the whole coupling space; "0/1500" bounds the measure of the
TM₁ region and does not exclude it. Real triplet flavons only. The residual-symmetry identification follows Pass 11225.
