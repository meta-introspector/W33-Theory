# Pass 11097 — in the SO(16)×SO(16) A8 models the fractional charges are not symmetry-protected

Producer: `analysis/w33_pass11097_fractional_fermion_masses.py`
Engine: `analysis/w33_so16_selection_rules.py` (selection rules; lattice test; minimal orders by integer programs)
Frozen results: `data/w33_pass11097_symmetry_test_104.json` (all 104), `data/w33_pass11097_orders_33.json`; three sample
models re-run from their field dumps (`data/w33_pass11097_so16_a8_sample_models.json`)
Certificate: `data/w33_pass11097_fractional_fermion_masses.json`
Regression: `tests/test_w33_pass11097_fractional_fermion_masses.py`

Every one of the 104 tachyon-free three-generation models of Pass 11095 carries **136–294 fractionally charged
fermion states**. These are colour singlets with non-integer charge and colour triplets with 3Q + t ≢ 0 mod 3; 72–222
of them are hidden singlets. In Z3×Z3 the analogous states proved unremovable, protected by an unbroken discrete
3-group (Passes 11024, 11087–11090). The question here is whether scalar VEVs that preserve the Standard Model can
give them all masses.

## Rules

A mass term ψᵢψⱼ·(monomial in scalar VEVs) must satisfy:
* gauge invariance;
* the Witten Z₂ (Σk even);
* the point group (Σl ≡ 0 mod 3);
* the space group (fixed-point classes ≡ 0 mod 3 in each torus);
* H-momentum, Σ Rᵃ = q_sh + oscillator. This is exact (Σ = 0) at cubic order and mod 3 at higher order.

The candidate VEVs are the Standard-Model-neutral scalars. **Variant A** uses only hidden singlets, so the hidden
gauge group stays unbroken. **Variant B** also uses hidden-charged ones, so the hidden group breaks; this is the most
permissive case.

The existence test asks whether the charge of ψᵢψⱼ lies in the Z-span of the VEV charges plus the periodicities.
This is valid without holomorphy: φ and φ* are separate listed fields. A maximum matching over the fractional
fermions then counts how many must stay massless. Self-conjugate representations may pair with themselves.

## A. The symmetry test (104 models)

| VEVs | rules | all fractional fermions can be massive | no hidden-singlet fractional state light |
|---|---|---|---|
| A (hidden unbroken) | gauge + Z₂W + PG + SG | 51 | 61 |
| A | **+ H-momentum** | **33** | 42 |
| B (hidden broken) | gauge + Z₂W + PG + SG | 66 | 67 |
| B | + H-momentum | 45 | 48 |

## B. The order at which they become massive (33 models, variant A, all rules)

For every allowed mass term the exact minimal order is computed:
* **Order 1**: a cubic ψψφ with a single scalar, under the exact H-momentum rule.
* **Higher orders**: integer programs (HiGHS). These agree 240/240 with the exact Smith-form solver
  `w33_exact_monomial_orders` on test targets.

The bottleneck order T* is the smallest order at which **all** fractional fermions can be massive simultaneously.

| T* | 1 | 3 | 4 | 5 | 7 |
|---|---|---|---|---|---|
| models | 6 | 15 | 6 | 4 | 2 |

Every one of the 33 is fully massive by order 7. For six of them a single scalar suffices, through cubic couplings.

## Reading

Unlike in Z3×Z3, **the fractional charges of the W(3,3) SO(16)×SO(16) models are not protected by any symmetry in
a third of the models**. With Standard-Model-neutral scalars at M_s-scale VEVs ε·M_s they get masses of order
εᵀ*·M_s: heavy for any reasonable ε. In the other 71 models, under the same rules, some fractional fermions stay
massless, protected by an unbroken symmetry exactly as in Z3×Z3.

Scope. This establishes that the couplings are allowed and at what order. It does not establish that the required
VEVs form a vacuum; in a non-supersymmetric theory that is decided by the (one-loop) scalar potential; see
Pass 11099. The coupling coefficients are not computed.
