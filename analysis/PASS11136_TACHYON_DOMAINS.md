# Pass 11136 — the tachyon domains are exact hyperbolic disks at the special points of Γ0(3)

Producer: `analysis/w33_pass11136_tachyon_domains.py`
Scan: `analysis/w33_pass11136_scan_domains.py` (22 points per model, including B = 1.0 and 1.5)
Frozen: `data/w33_pass11136_domain_samples.json`
Regression: `tests/test_w33_pass11136_tachyon_domains.py`

Setup: both Wilson-line tori at T = B + iy, family torus at ρ. Δ = −½ + p_R²/2, with b = B − nearest integer and
c = B − ⌊B⌋ − ½. Verified at every sample to 10⁻⁹; period 1 in B.

| instability | p_R² | domain (p_R² < 1) | hyperbolic description |
|---|---|---|---|
| neutral, √3 models | (b² + y²)/(√3 y) | b² + (y − √3/2)² < ¾ | **horodisk at the cusp 0**: Im(−1/T) > 1/√3 |
| neutral, golden models | (b² + y²)/(√3 y) + 1/(3√3 y) | b² + (y − √3/2)² < 5/12 | disk centred at the **Fricke point i/√3**, radius **2 ln φ** |
| charged, both | (c² + y² + 1/12)/(√3 y) + ⅓ | c² + (y − 1/√3)² < ¼ | disk centred at the **elliptic point (3 + i√3)/6**, radius **ln(2 + √3)** |

* The golden ratio of Passes 11123 and 11128 is e^{R/2} for the neutral disk around the Fricke point. Pass 11115 had seen
  an approximate Fricke duality y → 1/(3y) in model 10; here the Fricke point is the exact centre of the neutral domain.
* The √3 domain is the S-image of the large-volume cusp. This is the T-dual reading of Pass 11128, now exact in B.
* In the golden models the charged-disk level also contains SM-neutral states.
* **The neutral instability comes first on the whole B circle.** The smallest gap between the neutral and charged tops is
  at B = ½: y_N = √3/2 + √(1/6) = 1.2743 (golden; 1.5 for √3) against y_C = 1/√3 + ½ = 1.0774. A descending Wilson-line
  modulus therefore meets the SM-neutral instability first at every B, exactly. This strengthens Pass 11125, which had
  a grid of 0.1.

Scope:
* Two models sampled; the charged form is identical in both.
* Single-winding states only; higher windings matter only near the real axis.
* Equal T on both Wilson-line tori.
* Z3 Wilson lines generically reduce the torus duality group to a level-3 congruence subgroup. That is consistent with the
  Γ0(3) special points found here, but no duality-group derivation is attempted.
