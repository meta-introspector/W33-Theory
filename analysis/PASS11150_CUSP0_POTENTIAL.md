# Pass 11150 — near the cusp at 0 the potential is the Fricke mirror of large volume and runs into the tachyon disk

Producer: `analysis/w33_pass11150_cusp0_potential.py`
Scan: `analysis/w33_pass11150_scan_cusp0.py`
Frozen: `data/w33_pass11150_cusp0_beta_integrands.json`
Regression: `tests/test_w33_pass11150_cusp0_potential.py`

Setup: model 40521, both Wilson-line tori at iy, family torus at ρ. The small radii y = 0.20 and 0.17 are tachyon-free,
below the golden disk's lower edge at 0.2205 (Pass 11141). Five cosets suffice there (K = 5 and K = 7 agree to 10⁻⁸).

| y | Λ_β | Fricke image y′ = 1/(3y) | Λ_β(y′) | difference |
|---|---|---|---|---|
| 0.20 | 439.30 | 1.667 | 436.92 | 0.5% |
| 0.17 | 610.58 | 1.961 | 608.81 | 0.3% |

At fixed τ the partition function matches its Fricke image to 3–4·10⁻⁴ (y = 0.2) and 6–7·10⁻⁵ (y = 0.15).

* **Λ falls as y rises toward the disk** (dΛ/dy ≈ −5.7·10³). From the cusp-0 side the Wilson-line modulus rolls up into
  the tachyonic disk; toward y → 0, Λ grows, as the image of large volume.
* **Pass 11141's conclusion extends.** The tachyon-free region near cusp 0 in the six golden-type models is, to
  0.3–0.5%, the Fricke mirror of the large-volume runaway, and it has no critical point. No stabilisation of the
  Wilson-line modulus exists there either.
