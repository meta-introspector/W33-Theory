# Pass 11107 — the winding tachyons: explicit spectrum, critical radii, fractional charges

Producer: `analysis/w33_pass11107_winding_tachyons.py` (the enumerator; `--check` reruns the controls)
Data: `data/w33_pass11107_critical_radii_12.json`, `data/w33_pass11107_tachyon_quantum_numbers_12.json`
Certificate: `data/w33_pass11107_winding_tachyons.json`
Regression: `tests/test_w33_pass11107_winding_tachyons.py`

Pass 11106 found a level-matched tachyon in the one-loop integrand when the Wilson-line tori shrink. This pass builds it
state by state.

## The states

The tachyon lives in the Witten-twisted sector. It has the right-moving ground state r = −v0 (q_R = 0) and no
oscillators, so

    Δ = −½ + p_R²/2,     ℓ² + p_L² − p_R² = 1,     ℓ = π + V0 + A n  (π ∈ E8×E8).

Two projections apply:
* the β projection keeps states with e^{2πi(π·V0 + V0²)} = −1, the sign fixed by Z[1,1] = −H[1,1] in the one-loop engine;
* the Z3 projection keeps one combination per orbit of three rotated winding states.

The torus momenta are the Narain vectors of `w33_so16_one_loop.py`, with the family torus at T\* = ρ.

## Checks

| check | result |
|---|---|
| lowest level at Im T_WL = 1 | **Δ = −0.211325 = −(3−√3)/6**; the one-loop integrand's growth gives −0.205 |
| large radius (T_WL = T\* = 5i) | no tachyon: the 10D SO(16)² string is tachyon-free |
| Im T_WL = 2 | no tachyon in any of the 12 survivors |

## Critical radii (both Wilson-line tori at i·y)

| y_c | models |
|---|---|
| √3 ≈ 1.73 | 2, 14, 53, 102 |
| 1.51 | 10, 13, 15, 35, 57, 69, 78 |
| tachyon-free on the whole imaginary axis, tachyonic at the SU(3) point ρ | 77 |

**Correction to Pass 11106.**
* Its onset "Im T_WL ≈ 1.4" came from a growth-ratio test. That test cannot see a tachyon with |Δ| ≲ 0.05 under the
  evolving Kaluza–Klein tower. The enumerated onset is 1.51–1.73.
* Its Λ value at T_WL = 1.5i for model 2 lies inside the tachyonic region, where the integral diverges, so that value
  is void.
* The monotone fall of Λ holds on the tachyon-free segment 1.75 ≤ Im T_WL ≤ 3.

## Charges

Every winding tachyon is a colour singlet. Their electric charges:
* **In 10 of 12 models every tachyonic state has fractional electric charge:** ±1/3 or ±2/3. For example, in model 2
  there are SU(2) doublets with Y = ±1/6 and singlets with Y = ±1/3 and ±2/3. Any condensate therefore breaks U(1)_EM.
* In models 53 and 57, some degenerate tachyons are fully Standard-Model-neutral (Y = 0, SU(2) singlets). A condensation
  direction that preserves the Standard Model exists there, but nothing selects it over the charged ones.

## Reading

The small-radius instability of Pass 11106 is a fractionally charged winding tachyon. This is the same disease as the
twin route, whose tachyons were all half-charged or coloured (Pass 11094). The one-loop potential pushes the Wilson-line
moduli toward it (Pass 11106, 11110). Only models 53 and 57 leave room for a neutral condensate, and computing its
endpoint needs the tachyon potential, which is not done here.
