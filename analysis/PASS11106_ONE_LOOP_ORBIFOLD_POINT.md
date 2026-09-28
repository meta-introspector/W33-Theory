# Pass 11106 — the one-loop potential at the orbifold point: the family modulus is stabilised, the Wilson-line moduli roll into a winding tachyon

Engine: `analysis/w33_so16_one_loop.py`
Producer: `analysis/w33_pass11106_one_loop_orbifold_point.py`
Data:
* `data/w33_pass11106_survivor_shifts.json` (the 12 survivors' V0, V, W)
* `data/w33_pass11106_beta_integrands_model2.json`
* `data/w33_pass11106_twisted_integral_model2.json`
* `data/w33_pass11106_tachyon_onset_12*.json`
* `data/w33_pass11106_massless_bose_fermi_corrected.json`

Certificate: `data/w33_pass11106_one_loop_orbifold_point.json`
Regression: `tests/test_w33_pass11106_one_loop_orbifold_point.py`

## Background

Pass 11099 found the vacuum energy positive and growing with the volume at large volume: Λ₄ → (V/3)Λ₁₀. It left the
orbifold point, meaning small volume with all twisted sectors included, uncomputed. This pass computes the full one-loop
integral over the fundamental domain.

## The engine

G = Z2W × Z3 is cyclic of order 6. All 36 sectors are built. Wilson lines satisfy W1 = W2, W3 = W4 and W5 = W6; exactly
one torus has none (Pass 11105).

* **The 4 Witten-type sectors** (g, h ∈ {1, β}; β is trivial on T⁶) form (1/3)·Z[SO(16)² on T⁶]. They are Hamiltonian
  sums over the Narain lattice Γ₆,₂₂ with Wilson lines:
  * the Wilson-line tori are handled by decomposing into the cosets π·W mod 1, which gives E8 theta functions with
    characteristics;
  * the Wilson-line-free torus factorises as Γ₂,₂(T\*, ρ);
  * the vacuum phases are fixed by modular covariance, with Z[1,1] = −H[1,1].
* **The other 32 sectors** form one SL(2,Z) orbit per fixed point:
  * the seed block (the untwisted sector with insertion βθ_f, local shift V0 + V + Σ c_t W_t) is Γ₁(6)-invariant to
    10⁻¹³;
  * 24 coset images × 27 fixed points give every twisted sector;
  * the 8 pure-Z3 pairs vanish identically.

### Checks

* The twisted sum and the Witten sum are each modular invariant under S and T, to 10⁻¹².
* The integrator reproduces the 10D integral of Pass 11099: I = −725.97.
* Large volume gives τ₂⟨Z_β⟩ → (1/3)(−2112)V/τ₂³: −712 vs −704 at T = 2.5i.
* **The massless coefficient matches the orbifolder spectrum sector by sector (model 2):**

| sector | partition function | spectrum |
|---|---|---|
| untwisted | 122 | 104 listed + 14 (the 7 U(1) gauge bosons, not listed in the dump) + 4 (gravity) |
| Witten-twisted | −68 | −68 |
| θ, θ² without β | 174 | 87 + 87 |
| βθ, βθ² | −472 (still converging on a light level at τ₂ = 3) | −480 |

The total is −252. Correction to Pass 11099: its n_B − n_F omitted the U(1) gauge bosons. The corrected range over the
104 models is **[−480, −48]** (was [−500, −64]), still never zero, so that conclusion stands.

## Results for model 2 (Λ₄ = −½ M⁴ I)

The twisted sectors give I_tw = −74.26, independent of the moduli.

| moduli | Λ₄ / M⁴ |
|---|---|
| T_WL = 2i; T\* = ρ | **657.6** |
| T\* = 0.25 + i | 676.5 |
| T\* = i | 687.4 |
| T\* = 0.5 + 1.2i | 703.1 |
| T\* = 1.2i | 710.7 |
| T\* = 1.5i | 800.0 |
| T\* = 2i | 1006.1 |
| T\* = ρ; T_WL = 3i | 1464.9 |
| T_WL = 2.5i | 1022.9 |
| T_WL = 2i | 657.6 |
| T_WL = 1.75i | 501.7 |
| T_WL = 1.5i | 355.9 |

* **The family modulus is stabilised at its self-dual point.** V(T\*) is exactly SL(2,Z)-invariant, because the
  twisted sectors do not depend on T\*. Λ is lowest at T\* = ρ, the SU(3) point, rises along the arc to i (a saddle),
  and rises upward.
* **The Wilson-line moduli roll inward.** Λ decreases monotonically as their radius shrinks, and it stays positive.
* **At small radius a tachyon appears.** A level-matched winding/momentum state of the Witten-twisted sector becomes
  tachyonic. It has nonzero torus momentum, so the zero-momentum orbifolder spectrum cannot see it.
  * Model 2 is tachyonic for Im T_WL ≲ 1.4: m²/4 ≈ −0.13 at 1.2, −0.21 at 1.0, and −0.23 at 0.87.
  * All 12 survivors are tachyon-free at Im T_WL = 2.
  * Every survivor is tachyonic at some smaller radius: 4 already at 1.2 (models 2, 14, 53, 102), 11 at 1.0, and
    model 77 only at the SU(3) point.
  * The fully symmetric point, with all T at the SU(3) point, is tachyonic.

## Reading

The positive one-loop potential cannot be balanced inside the region where the models are tachyon-free:
* the Wilson-line moduli run into the tachyonic region;
* the dilaton tadpole, since Λ > 0, remains;
* only the family modulus is stabilised, at T\* = ρ.

Toroidal compactifications of the O(16)×O(16) string are known to develop tachyons at small radius (Ginsparg and Vafa,
Nucl. Phys. B289 (1987) 414; Itoyama and Taylor 1987). What is new here is that the W(3,3) models reach that region
dynamically. The end point of tachyon condensation is not computed.

Scope: one model (2) for the potential, 12 for the tachyon onset. The WL-torus moduli are taken equal; the diagonal
and T\* directions are sampled, not minimised over the full moduli space.
