# Pass 11141 — every symmetric point of the Wilson-line modulus is tachyonic or a runaway

Producer: `analysis/w33_pass11141_no_symmetric_stabilization.py` (reads the certificates of 11136, 11139 and 11146)
Regression: `tests/test_w33_pass11141_no_symmetric_stabilization.py`

A potential invariant under the duality group Γ0(3) (Pass 11139) is automatically stationary at the group's elliptic
point and extremal toward its cusps. Γ0(3) has genus 0, one elliptic point of order 3 at (3 + i√3)/6, and two cusps (∞ and
0).

| symmetric point | what is there |
|---|---|
| elliptic point (3 + i√3)/6 | **the exact centre of the charged tachyon disk**: p_R² = ⅔, Δ = −⅙, in every model |
| cusp ∞ | large-volume runaway: Λ grows with the volume, so the moduli roll inward |
| cusp 0 (diagonal) | tachyonic at Im T = 0.2 in 15/21. In the other 6 (46043, 24165, 40521, 5904, 5285, 2080) the tachyon ends at the golden disk's lower edge (0.2205): the near-cusp region is tachyon-free, the approximate Fricke image of large volume, and by that approximate symmetry a runaway into the disk from below |
| Fricke point i/√3 | only an approximate symmetry point, and tachyonic (centre of the golden neutral disk) |

The tachyon-free part of the fundamental domain is a neighbourhood of the cusp at ∞ above the disks, plus, in six
models, a neighbourhood of the cusp at 0 below the golden disk. (The single-torus zero-radius limit is the U(16) string,
11140, which is a different direction.) There Λ_β is
monotone in Im T (11122) with a near-flat B direction (11129), so it has no critical point. **Stabilising the Wilson-line
modulus at a self-dual point is excluded in these vacua.** This is sharper than, and consistent with, Fraiman et al.'s
finding for circles (arXiv:2307.13745) that points of maximal enhancement are unstable or lie on the edge of tachyonic
regions.
