# Pass 11121 — bridge: the other track's Hesse / H27 / Clifford / 243 objects in the explicit string vacua

Producer: `analysis/w33_pass11121_hesse_gauge_bridge.py` (the bridge table is frozen in the certificate)
Certificate: `data/w33_pass11121_hesse_gauge_bridge.json`
Regression: `tests/test_w33_pass11121_hesse_gauge_bridge.py`

## Object by object

| other track | in the W(3,3) SO(16)×SO(16) vacua | status |
|---|---|---|
| external-A2 qutrit H27 (`2026-09-21_physical_external_a2_h27.md`) | the family Δ(54) of every model; X = fixed-point translation, Z = space-group phase (Pass 11105) | **realised** |
| H27 centre = FI Z3 (`2026-09-21_physical_fi_is_h27_center.md`) | centre = orbifold twist, a gauge element in 104/104 models, exactly the anomalous-U(1) Z3 in models 4, 25, 29 (Pass 11116) | realised as a gauge element; literally the FI Z3 in 3/104 |
| Clifford-648 = W33 point stabiliser (`2026-09-21_physical_a2_clifford648_w33_bridge.md`; 1296 with similitudes) | only H27:⟨−I⟩, index 12; quadratic phases are forbidden by the space-group rule (Pass 11105) | partial (−I only) |
| two-qutrit 3^(1+4), commutation graph W(3,3) (`2026-09-21_e8_trinification_two_qutrit_pauli243.md`) | not on generations (faithful degree 3^k forces k = 1), not on generation ⊗ colour, not on hidden SU(3) (Passes 11111, 11116) | **not realised** |
| matter 81 = 9 × H9 (`2026-09-21_e8_matter81_pauli243_restriction.md`) | needs the E6 trinification centre; these vacua carry SM × U(1)ⁿ × hidden | not realised |
| Hesse configuration / pencil (Passes 11061–11066; the 2026-09-23 Hesse notes) | the up-quark mass determinant over the Δ(54) Higgs triplet is a Hesse-pencil cubic, and the family modulus picks the member (Pass 11114) | **realised: the geometry of the quark masses** |
| AG(3,3), collinearity | θ³ couplings are collinear triples; θ × θ² × untwisted requires the same point (Pass 11105) | realised |

## The Hesse parameter as a function of the family modulus

λ(T) = (s³ + 2d³)/(3sd²), where s and d are the A2 theta functions of T.

| T | λ | curve |
|---|---|---|
| ρ (fixed by ST; the one-loop minimum) | ω² | singular member (a triangle), j = ∞ |
| i (fixed by S) | **1 + √3** exactly | j = **1728** = j(i) |
| → i∞ (the cusp) | → ∞ | the fourth singular member |

Away from the fixed points, j_curve(T) is not j(T), j(3T) or j(T/3) (checked at 1.2i, 1.5i, 2i and 0.3 + 1.1i). The
special values are the ones the fixed points force, not a hidden modular identity.

## For the other track

The objects that live in these vacua are the one-qutrit ones and the Hesse pencil:
* H27 with its centre, which here is a gauge element;
* Δ(54) = H27:⟨−I⟩;
* the Hesse pencil, as the locus of up-quark mass degeneracy over Higgs directions.

The two-qutrit 243 group and its W(3,3) commutation geometry act on no matter. W(3,3) appears only through the twist's
underlying E8 geometry (Passes 1020–1021, 68df33f).
