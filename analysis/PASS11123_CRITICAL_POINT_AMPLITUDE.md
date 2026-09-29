# Pass 11123 — the critical radius: √3 and φ²/√3, no enhanced symmetry, and no tree-level contact quartic

Producer: `analysis/w33_pass11123_critical_point_amplitude.py`
Data: `data/w33_pass11123_critical_point_pairs.json`, `data/w33_pass11123_radion_slopes.json`
Certificate: `data/w33_pass11123_critical_point_amplitude.json`
Regression: `tests/test_w33_pass11123_critical_point_amplitude.py`

This pass studies the six survivors whose instabilities are all Standard-Model-neutral (Pass 11119), with both
Wilson-line tori at i·y, B = 0 and the family torus at ρ.

## A. Exact critical radii

Solving p_R²(y) = 1 for the neutral state's own Narain vector gives two values:

| models | y_c |
|---|---|
| 36621, 17224 | **√3** = 1.732050808 |
| 46043, 24165, 40521, 5904 | **φ²/√3** = 1.511522628, where φ is the golden ratio (3y_c² = φ⁴) |

At y_c every neutral state has P_L² = 2 and p_R² = 1 exactly: a left current times a right weight-½ operator.

## B. Not an enhanced-symmetry point

No pair of the massless states sums to an untwisted vector with p_R = 0 and P_L² = 2. The pair sums have
(P_L², p_R²) = (0, 0) for conjugate pairs, (3, 3) or (5, 1) at √3, and (φ², φ²) or (5.382, 1.382) at φ²/√3. So:
* no new massless gauge bosons appear;
* there is no D-term-like quartic;
* there is no moduli trapping of the Kofman–Linde–Liu–Maloney–McAllister–Silverstein type.

## C. The tree-level four-point amplitude T T → T T at y_c

Conventions: α′ = 2, two vertices in the (−1) picture and two in the (0) picture. The right-moving fermion contraction
gives k₂·k₃ + p_R₂·p_R₃. The integral is the complex beta function of Kawai–Lewellen–Tye, checked numerically to 10⁻¹⁰.

    A = Cπ Γ(3 − s/2) Γ(−t/2) Γ(−u/2) / [Γ(s/2 − 1) Γ(2 + t/2) Γ(2 + u/2)]

It is symmetric under t ↔ u, as Bose symmetry requires. At low energy,

    A = 4Cπ (1/t + 1/u) + 3Cπ s²/(tu) + O(E²):

pure massless exchange, with **no analytic contact term**. At tree level nothing stops the condensate at a small VEV.

## D. The radion coupling

dΔ/dy at y_c is **1/(2√3) = 0.2887** in the √3 models and 0.2466 in the φ models. Both are positive, so a condensate |T|²
lowers the energy as the torus shrinks. T and the radius run together, and the endpoint lies outside the perturbative
neighbourhood of y_c. This is the localized/winding-tachyon regime of Adams, Liu, McGreevy, Saltman and Silverstein, "Things fall apart" (hep-th/0502021), and of Horowitz (hep-th/0506166), where the condensate pinches off the cycle and changes the topology.

## Reading

The neutral exit is real, but it does not end in a nearby controlled vacuum. Condensation drives the Wilson-line torus
to small size, beyond tree-level control. The next step would need the combined (T, radion) potential beyond tree
level, or a known endpoint of winding-tachyon condensation.
