# Pass 11226: an exact area law (Ryu–Takayanagi) in networks of the substrate's perfect tick

Producer: `analysis/w33_pass11226_holographic_area_law.py`
Certificate: `data/w33_pass11226_holographic_area_law.json`
Regression: `tests/test_w33_pass11226_holographic_area_law.py`

**Background.**
* The paper's scorecard lists "gravity, a dynamical space-time, the cosmological constant" as OPEN.
* The one exact bridge from quantum information to gravity is the Ryu–Takayanagi (RT) law: the entropy of a boundary
  region equals the area of a minimal surface.
* Perfect-tensor networks realise it exactly as a discrete min-cut law (Pastawski, Yoshida, Harlow and Preskill 2015).
* The geometry supplies perfect tensors:
  * the canonical two-qutrit tick p of Pass 11193, whose Choi state is AME(4,3);
  * AME(6,3) states;
  * the AME(10,3) Glynn state of Pass 11224.
* The corpus mentions holographic codes and RT only in legacy text and toy graph min-cuts. No network of the geometry's
  own tensors had been built and checked.

**Method.**
* Each tensor is a qutrit stabilizer state.
* Bonds are contracted by projecting leg pairs onto the qutrit Bell pair, which is exact F₃ linear algebra. Bulk legs
  are fixed in a product state.
* For every contiguous boundary interval A:
  * S(A) = rank(L|_A) − |A| trits, computed from the boundary Lagrangian;
  * the minimal cut separating A from the rest of the boundary is computed by max-flow, with unit-capacity bonds.

## Results

| network | tensors | boundary legs | intervals | S(A) = min-cut |
|---|---|---:|---:|---:|
| HaPPY {5,4} patch | AME(6,3), bulk legs fixed | 25 | 300 | **300** |
| flat {4,4} patch | AME(6,3), bulk legs fixed | 21 | 210 | **210** |
| {4,5} patch | **the perfect tick p alone** (AME(4,3)) | 20 | 200 | **200** |
| {4,5} patch, control | SUM (not perfect) | 20 | 200 | 96 (104 below the cut) |

* The network built only from the substrate's own perfect two-qutrit tick satisfies the discrete area law exactly on
  every interval.
* Replacing p by the non-perfect SUM gate breaks it on more than half the intervals.
* So it is the perfect (AME) property of the geometry's tick that carries the area law.

## What this does and does not do for the scorecard

* **Gravity: HOSTED, not derived.**
  * An exact entanglement = area law holds in networks of the geometry's tick.
  * Einstein's equations follow from the entanglement first law only for holographic conformal field theories. Toy
    perfect-tensor codes do not supply that dynamics, so no Einstein equation is claimed.
* **Cosmological constant: not decided.**
  * The plan was to contrast a hyperbolic {5,4} patch with a flat {4,4} patch, to see whether the area law forces
    negative curvature (Λ < 0). Both patches satisfy RT exactly, so small patches do not decide the sign. That claim is
    withdrawn.
  * Perfect-tensor holography is built on anti-de Sitter-like tilings, while the observed Λ is positive. Nothing here
    addresses that.
* **Dynamical space-time: OPEN.** These are static networks.
  * The bridge to this session's arrow results is that the same perfect tick has A = 0 (Pass 11223): no intrinsic arrow,
    yet an exact area law.
  * The tick that carries the arrow, the three-qutrit V K V of Pass 11223, is a natural candidate for a dynamical
    (circuit) version. That is left open.
