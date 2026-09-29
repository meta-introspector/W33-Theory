# Pass 11118 — the Wilson-line flow of model 57: the neutral boundary wins from equal radii, but only by a narrow margin

Producer: `analysis/w33_pass11118_wilson_line_flow.py`
Data: `data/w33_pass11118_grid_integrands_model57.json` (36 one-loop integrals, K = 4 windings)
Certificate: `data/w33_pass11118_wilson_line_flow.json`
Regression: `tests/test_w33_pass11118_wilson_line_flow.py`

Model 57's neutral winding tachyons live on Wilson-line torus 1 and its charged ones on torus 2. Both appear at the same
critical radius, Im T = 1.512 (Passes 11107, 11112).

## Method

* Λ(y1, y2) was computed on the 6 × 6 grid y ∈ {1.55, 1.7, 1.85, 2.0, 2.3, 2.6}, with B = 0 and the family torus at ρ,
  then interpolated with a bicubic spline.
* The flow uses the moduli-space metric: the kinetic term is dy²/(4y²) per torus, giving dy_i/dt = −4y_i² ∂Λ/∂y_i.
* It is integrated until it leaves the grid at y = 1.55, just above the critical radius.

## Results

* **Asymmetry everywhere.** All 15 antisymmetric differences Λ(y_small, y_large) − Λ(y_large, y_small) are negative.
  The potential is always lower with the neutral torus the smaller one.
* **No interior minimum**, as before.
* **Equal starting radii** (2.6, 2.3 and 2.0 on both tori): the flow reaches the **neutral** boundary first in every
  case.
* **The separatrix.** Starting with torus 2 smaller by δ, the neutral boundary still wins for δ ≤ 0.010 (≤ 0.015 at the
  top of the grid). The charged one wins from δ = 0.015–0.020 on:

| start (2.6, 2.6−δ) | δ = 0 | 0.005 | 0.01 | 0.015 | 0.02 | 0.04 |
|---|---|---|---|---|---|---|
| first boundary | neutral | neutral | neutral | neutral | charged | charged |

## Reading

The SM-preserving exit of Pass 11112 is real but narrow. It is reached from symmetric initial conditions and from any
start within about 1% of them on the charged side. A larger initial asymmetry toward the charged torus sends model 57
into the fractionally charged tachyon. The vacuum's fate therefore depends on the initial Wilson-line moduli, which is
a cosmological question. The potential tilts the outcome toward the neutral side without forcing it.

Scope: an overdamped (gradient) flow, i.e. the direction of steepest descent, not a cosmological trajectory with
kinetic energy. The family modulus is fixed at ρ, and the twisted sectors add a moduli-independent constant.
