# Pass 11115 — the duality turnaround of the Wilson-line potential never lands in a tachyon-free window

Producer: `analysis/w33_pass11115_axis_duality_and_window.py`
Data: `data/w33_pass11115_axis_tachyons.json`, `data/w33_pass11115_axis_integrands_model77.json`
Certificate: `data/w33_pass11115_axis_duality_and_window.json`
Regression: `tests/test_w33_pass11115_axis_duality_and_window.py`

## The idea

At large volume Λ grows like the volume in any T-dual description. So along the zero-B-field axis of the Wilson-line
tori, V(y) must turn around somewhere at small radius. If a model were tachyon-free there, its Wilson-line moduli would
be stabilised.

## Tests

The explicit enumerator of Pass 11107 was run with windings and momenta up to 7 and cached torus tables. The one-loop
engine used K = 6 windings below Im T = 0.9, converged to 10⁻⁹ at Im T = 0.5.

| model | tachyon levels on the axis, y = 0.3 … 1 |
|---|---|
| 10 | approximate **Fricke duality** y → 1/(3y): paired levels agree within 0.03; deepest at the self-dual radius 1/√3 (Δ = −1/6). The duality extremum lies inside the tachyonic region |
| 2, 53 | tachyonic throughout, increasingly so |
| 57 | tachyonic throughout |
| **77** (a full survivor) | tachyon-free at y = 0.87, 0.9, 1.0; tachyonic at 0.85 and 0.8. The first tachyon has charge ±1/3. Onset ≈ √3/2 |

Model 77's potential on its window:

| y | 2.5 | 2.0 | 1.5 | 1.25 | 1.1 | 1.0 | 0.9 |
|---|---|---|---|---|---|---|---|
| Λ_β/M⁴ | 993.7 | 633.2 | 350.8 | 238.6 | 180.9 | 146.7 | 116.0 |

It is monotone all the way down to the onset.

## Reading

The turnaround that duality requires does not happen inside a tachyon-free window in any model tested:
* model 10's duality extremum sits at its tachyonic self-dual radius;
* model 77, the survivor that is safe furthest down the axis, rolls into a fractionally charged tachyon.

The Wilson-line moduli are not stabilised in the tachyon-free region. The only SM-safe exit is the neutral corridor of
models 53 and 57 (Pass 11112).
