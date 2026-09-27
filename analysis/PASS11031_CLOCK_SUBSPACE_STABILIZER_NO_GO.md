# Pass 11031 — the coarse 3D clock has trivial stabilizer in the exact signed action

Producer: `analysis/w33_pass11031_clock_subspace_stabilizer_no_go.py`
Certificate: `data/w33_pass11031_clock_subspace_stabilizer_no_go.json`
Regression: `tests/test_w33_pass11031_clock_subspace_stabilizer_no_go.py`

Pass 11029 showed that the full signed GL₂(3) orbit of the three coarse clock-augmentation generators spans all 24 noncentral coordinates. This pass asks the sharper subgroup question.

For every one of the 48 exact signed group elements, we test whether it preserves the three-dimensional augmentation span setwise.

Result:

- stabilizer order = 1;
- the identity is the only preserving element;
- every one of the other 47 elements raises the joint span of the clock subspace and its image from rank 3 to rank 6;
- central −I already produces rank 6;
- full orbit closure remains rank 24.

Therefore there is no nontrivial subgroup of the current exact signed GL₂(3) action that can preserve this particular coarse 3D clock as an invariant linear image. A linear dynamical selector for this embedding must completely break the signed symmetry, or use a quotient/readout map instead of invariant-subspace selection.

Boundary: this is exact for the frozen signed H27 embedding. It does not forbid a different 3D field, nonlinear order parameter, or a jointly transformed gauge/readout convention.
