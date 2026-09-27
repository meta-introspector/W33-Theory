# Pass 11039 — the 4×9 Bell-relative atlas has exact rank 21

Producer: `analysis/w33_pass11039_bell_relative_tomography_rank.py`
Certificate: `data/w33_pass11039_bell_relative_tomography_rank.json`

Fix the Bell line. Its complement is four canonical nine-point sectors. The 27 Lagrangian lines transverse to the Bell line each contain exactly one point from each sector.

The 36×27 point-line incidence matrix has:
- row weight 3;
- column weight 4;
- rank 21 over Q and over F2, F3, F5, F7, F11;
- right kernel dimension 6;
- left kernel dimension 15;
- Gram spectrum `12^1, 6^12, 3^8, 0^6`.

Adding complete nine-point response sheets grows the rank as
`9 -> 15 -> 19 -> 21`,
so the successive increments are `9,6,4,2`.

Thus the proposed 4×9 measurement layout is an overcomplete response frame, not 36 independent observables. It contains 15 exact response redundancies and leaves a six-dimensional context kernel unless extra gauge/normalization information is supplied.

This is a clean target for a future CP-completed W33 tomography experiment.
