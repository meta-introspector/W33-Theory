# Pass 11029 — canonical twirl dynamics does not select the 3D clock

Producer: `analysis/w33_pass11029_exact_twirl_dynamics_firewall.py`
Certificate: `data/w33_pass11029_exact_twirl_dynamics_firewall.json`
Regression: `tests/test_w33_pass11029_exact_twirl_dynamics_firewall.py`

The exact signed (GL_2(3)) action on the minimal 24-dimensional clock carrier
admits a canonical finite-group twirl

[
mathcal T(X)=
rac1{48}sum_g U_g XU_g^T .
]

Because every (U_g) is orthogonal, the Kraus family (U_g/sqrt{48}) makes
(mathcal T) completely positive, trace preserving and unital. It is also
idempotent: the unique trace-preserving conditional expectation onto the
commutant of this finite representation.
The exact character calculation gives

[
dimoperatorname{Fix}_V(G)=2,
qquad
dimoperatorname{End}_G(V)=19.
]

Hence the vector Reynolds average has rank two, while the operator twirl fixes
a 19-dimensional algebra.

The associated finite Markov generator

[
mathcal L=mathcal T-I
]

has the exact CPTP semigroup

[
e^{tmathcal L}
=e^{-t}I+(1-e^{-t})mathcal T .
]

So a mathematically canonical dissipative dynamics exists — but it does not
produce the coarse clock field.
The naive clock augmentation has dimension three. The central signed element
(-I) takes it outside itself: adjoining its image raises the span from three
to six dimensions. Closing the same three generators under the full exact
signed group gives all 24 dimensions.

Therefore

[
oxed{	ext{no exact signed-}GL_2(3)	ext{-equivariant idempotent has the
3D clock augmentation as its image}.}
]

The reason is elementary and decisive: the image of an equivariant idempotent
must be invariant, and this particular 3D subspace is not.

This sharpens the open dynamics problem left by Passes 10975–10976. A physical
clock-selection mechanism must reduce/break the exact signed symmetry, employ
an explicit coarse quotient/readout, or add dynamics not captured by the
canonical group twirl. The theorem is finite-dimensional; it does not claim a
microscopic bath, rate, or laboratory relaxation process.
