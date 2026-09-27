# Pass 11029 — canonical twirl dynamics does not select the 3D clock

Producer: `analysis/w33_pass11029_exact_twirl_dynamics_firewall.py`
Certificate: `data/w33_pass11029_exact_twirl_dynamics_firewall.json`
Regression: `tests/test_w33_pass11029_exact_twirl_dynamics_firewall.py`

The exact signed GL₂(3) action on the minimal 24-dimensional clock carrier
admits the canonical finite-group twirl

T(X) = (1/48) Σ_g U_g X U_gᵀ.

Every U_g is orthogonal, so K_g = U_g/√48 is a Kraus family. Therefore T is
completely positive, trace preserving, unital, and idempotent.

The exact fixed dimensions are

dim Fix_V(G) = 2,

dim End_G(V) = 19.

Thus the vector Reynolds projector has rank two, while the operator twirl is
a conditional expectation onto a 19-dimensional commutant algebra.
The Markov generator

L = T − I

has the exact CPTP semigroup

exp(tL) = exp(−t) I + (1 − exp(−t)) T.

So there is a canonical finite-dimensional dissipative dynamics associated
with the exact symmetry.

But it does not select the coarse three-dimensional clock. The naive clock
augmentation has dimension 3. Acting with the central signed element −I and
adjoining its image raises the span to dimension 6. Closing the same three
generators under the full exact signed group spans all 24 dimensions.

Hence no exact signed-GL₂(3)-equivariant idempotent can have precisely the
naive 3D clock augmentation as its image: the image of an equivariant
idempotent must itself be invariant.

A physical mechanism selecting the coarse 3D clock must therefore reduce or
break the exact signed symmetry, use an explicit quotient/readout map, or add
dynamics beyond the canonical group twirl.

This is a finite representation/channel theorem. It is not a microscopic
open-system derivation, a measured relaxation rate, or evidence that nature
implements this twirl.
