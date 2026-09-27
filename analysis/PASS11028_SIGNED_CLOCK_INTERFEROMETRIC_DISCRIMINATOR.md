# Pass 11028 — a phase-scanned interferometer sees the hidden signed clock

Producer: `analysis/w33_pass11028_signed_clock_interferometric_discriminator.py`
Certificate: `data/w33_pass11028_signed_clock_interferometric_discriminator.json`
Regression: `tests/test_w33_pass11028_signed_clock_interferometric_discriminator.py`

The central support element (-Iin GL_2(3)) is invisible on the four
projective clock labels. A purely coarse four-clock model therefore treats it
as the identity.

The exact 24-mode signed carrier does not. Let (v_d) be the uniform
six-mode state on clock direction (d), and let (z) be the exact signed
central operation. The executable calculation gives

[

rac{langle v_d,zv_d
angle}{langle v_d,v_d
angle}
=
left(
rac13,-
rac13,
rac13,-
rac13
ight).
]
The signs are gauge/direction dependent, but the phase-scanned visibility is
not:

[
oxed{

rac{|langle v_d,zv_d
angle|}{langle v_d,v_d
angle}=
rac13
quad	ext{for all four clock directions}.}
]

Thus an equal-arm identity-versus-(z) interferometer predicts

[
P_{max}=
rac23,qquad P_{min}=
rac13
]

for the exact signed carrier, versus

[
P_{max}=1,qquad P_{min}=0
]

for the projective-identity model.
Operationally (z) is a signed permutation of the 24 noncentral modes:
twelve disjoint swaps, with six swap-pairs carrying the negative sign.
A coherent photonic network can represent that ideal operation with mode
routing plus calibrated (pi)-phase shifts.

The statistic is deliberately permutation-blind: ordinary (S_4) relabeling
can interchange the four fibres and flip which ones carry the positive versus
negative overlap, but the scanned visibility remains (1/3).

There is an independent numerical cross-check. Pass 10952 found that the
qutrit Clifford representation of the same abstract central element has

[
|operatorname{tr}U_{-I}|/3=1/3.
]

That equality is **not** an identification of representations; Pass 10952
already proves the qutrit Clifford is different from the signed 24-mode
carrier. It is simply a second operational appearance of the same contrast
number. Real hardware still requires loss, phase-noise and detector-visibility
calibration.
