# Passes11443–11447: preserved gauge, Lorentzian matter and a UV obstruction

Reservation `bad2527a6`. [Producer](w33_pass11443_11447_preserved_gauge_lorentz_recovery.py),
[certificate](../data/w33_pass11443_11447_preserved_gauge_lorentz_recovery.json),
[independent tests](../tests/test_w33_pass11443_11447_preserved_gauge_lorentz_recovery.py).

All five requested directions were investigated with explicit maps or negative
controls. The ultimate global chiral measure, dynamical Lorentzian gravity,
full joint vacuum and noisy native recovery remain unfinished. The strongest
positive result is a lower-energy native candidate preserving eight continuous
generators; the strongest obstruction is the elementary heavy spectrum's
loss of colour asymptotic freedom. Neither establishes a TOE.

## Ownership and intake

Built on published11438–11442 (`6b7fa867f`, receipt `ca849a0ee`) after
integrating parallel inventory `b71eff8fc` through `cea396cfc`. Searches used
RESULTS_INDEX, docs/index.html, the paper and prior scripts/certificates, with
result strings for the beta coefficient, heavy multiplicity, relaxed quartic,
Lorentzian Schur map and spectroscopy. The native tensors, coercive action,
Dirac masses, overlap theory and Steane code are prior-owned. General scalar
finite elements, chiral anomaly cancellation, one-loop running and erasure
correction are established methods; the claims below concern their explicit
application and limits in the supplied repo models.

##11443: native CP with an unbroken subgroup

Take the eight-generator stabilizer of the prior single-field native vacuum.
Its common complex fixed space in the81-field representation has dimension27.
Restrict **both** fields to that space and minimize the same supplied11438
coercive action, without adding another term. The candidate has

- energy `-2.949259565341485`, below the old `-2.935743239533`;
- full324-gradient norm `1.29e-6`;
- common orbit rank78, hence eight continuous stabilizer generators;
- cross-pin CP flux `.00070093647833`;
- fixed-space residual `8.16e-16`, Lie bracket closure residual `1.38e-15`.

Compact invariance makes the gradient lie in the common fixed space; the
restricted criticality is therefore meaningful in the full action. A full
finite-difference324 Hessian gives246 normal directions:242 are positive above
`1e-3`, with smallest positive eigenvalue about `.4592`; four are unresolved
at numerical resolution. The first four eigenvalues are approximately
`[-3.62e-7,-2.65e-7,9.08e-12,1.78e-7]`. Smaller-step directional errors are of
order `1e-7`. There is no resolved negative normal mode, but **this does not
prove an isolated stable or global vacuum**. The eight-generator embedding
must be identified before it is called physical colour or an SM gauge vacuum.

The old four soft modes supply a warning. Eliminating the quadratic massive
response using the inverse positive normal Hessian reduces their apparent
quartic coefficients from roughly `5.8–16.5` to `3e-5–.0037` in the tested
amplitudes. Unrelaxed line stiffness is consequently unsuitable evidence for
coupled stability. These finite-amplitude, finite-precision scans do not prove
exact vanishing, positivity of the mixed quartic, or a flat moduli space.

##11444: what a gauge-orbit section can and cannot establish

For a Weyl frame transforming as `E_q(A^G)=G_q E_q(A)`, the ordinary chiral
output frame has two spin components per lattice site. The determinant cocycle
is `exp(2 i q sum_x alpha_x)`. The stored hypercharge multiplet cancels this
factor through its zero linear anomaly. It also has zero cubic anomaly.

The negative control `[1,1,-2]` cancels the orbit factor but has cubic anomaly
`-6`. Thus constructing this nonlocal orbit section does **not** construct the
local current, an anomaly-safe global section on the quotient, or sector
gluing. A nontrivial sector is instantiated separately with actual periodic
unit-flux links at L34. Charged flux integrals equal q; maximum q6 plaquette
distance from unity is `.03261024 < 1/30`. No full34^4 overlap diagonalization
or inter-sector transition current is claimed. Luescher's theorem remains an
external construction whose required current has not been implemented here.

##11445: an off-shell Lorentzian matter map

Use the actual curved-boundary simplex coordinates, declaring the metric
`eta=diag(-1,1,1,1)`. Its edge Gram has signature `(3,1)`. For affine scalar
shape gradients B and absolute simplex volume V, define

`K = V B^T eta B - m^2 V (ones+I)/30`.

This is the explicit quadratic scalar action `f^T K f/2`. Subdivide around an
interior point, assemble the six-variable action, and eliminate its scalar by
a stationary Schur complement. Two interior positions and a Lorentz boost
provide independent realizations. At m0 the boundary action equals the coarse
one to about `1.5e-15` for arbitrary boundary fields; volumes add exactly.
At m.2 the Schur action differs from the naive coarse massive action by about
`3e-6`, while boost covariance remains at machine precision. This gives an
actual off-shell matter coarse map, including a massive refinement correction.
Signature and boundary geometry are supplied. Lorentzian gravitational boost
angles, light-cone branches, varying geometry and Newton coupling remain open.

##11446: the heavy inventory changes the ultraviolet question

Re-read the11423 certificate's representations, rather than counting colour
states as independent flavours. Q contributes480 Dirac colour triplets;
U,D contribute480; F_u,F_d contribute6. Add the six light quarks:

`N_f = 972`, `b0 = 11 - 2 N_f/3 = -637`.

Normalization is `beta(g)=-b0*g^3/(16*pi^2)`, with `T(fund)=1/2`. This conclusion
assumes every stored triplet is elementary and no additional coloured scalar
or gauge sector changes the coefficient. The producer uses all two486
singular spectra as thresholds and integrates inverse g² between thresholds.
The extension is not asymptotically free in that interpretation.

Above all thresholds, the one-loop extrapolation is
`log(Lambda/mu)=8*pi²/(637 g(mu)²)`. Supplied example couplings `.1,.2,.5,1`
give ratios about `241617,22.17,1.642,1.132`. These are not measured couplings;
perturbation theory fails near the extrapolated pole. A successful joint
pin/link/Majorana vacuum cannot by itself repair this ultraviolet problem.
The full joint vacuum remains open; this pass prioritizes the concrete gauge
running obstruction rather than supplying another arbitrary vacuum condition.

##11447: operational erasure recovery and optimal placement

For one or two **known** erased sites, reset and independently Pauli-twirl the
sites. This produces the erasure marginal channel. Measure all six Steane CSS
stabilizers; every one of the `4^e` Pauli branches on the known mask has a
unique syndrome, so its inverse restores the encoded state deterministically.
Tests use the actual128-state code frame, all21 two-site masks and all branches,
as well as a direct reset/twirl density-matrix identity. Prior11442 owns the
complete abstract channel; this is a process using syndrome primitives.
Additional unknown faults require a separate audit; `2t+e<3` does not cover
one erasure plus one unknown fault.

The actual native clock has a two-dimensional dark space equal to the logical
frame. Uniform coherent clock interrogation has zero-phase acceptance
`|sin(M theta/2)/(M sin(theta/2))|²` on a leaked eigenmode. At M1024 the maximum
false acceptance is `.000130785`; the logical residual is `1.07e-16`.
Controlled clock, coherent interrogation/readout and reset remain supplied.
This is a detector design for an isolated cell, not a compiled noisy threshold.

The supplied80-copy A80 architecture admits a global placement optimum for the
previous flagged protocol's route count. Put the flag next to the syndrome
hub and the three highest-weight data columns at the other three neighbours.
Remaining data lie at distance2. The degree4 lower bound is attained:
`6*(7+7*5)+36=288` routed CNOTs; conditional extraction uses372, versus the old
1836 and2220. Conditional CNOT-only ticks fall to `40040220`. Preparation,
measurement, clock interrogation, noisy SWAP propagation and recovery are
excluded; no threshold is inferred from the resource improvement.

## External checks used

- [Luescher, Abelian chiral gauge theories](https://arxiv.org/abs/hep-lat/9811032): local/global measure conditions and all-sector theorem.
- [PDG QCD review, section9.1.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf): fundamental normalization and threshold-dependent one-loop running.
- [Jercher and Steinhaus, Lorentzian Regge cosmology with scalar matter](https://arxiv.org/abs/2312.11639): prior Lorentzian matter/gravity setting; no claim that our scalar map constructs their gravitational dynamics.
- [Aliferis and Terhal, fault-tolerant computation with leakage](https://arxiv.org/abs/quant-ph/0511065): leakage reduction requires operations beyond an abstract erasure witness.
