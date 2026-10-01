# A D-flat kinetic bridge and an explicit mass operator: the simplest loop is a saddle

Producer: `w33_20261001_dflat_cartan_mass_loop.py`.
Certificate: `data/w33_20261001_dflat_cartan_mass_loop.json`.

This extends the exact Cartan/global polynomial selector with actual canonical
field kinetics and a named conditional Yukawa operator. It separates the
numerical gauge normalization from exact design identities and interval signs.

## Full-plane kinetic normalization

The raw signed embedding was not canonically D-flat. A Newton descent along
noncompact E6 x SL3 uses all86 compact Hermitian generators, with Hessian
4 Re<T_a V,T_b V>. Starting at the actual Cartan T direction it reaches moment
norm below1e-14 in six steps. The certificate stores both finite transformation
matrices, not just a selected norm or abstract equivalence.

After applying these same maps to the **entire** Cartan plane, all basis and
cross moments vanish to numerical precision. The resulting ambient Gram is

    c diag(2^(2/3), 2^(4/3), 1),  c=11.5553820672...,

which is c/3 times the prior exact standard-unitary G26 pullback. Thus the
abstract unitary metric has a numerical full-field D-flat realization.
The finite E6 map preserves the signed cubic and the full mass pencil by
congruence to about1e-14. This is a floating witness, not an algebraic proof
of the displayed transformation. F-flatness, the complete action and the
physical radial scale are additional questions. Earlier reports' metric
compatibility boundary is repaired numerically here, not removed by fiat.

## A mass map with its species specified

Use two distinct chiral fermion multiplets chi and psi and the actual signed
operator M(Phi) from d_E6 tensor epsilon_SL3:

    L_Y = y chi_i^alpha epsilon_alpha_beta M(Phi)_ij psi_j^beta + h.c.

Canonical kinetic terms make its singular values the conditional Dirac mass
ratios. The gauge/matter content must still be completed consistently; this
is not a Standard-Model assignment. A subtle distinction matters: M is skew,
and the underlying trilinear tensor is alternating. Its contraction with
**one commuting scalar field three times vanishes**. It is not a scalar
Hessian or one-field cubic superpotential.

Five real/complex Cartan directions replay the following mirror-spectrum law
to below1e-12. At unit field radius, it is

    12 MUB mirror probabilities, each twice;
    9 SIC probabilities times2/3, each six times;
    3 zero modes.

At T the positive singular values have seven groups, with multiplicities
6,18,6,6,18,6,18 in descending order. The masses are |y|r times these
numbers; neither scale nor particle assignments are predicted. The spectrum
law is a numerical match to the actual normalized field pencil. The subsequent
exact calculations concern this declared mirror mass law.

## Exact low moments: the divergences do not align the vacuum

The SIC and complete MUB sets already exist (Pass11269). Exact arithmetic in
Q(omega) verifies their first/second design operators:

    SIC: sum P = 3I; sum P tensor P = (3/4)(I+Swap),
    MUB: sum P = 4I; sum P tensor P = I+Swap.

The stated mass law consequently has

    sum m² = 20 r²,    sum m⁴ = 8 r⁴.

The quadratic and quartic moments are radial. At fixed radius the subtraction
constant and log(mu) in the one-loop potential are angular constants. Varying
that scale or the common Yukawa magnitude therefore cannot tune the angular
vacuum into existence.

## Outward interval saddle test

For F=sum multiplicity*m⁴ log(m²), the T point has two negative and two
positive projective Hessian directions. Its numerically evaluated eigenvalues
are approximately -3.005465,-3.005465,0.776142,0.776142. Exact-phase tangent
vectors and integer-thousandth Rayleigh witnesses are evaluated with60-digit
outward interval arithmetic; the negative upper bound and positive lower
bound are separated from zero. Reversing the common boson/fermion sign still
leaves a saddle.

Thus this simplest shared mass-loop spectrum does not dynamically produce
the constructed polynomial T minimum. Additional scalar/gauge sectors,
other covariants or UV interactions must be built and tested; they cannot be
assumed to repair it. This does not refute the exact tree-level global selector.

## Why moving the one background cannot lift all three zero modes

Per the exact Pass11260 semisimple certificate, the graded adjoint has a direct
sum of kernel and image in grade one: the commuting Cartan plane has dimension3,
and the G0 orbit tangent has dimension78. The differential of
G0 x Cartan -> g1 is therefore onto. Its image is Zariski dense. Every member
of the Cartan has at least that commuting three-plane in its kernel, and rank
is invariant under G0. The closed condition rank(M)<=78 therefore contains a
dense subset and holds on **all**81-component backgrounds. Rank78 is attained.

This is a consequence of the prior exact Cartan result, not a newly claimed
rank classification. Its mass-model consequence is useful: this one Yukawa
tensor always leaves at least three null modes. A different covariant or
additional fields are needed to lift them; changing Phi alone cannot do it.
No identification of those null modes with three observed families is made.

## An explicit dimension-five light-mass mechanism

A named compact-gauge-invariant nonsupersymmetric EFT operator does lift one
of the null modes:

    (c5/Lambda) (Phi†chi)^alpha epsilon_alpha_beta (Phi†psi)^beta + h.c.

Its matrix is B5=conj(Phi)conj(Phi)^T. For every background, skew symmetry
and M Phi=0 give M†B5=B5†M=0 and B5†B5=||Phi||² Phi Phi†. Thus all old
nonzero singular masses are unchanged, one new mass is exactly
**|c5| r²/Lambda**, and generic rank78 becomes79. Two null modes remain.
The finite witness replays this directly on the actual balanced81-dimensional
operator. This is a controlled hierarchical primitive: the light/heavy ratio
carries r/Lambda. No measured scale, coefficient or particle assignment is
predicted. It is nonholomorphic and cannot be silently treated as a
supersymmetric superpotential or string-allowed coupling.

The new mass depends only on radius, so it does not repair the angular
Coleman–Weinberg saddle. Additional covariants are still needed for the
other two zero modes and for stable quantum alignment.

## The same named mass map can now be a finite gravity input

Define the self-adjoint internal Dirac operator explicitly as
D_F=[[0,A],[A†,0]], where A=yM+(c5/Lambda)B5. Its Hilbert dimension is162,
with each81-dimensional squared singular spectrum duplicated. For the stated
mirror model its moments are

    F2=40|y|²r² + 2|c5|²r⁴/Lambda²,
    F4=16|y|⁴r⁴ + 2|c5|⁴r⁸/Lambda⁴.

Its positive heat factor2 Tr exp(-tA†A) is a fully specified finite input to
the curved isovolume gravity benchmark. This is an operator-level connection,
not an identification of162 with an observed particle count or with the
separate40/480-mode fixtures. The EFT cutoff and the gravitational heat cutoff
remain distinct parameters. It supplies neither spacetime emergence nor the
measured Newton/cosmological scale.

External checks: [Luty–Taylor on complex gauge orbits and vacua](https://arxiv.org/abs/hep-th/9506098),
[Coleman–Weinberg](https://doi.org/10.1103/PhysRevD.7.1888).
Prior owners: Pass11260/61 Cartan, Pass11269 dictionary/design objects,
Pass11267 messenger loop, and the dated exact global selector.
