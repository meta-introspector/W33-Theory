# 2026-10-01 — the balanced G26 chamber center is a qutrit T-magic orbit

Producer: `analysis/w33_20261001_balanced_g26_tmagic_vacuum.py`

Certificate: `data/w33_20261001_balanced_g26_tmagic_vacuum.json`

The mirror-separation theorem gives two products on the same 21-wall
arrangement: A for the nine SIC mirrors and B for the twelve
stabilizer/MUB mirrors, with (A B)^5 proportional to P J^2.

This makes the coefficient-free equal-wall projective barrier natural:

    V_bal(v) =
      - sum_(9 SIC) log |<phi|v>|^2
      - sum_(12 MUB) log |<s|v>|^2,     ||v|| = 1.

It is a repulsive chamber-center functional, not a fitted scalar potential.

## Exact local theorem

At the qutrit T state

    |T3> = (1, zeta_9, zeta_9^-1) / sqrt(3),

the stationarity equation is exact in Q(zeta_36):

    sum_n n <n|T3> / |<n|T3>|^2 = 7 T3.
In an exact four-real-dimensional projective tangent basis the Hessian is

    [ 42   0 -12   0 ]
    [  0  42   0  12 ]
    [-12   0  42   0 ]
    [  0  12   0  42 ],

so its spectrum is exactly (30, 30, 54, 54). The point is a strict local
minimum of the equal-wall barrier.

The unit-state basic invariant coordinates are also exact:

    u6 = 0,     u12 = 0,     u18 = -1/729.

The standard unitary G26 reflection representation supplies a canonical
Hermitian form up to scale. Pulling it through the exact Pass 11261 coordinate
map gives

    3 diag(a^2, a^4, 1),     a^3 = 2,

or 3 I in rescaled coordinates (a x, a^2 y, z).

This Hermitian metric is distinct from the real grade-pair metric extracted
directly from the E8 bracket in the companion Cartan certificate. They encode
different additional structures; symmetry alone does not identify them.
## Finite orbit and the TM1 seed

The projective qutrit Clifford group has order 216. The T3 ray has stabilizer
order 3 and orbit size 72. This 72-ray resource orbit is not new to the repo:
Pass 416 already exhausted it as the input/output correction family for the
five-qutrit distillation search. The qutrit T-port also uses exactly this
T_PLUS ray as its Hesse-SIC factory-audit reference. The increment here is
that the G26 mirror chamber selects that known computational resource orbit
dynamically. The port still requires a separate entangled T-Choi resource for
actual gate injection; selecting T_PLUS does not manufacture that Choi pair.

Every one of the 72 rays selects exactly one of the four qutrit MUBs by the
purity pattern

    (1/3, 5/9, 5/9, 5/9),

and each MUB is selected by exactly 18 rays.

Pass 11262 identifies those four MUBs with the four tetrahedral phi axes.
Thus the T-magic chamber orbit canonically chooses one tetrahedral vertex.
Its residual C3 stabilizer cycles the remaining three vertices. A further
ternary selector therefore chooses an ordered MUB pair, equivalently one of
the incident oriented tetrahedral chords; Pass 11267 classifies every such
incident chord as TM1.

## Chirality firewall

Complex conjugation of T3 is a Clifford permutation of the same orbit and has
exactly the same G26 quotient coordinates. Consequently no scalar
G26-invariant potential on the Cartan quotient can select which conjugate/time
orientation occurs.
The fixed cubic/SUM postselection interface makes the same point
operationally: CP is not constant on the 72-ray vacuum orbit. A secondary
frame/interface choice is required. The companion U81 certificate supplies an
exact odd variable — its two cocycle orientations are exchanged by qutrit
conjugation — but an energetic coupling selecting one sign remains open.

The recurring cyclotomic level

    ((1 + 2 cos(2 pi / 9)) / 3)^2 = 0.712386014201...

appears both as the largest T3 MUB probability and as the exact two-qutrit
reversal-fidelity level. This is a shared algebraic invariant, not yet a
derived physical identity.

## Polynomial-wall no-go

The lowest simple homogeneous polynomial that can weight the two separated
mirror products independently is

    V24 = lambda_S ||v||^6 |A|^2 + lambda_M |B|^2.

For nonnegative coefficients it has zero-energy SIC-intersection-MUB lines,
so it drives the field onto the discriminant rather than selecting a regular
vacuum. The repulsive logarithmic barrier avoids this failure mode.

## Boundary

The stationary point, quotient coordinates and Hessian are exact. The
72-ray orbit/MUB census is a finite numerical group enumeration. Global
optimality of the logarithmic barrier is not proved. No Standard-Model mass,
mixing angle, physical radial scale, or energetic chirality choice is derived.
