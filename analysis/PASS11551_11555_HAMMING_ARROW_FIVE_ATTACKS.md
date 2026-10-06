# Passes 11551–11555 — five attacks from the Hamming-arrow breakthrough

**Status:** all five requested tracks executed with exact finite certificates.

These passes start from Passes 11547–11550 and deliberately separate positive
unifications from negative firewalls.  They do **not** declare the TOE solved:
the continuum, physical action, vacuum selection and measured couplings remain
open unless explicitly derived.

## Pass 11551 — the arrow cubic is an Albert diagonal norm

Pass 11550 identified the temporal arrow on the Hamming null cube with

[
chi(z)=z_1z_2z_3in{pm1}
]

after the usual (1leftrightarrow+1,;2leftrightarrow-1) lift.

The repository already owns an explicit clock Albert algebra in Pass 10950 and
an exact clock/native cubic map in Pass 11471.  Reconstructing that actual
Jordan product and its canonical three primitive idempotents
(e_1,e_2,e_3), Pass 11551 checks all (27) ternary triples and finds

[
oxed{N_{m Albert}(z_1e_1+z_2e_2+z_3e_3)=z_1z_2z_3}
]

**over the integers before reduction modulo 3**.

Thus the history-arrow cubic is not just another degree-three polynomial: it is
the diagonal restriction of the committed Albert determinant/cubic norm.

This is the strongest E6-side weld found in this run.  It is still a
three-dimensional diagonal restriction, not an identification of the entire
27-event history carrier with the full 27-dimensional Albert algebra.

## Pass 11552 — the naive all-null Hamming tower fails as a continuum limit

For (H(n,3)), every nonzero coordinate contributes (1) to
(sum_i d_i^2pmod3), so the finite quadratic-null relation is exactly

[
{m wt}_H(d)equiv0pmod3.
]

The naive (n)-dimensional generalization therefore joins **all**
distance-(3,6,9,ldots) shells.

Using the ternary Krawtchouk spectrum and a roots-of-unity filter, its degree is

[
k_n=rac{3^n+2operatorname{Re}(isqrt3)^n}{3}-1,
]

while every nontrivial character eigenvalue obeys

[
|lambda|le 1+2,3^{n/2-1}.
]

Hence

[
rac{max_{lambda
e k_n}|lambda|}{k_n}	o0,
]

so the normalized absolute spectral gap tends to **one**, exponentially.
That is expander-like behavior, not a diffusive mesh limit.

The distance-three truncation behaves differently.  Its weight-one Fourier
branch has the exact normalized gap

[
oxed{rac{9}{2n}},
]

and this is the ordinary spectral gap for every (nge5) checked through
(n=40).  It therefore supplies a genuine gap-closing family, but increasing
(n) changes the ambient Hamming dimension.  It is not yet a fixed
three-space refinement.

## Pass 11553 — the arrow is the unique PSp-fixed class in (H^1(-;mathbf F_3))

Build the temporal 2-complex from

- 27 events,
- 108 null edges,
- 36 temporal triangles.

Let the arrow edge field be the cubic Hamming sign.

Two different algebraic statements hold simultaneously.

### Chain side

Over the integers,

[
partial_1,chi=0,
]

and because every edge belongs to exactly one temporal triangle,

[
oxed{chi=partial_2 C}
]

for a (pm1) sum (C) of all 36 triangles.

So the arrow is an integral **boundary** as a 1-chain; it is not a new
homology class.

### Cochain side

The oriented circulation around every triangle is (+3) or (-3).  Therefore

[
deltachi=0pmod3.
]

The temporal complex has

[
dim H^1(-;mathbf F_3)=46.
]

After constructing the induced action of generators of
(3^3{:}W(D_3)) on this 46-dimensional quotient, the fixed subspace is

[
oxed{dim H^1(-;mathbf F_3)^{PSp}=1}
]

and the cubic arrow spans it.

So the global arrow has a precise equivariant-topological meaning: it is the
**unique PSp-invariant mod-3 cohomology line**.

The same complex has (H^2(-;mathbf F_3)=0), so this does not by itself
produce a Chern class.  A Maslov or transgression interpretation requires
additional structure.

## Pass 11554 — the tetracode sign lift is genuinely non-split

Pass 10946 had to retain oriented representatives of the four points of
(mathbf P^1(mathbf F_3)) to recover the repository tetracode **exactly**.

The relevant group extension is

[
1longrightarrow{pm I}
longrightarrow GL(2,3)
longrightarrow PGL(2,3)cong S_4
longrightarrow1.
]

Pass 11554 chooses an explicit section (S_4	o GL(2,3)), constructs its
(C_2)-valued factor set, and asks whether that cocycle is a coboundary.

The resulting (mathbf F_2) linear systems have ranks

[
23quad	ext{and}quad24
]

before and after adjoining the cocycle column.  Therefore the section cannot
be gauge-corrected into a homomorphism:

[
oxed{	ext{the extension is non-split}.}
]

This exactly explains why the tetracode representative signs cannot be
removed by a globally (S_4)-equivariant convention.

Also, objectwise,

[
det(A)=operatorname{sgn}(pi_A),
]

so (SL(2,3)) is the preimage of (A_4).

A crucial firewall is now explicit: the central (-I) here sends
(vmapsto-v) and is killed by the Veronese map (vv^{mathsf T}).  It is
**not** the same (C_2) as the later nonsquare symmetric-matrix scaling
(Smapsto-S) that swaps the two temporal/Weil sheets.

External finite-group references also distinguish (GL(2,3)) from the
binary octahedral group: they are two non-isomorphic double covers of (S_4).

## Pass 11555 — the unique local chiral operator and its exact square

Let (A) be the null/distance-three adjacency and define the oriented local
operator

[
D_{xy}=
egin{cases}
chi(y-x),&xsim y,\
0,&	ext{otherwise}.
end{cases}
]

The eight directed null steps form two PSp orbits, the two cubic-sign
tetrahedra.  Therefore every translation- and PSp-invariant nearest-null
kernel is two-dimensional.

Reversal decomposes this space canonically:

- even part: one-dimensional, generated by (A);
- odd part: one-dimensional, generated by (D).

Thus (D) is the **unique** nearest-null PSp-invariant chiral kernel, up to
scale.

With (L=8I-A), direct integer matrix multiplication gives

[
oxed{3D^2=L(L-6I)(L-12I)}.
]

Consequences:

[
operatorname{rank}D=8,qquad dimker D=19,
]

and

[
P_3=-rac{D^2}{27}
=-rac{L(L-6I)(L-12I)}{81}
]

is the exact projector onto the eight-dimensional Hamming-weight-three
Fourier sector.

Define

[
J=-rac{D}{3sqrt3}.
]

Then

[
J^2=-P_3,
]

so (J) is a complex structure exactly on that eight-dimensional sector and
zero on its orthogonal complement.

Fourier transformation gives

[
widehat D(k)=
egin{cases}
0,&{m wt}(k)<3,\
-i,3sqrt3,chi(k),&{m wt}(k)=3.
end{cases}
]

By Pass 11550, the normalized Weil phase is (ichi).  Therefore

[
oxed{J	ext{ acts by exactly the stored Weil phase on all eight null modes}.}
]

The most general local Hermitian operator built from these symmetries is

[
K=m^2I+alpha L+eta,iD.
]

Its sectors are

[
egin{array}{c|c|c}
{m weight}&{m multiplicity}&K	ext{ eigenvalue}\
0&1&m^2\
1&6&m^2+12alpha\
2&12&m^2+6alpha\
3&4+4&m^2+9alpha+3sqrt3,eta,chi
end{array}
]

so (eta) is the unique nearest-null PSp-invariant coefficient that splits
the old (4+4) Weil/chiral sector.

No value of (m,alpha,eta), kinetic term, or continuum scaling is inferred.

## Combined synthesis

The five tracks now give a compact exact chain:

[
oxed{
N_{m Albert}|_{m diagonal}
=
z_1z_2z_3
=
chi
}
]

on the finite null directions,

[
oxed{
[chi]	ext{ spans }H^1(mathcal H;mathbf F_3)^{PSp}
}
]

and the corresponding local skew operator satisfies

[
oxed{
3D^2=L(L-6I)(L-12I),qquad
-rac{D}{3sqrt3}igg|_{m null}
=	ext{Weil phase}.
}
]

Meanwhile two tempting shortcuts are now blocked:

1. the naive full-null Hamming tower is expander-like, not a continuum mesh;
2. the tetracode sign choice is controlled by a non-split (GL(2,3)) double
   cover and cannot be erased by an (S_4)-equivariant convention.

## External prior-art checks

- Standard Hamming-scheme relation eigenvalues are the ternary Krawtchouk
  values; the new content is the application to the repository's null tower.
- The Albert algebra carries the classical cubic norm/determinant; the new
  content is the exact identification with the repository history arrow on
  the committed clock frame.
- (GL(2,3)) and the binary octahedral group are distinct double covers of
  (S_4); the new content is the explicit cocycle obstruction attached to
  the repository tetracode representatives.

## Files

- `analysis/w33_pass11551_hamming_arrow_albert_diagonal_cubic.py`
- `analysis/w33_pass11552_hamming_refinement_continuum_firewall.py`
- `analysis/w33_pass11553_history_arrow_equivariant_cohomology.py`
- `analysis/w33_pass11554_tetracode_orientation_double_cover.py`
- `analysis/w33_pass11555_temporal_local_dynamics_dirac.py`

with matching frozen JSON certificates under `data/`.


## Additional controls from the correlated producer

A second, independently structured producer recomputes all five passes together
and adds three cross-checks that were not assumed in the individual scripts.

### The cohomology line is characteristic-three

The PSp-fixed cohomology dimensions were recomputed over several prime fields:

[
dim H^1(mathcal H;mathbf F_p)^{PSp}
=
egin{cases}
1,&p=3,\
0,&p=2,5,7.
end{cases}
]

Thus the arrow line is not a generic invariant surviving every coefficient
field; it is specifically tied to the ternary substrate.

### The two clock/time double covers are inequivalent

The tetracode lift is the non-split extension
[
2.GL(2,3)	o S_4,
]
whereas the signed-permutation history group decomposes
[
O(3,3)=W(D_3)	imeslangle-Ianglecong S_4	imes C_2.
]

So the code-coordinate orientation lift and the temporal cone-sheet lift are
not merely two gauges for the same cover.

### The order-96 compatibility lift is not the tomotope

Taking both (C_2) lifts simultaneously yields the natural fiber-product
group
[
GL(2,3)	imes C_2
]
of order (96).  Its element-order spectrum was checked against the
repository's archived labelled-tomotope group from Pass 2430 and does not
match.

This closes a tempting (96=96) count-match: the compatibility lift is **not**
the tomotope symmetry group.

Correlated executable:
`analysis/w33_pass11551_11555_hamming_arrow_five_attacks.py`

Correlated certificate:
`data/PART_W33_PASS11551_11555_HAMMING_ARROW_FIVE_ATTACKS.json`
