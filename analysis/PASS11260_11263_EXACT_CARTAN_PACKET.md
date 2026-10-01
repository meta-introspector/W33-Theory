# Passes 11260–11263 — the signed E8 Cartan plane is explicit

This packet closes the embedding blocker left by Pass 11255 and follows the
result through the G26/qutrit coordinates, the line-stabilizer flavour layer,
and the complete 248-mode diagnostic. It also records why neighboring Pass 11264
remains deferred rather than silently inheriting a claim.

## 11260 — an exact semisimple Cartan three-plane

In the repository's signed `3 tensor 27` gauge, the twelve-support vector

```
(3,-1) (7,-2) (13,-2) (14,-1) (22,2) (23,1)
(27,1) (42,-2) (65,2) (72,2) (79,1) (80,1)
```

has exact cubic-Jacobian rank 78. Its centralizer in grade one is the displayed
three-dimensional commuting plane in the certificate. For `N=18 ad(v)`, exact
integer elimination gives `rank(N)=rank(N^3)=234`. The three diagonal blocks of
`N^3` have ranks `78,78,78` and the same nonzero characteristic polynomial.
The only repeated factors are two square-free degree-nine polynomials; each has
its full 27-dimensional kernel in each block. Hence `N^3`, then `N`, is
diagonalizable. The centralizer has the minimal Vinberg dimension three, so it
is a genuine Cartan subspace.

This corrects the *object* used in Pass 11218 without restoring its retracted
nilpotent plane.

## 11261 — the exact standard G26 coordinate map

Let `(x,y,z)` multiply the new integer Cartan basis. The fixed complementary
Pfaffian is

```
2^17 C3 N9 S9^3,
```

with factor degrees and multiplicities `(3,1),(9,1),(9,3)`. The identity is
checked on all 820 points of a degree-39 unisolvent grid.

Put `a^3=2`, let `t` be a primitive eighteenth root
(`t^6-t^3+1=0`), and set `omega=t^6`. A scaled ternary Fourier transform,
written explicitly in the certificate, obeys

```
X Y Z = C3,
u9(X,Y,Z) = 3(omega-omega^2) S9,
product(12 standard stabilizer mirrors) = -27 C3 N9.
```

Thus the E8 restricted Pfaffian and Pass 11269's G26 Jacobian live on exactly
the same 21 hyperplanes. They are different relative divisors: the Pfaffian
weights the 12 stabilizer mirrors once and the 9 SIC mirrors three times; the
Jacobian weights them two and one times. The degrees 6,12,18 invariants now
pull back to the signed E8 tensor by an explicit linear substitution.

## 11262 — why the flavour group is S4 only after time reversal

The twelve stabilizer mirrors split into the four qutrit MUBs. The complex-
linear G26 reflections induce exactly `A4` on these four bases. Complex
conjugation induces the odd transposition `(2 3)`; adjoining it closes the
permutation image to `S4`.

This supplies an objectwise dictionary for the line-stabilizer model:

- the four tetrahedral flavon axes label the four MUBs;
- the twelve vacua `±(d_i-d_j)` label the twelve oriented MUB pairs;
- odd MUB permutations, and therefore the full line-stabilizer `S4`, require
  the antiunitary completion.

The result also explains the scope found in Pass 11254. Epsilon controls the
breaking of the TM1 column, while theta13 and delta retain independent Yukawa
data. Pass 11267's bosonic messenger naturally produces the required
`|epsilon|` window; it does not turn the angles into functions of epsilon alone.

The parallel Pass 11266 magic-gate count sharpens the same dictionary from the
computational side. For one qutrit, all 18 one-gate arrow violations lie in the
unit-shear parabolic that fixes one MUB axis and translates transversely; its
parity-twisted coset has no violations. Under the four-MUB tetrahedron this
localizes the minimal time-reversal obstruction at an *oriented stabilizer
stratum*. This is an exact classification of the 216 Clifford prefixes, not yet
a theorem for two or more qutrits: Pass 11266 explicitly finds that the naive
fixed-axis rule fails there.

## 11263 — 78 exact Z3 spectral triplets

With the stored E8 coordinates declared Euclidean-orthonormal, the positive
quadratic bracket energy

```
V(delta)=1/2 ||[v,delta]||^2
```

has Hessian `(18 ad(v))^T(18 ad(v))/18^2`. Its complete 248 eigenvalues are in
the certificate. The `g0,g1,g2` blocks each have rank 78 and respectively
`8,3,3` zero modes. In particular, the full matter-81 block has precisely the
three Cartan directions as its kernel.

The basis-independent statement is sharper: `(18 ad(v))^3` has the same 78
nonzero eigenvalues in every grade. All 234 nonzero modes therefore form 78
exact Z3 spectral triplets, and every positive-power root-of-unity graded trace
cancels. The residual graded index is five, entirely in the zero modes:
`8+3 omega+3 omega^2=5`. No ordinary boson/fermion signing of the three grade
blocks cancels moments zero through six. The Z3 identity is not promoted to a
cosmological-constant calculation.

## 11264 — Z6-I R charges: input restored, exhaustive audit deferred

The raw 87-model R-charge archive is present and passes the frozen shape/hash
checks (30,980 fields). The producer intentionally requires either a cached
exhaustive run or an explicit `--run`; no cached exhaustive result was found in
this checkout, and the producer estimates tens of CPU-hours for a fresh run.

Therefore **Pass 11264 is not claimed in this publication packet**. The restored
input and producer remain local until the full 6,695,116-vacuum audit is frozen
as a certificate. Passes 11260–11263 are independent of this deferred front.

## Evidence

- `analysis/w33_pass11260_exact_semisimple_cartan.py`
- `analysis/w33_pass11261_g26_cartan_coordinate_map.py`
- `analysis/w33_pass11262_g26_mub_s4_yukawa_bridge.py`
- `analysis/w33_pass11263_full_graded_spectrum.py`
- `data/w33_pass11260_exact_semisimple_cartan.json`
- `data/w33_pass11261_g26_cartan_coordinate_map.json`
- `data/w33_pass11262_g26_mub_s4_yukawa_bridge.json`
- `data/w33_pass11263_full_graded_spectrum.json`

External checks used: Reeder–Levy–Yu–Gross, *Gradings of positive rank on
simple Lie algebras*, Table 21; Nilles–Ramos-Sánchez–Vaudrevange–Wingerter,
*The Orbifolder* and its public 1.2.1 implementation.
