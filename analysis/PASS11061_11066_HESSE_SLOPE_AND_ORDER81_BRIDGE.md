# Passes 11061–11066 — Hesse-slope finite compiler and split/non-split 81 bridge

This packet links three threads that had remained adjacent: the four Hesse/qutrit clock directions, the ten local triangle factors of the finite cubic weld, and the chamber Sylow-3 group U81.

## 1. The ten local factors are 4 x 2 + 2

The ten order-three Cayley directions supporting the diagonal cubic generator decompose exactly as

```text
P1(F3) x {external +, external -}
        plus
{central +, central -}.
```

So the eight noncentral factors are four Hesse directions with two external slopes each, and the final two are the two pure-center/external FI slopes.

Their symplectic-zero graph is `K2 join (4 K2)`: the two central factors commute at the direction level with everything, while the outer eight occur as four same-Hesse pairs. This explains the previously empirical 21/24 split of factor pairs.

## 2. One triangle foliation can already cover the full 54D quotient

Every individual factor is a rank-54 skew operator made from 27 disjoint weighted 3x3 blocks. Surprisingly, some single factors are already transverse to the compatible S1 sector:

- plus orientation: factors 7 and 8;
- minus orientation: factors 6 and 9.

For these factors the raw rank is 54 and the projection onto S2+L also has rank 54, with S2 and L projections each rank 27. Therefore the factor image maps isomorphically to the complete symmetry-retyped quotient.

Because a factor is skew, its finite Cayley gate preserves that image for every nonzero real parameter. Algebraically this reduces a full retyping transducer to one foliation of 27 parallel three-mode rotations. It is a sparse transducer, not 54 independent control parameters and not the same finite gate as the full diagonal weld.

## 3. Five redundant two-slope routes

Pairing the two branches over each Hesse direction and pairing the two central branches gives five independent route families. Every route projects with rank 54 onto S2+L for both FI orientations.

The raw-rank pattern is

```text
66, 54, 54, 72, 72
```

and hence the S1-overlap pattern is

```text
12, 0, 0, 18, 18.
```

Two Hesse pair-routes are exact transverse 54-planes.

## 4. The outer eight are tetracode-indexed

Choose the four projective Hesse points in the order

```text
(1,0), (0,1), (1,1), (1,-1).
```

Evaluating a nonzero covector on these points gives the ternary tetracode word

```text
(a,b,a+b,a-b).
```

Each projective direction is the unique zero coordinate of exactly one opposite word pair. Matching the two external branches above that direction with the two signs gives an exact bijection

```text
8 noncentral weld factors <-> 8 nonzero tetracode words
```

after the displayed coordinate/sign convention. The two central FI factors remain outside this tetracode eight.

This is an indexing theorem, not an assertion that factor matrices are the E8 tensor sectors.

## 5. There are two fundamentally different order-81 groups

The scheduler group used by the finite compiler is

```text
K81 = H27 x C3_external.
```

It has center 9, derived subgroup 3, class 2 and exponent 3. It is the split extension

```text
1 -> C3_external -> K81 -> H27 -> 1.
```

The chamber group U81 has center 3, derived subgroup 9, lower-central sizes

```text
81 -> 9 -> 3 -> 1,
```

class 3 and exponent 9. Quotienting by its center gives a nonabelian exponent-3 group of order 27 with center and derived subgroup of order 3, hence H27. Thus

```text
1 -> C3 -> U81 -> H27 -> 1
```

is a non-split central extension.

So K81 and U81 are not two gauges for the same group. They are the split and non-split 81-lifts of the same Heisenberg memory group.

## 6. The exact Jennings signature of the weld

Over F3,

```text
Hilb(K81) = (1+t+t^2)^3 (1+t^2+t^4)
degrees    = 1,1,1,2
```

while the previously certified chamber controller has

```text
Hilb(U81) = (1+t+t^2)^2 (1+t^2+t^4) (1+t^3+t^6)
degrees    = 1,1,2,3.
```

Passing from the split scheduler to the chamber group therefore replaces exactly one degree-1 ternary factor by the degree-3 highest-root factor. The total module dimension remains 81, but the last nonzero augmentation power moves from 10 to 14.

That is the sharp algebraic version of the intuitive statement that the external ternary scheduler coordinate becomes a class-three central chamber coordinate. Jennings degree is algebraic memory depth, not physical elapsed time.

## Current frontier

The strongest next target is now explicit rather than numerical: construct a concrete intertwining/deformation mechanism from the split K81 carrier to the non-split U81 chamber extension, and determine whether the single-foliation 54D transducers can be embedded into that class-three update law without losing the certified S2+L transversality.
