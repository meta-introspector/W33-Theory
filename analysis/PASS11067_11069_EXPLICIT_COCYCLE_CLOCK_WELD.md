# Passes 11067–11069 — explicit cocycle clock weld

The split/non-split comparison of Pass 11065 can now be written as one exact multiplication law on the same 81 labels.

## Exact root-coordinate law

Use the C2 normal form

```text
u(a,b,c,d)=x0(a)x1(b)x2(c)x3(d).
```

Symbolic matrix multiplication gives

```text
u(a,b,c,d) u(A,B,C,D)
 =
u(a+A,
  b+B,
  c+C-A b,
  d+D+A^2 b-2 A c).
```

The first three coordinates are exactly the frozen H27 multiplication law. Therefore the chamber group is a central C3 extension of H27 with cocycle

```text
kappa((a,b,c),(A,B,C)) = A^2 b - 2 A c  mod 3.
```

The cocycle identity was checked on all 27^3 triples. A direct coboundary linear system has coefficient rank 25 and augmented rank 26 over F3, so this cocycle is not a coboundary.

Hence K81 and U81 can be placed on the same set H27 x F3:

```text
K81: (h,d)(h',D) = (hh', d+D)

U81: (h,d)(h',D) = (hh', d+D+kappa(h,h')).
```

This is the explicit algebraic weld that was missing from the split/non-split theorem.

## A ternary deformation parameter

Scale the cocycle by s in F3:

```text
(h,d) star_s (h',D)
  = (hh', d+D+s kappa(h,h')).
```

Then:

- s=0 is the split scheduler K81;
- s=1 is the class-three chamber U81;
- s=2=-1 is the same non-split extension with reversed central orientation.

The two nonzero laws are isomorphic by

```text
(h,d) -> (h,-d).
```

Their invariants are exactly the U81 invariants: center 3, derived 9, lower-central sizes 81->9->3->1, and 36 elements of order 9. The split member has center 9, derived 3 and exponent 3.

Thus the split/nonsplit jump is controlled by one discrete cocycle coefficient.

## Primitive -> area -> central tick

Let

```text
e0=(1,0,0,0),
e1=(0,1,0,0).
```

The exact commutator ladder is

```text
[e0,e1]_s = (0,0,1,s),

[e0,[e0,e1]_s]_s = (0,0,0,-s),

[e0,[e0,[e0,e1]]] = 1.
```

At s=0, the first commutator produces only the H27 phase/area coordinate and the next bracket vanishes.

At s=+/-1, the second commutator produces a pure highest-root central coordinate, with sign reversed when the cocycle orientation reverses.

This is exactly the group-theoretic class jump from two-step Heisenberg memory to the three-step chamber update law. It also matches the independent Jennings signature of Pass 11066: degrees 1,1,1,2 become 1,1,2,3.

The phrase "central tick" is structural. The theorem is finite algebra; a physical arrow of time would still require a dynamical reason to select a nonzero cocycle orientation.
