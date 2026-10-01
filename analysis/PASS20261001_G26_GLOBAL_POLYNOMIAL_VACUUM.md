# A regular polynomial potential globally selects the T-magic orbit

Producer: `w33_20261001_g26_global_polynomial_vacuum.py`.
Certificate: `data/w33_20261001_g26_global_polynomial_vacuum.json`.

Pass 11256 already owns controlled G26 invariant potentials; the balanced
barrier already owns the exact local T-magic minimum. This extension removes
logarithmic singularities and upgrades **local selection to exact global
selection in a constructed potential**, without fitting an angular target.

For positive lambda, alpha, beta and r0, take

    V(v) = lambda (v†v-r0²)² + alpha |u6(v)|² + beta |u12(v)|².

Every term is nonnegative. The T ray has u6=u12=0, so the global value is zero.
The entire global minimizer set consists of exactly 72 projective rays at
norm r0, with an unfixed common phase. Their location does not depend on the
positive coefficients. They are the known projective G26/Clifford T orbit.

## Exact complete-intersection proof

In the chart z=1 put A=x³ and B=y³. The exact Groebner basis is linear in B
and has a square-free degree-eight eliminant in A:

    (A²+A+1)(A³-3A²-24A-1)(A³+24A²+3A-1).

The coefficient of B in the other basis element is the nonzero constant
1539. Exact coprimality checks exclude coordinate-zero points and points at
infinity. Thus there are eight simple cube-ratio solutions, each with nine
unramified cube-root lifts: **8×9=72 reduced projective points**.

No common zero has u18=0. Scale any zero to the T value of u18; the prior
basic invariant ring and finite-group orbit separation then identify it
with the T orbit. This argument does not use a numerical orbit count.
Transversality gives positive projective Hessian
2 Re(J† diag(alpha,beta) J) at every zero.

The exact gradient Gram at unit T is diag(16,100). In the declared Euclidean
real kinetic coordinates the full Hessian spectrum at radius r0 is

    0, 8 lambda r0², 32 alpha r0^10 (twice), 200 beta r0^22 (twice).

The zero is the common phase. These doublet curvatures are exact model
fluctuations; the unspecified couplings prevent interpreting them as measured
particle masses.

## What this supplies and what it leaves open

This is a regular scalar effective-potential candidate on the exact three-complex-dimensional Cartan
dictionary. An explicit invariant extension and its couplings to the other
78 grade-one directions have not been built. Its degree-12 and degree-24 interactions need a cutoff or UV
completion in four dimensions. The positive couplings, radius and kinetic
normalization are inputs: angular vacuum location is fixed, physical masses
are not. The common phase must be gauged or lifted. The selected orbit does
not choose CP/time orientation, and a scalar vacuum is not a prepared qutrit
or a measured magic-state factory. The zero of this constructed potential
does not predict the cosmological constant; an additive constant and quantum
corrections remain unconstrained.

External mathematical checks: [Milne, Elliptic Curves, Theorem 1.18](https://www.jmilne.org/math/Books/ectext6.pdf) states the projective intersection count; [finite-group invariant theory](https://people.maths.ox.ac.uk/drutu/tcc6/onishchik-vinberg.pdf) supplies orbit separation. The certificate gives the explicit elimination proof rather than relying on the degree count alone.
