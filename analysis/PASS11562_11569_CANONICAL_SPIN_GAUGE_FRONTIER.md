# Passes 11562-11569 — canonical spin / connection / cohomology / exceptional gauge frontier

This packet executes the five post-11561 frontier attacks plus three additional probes.  It contains both constructions and no-go results; the latter are frozen because they remove genuine ambiguity from the TOE program.

## 11562 — the Albert event Cl3 is valid, but not canonical

The committed Peirce-16 carries exact Cl(9).  Any chosen three orthonormal gamma directions generate an 8-dimensional Cl3 algebra.  For the Pass11558 choice:

- Spin(9) has dimension 36.
- The stabilizer of the event three-plane is Spin(3)xSpin(6), dimension 18.
- The centralizer of the internal event Spin(3) is Spin(6), dimension 15.
- There are 84 coordinate three-planes, and Spin(9) acts transitively on oriented orthonormal three-planes.

Therefore Cl(9)/Spin(9) alone cannot canonically select the event Cl3.  A new selector is still required.

## 11563 — the spin connection is kinematically unique

Let E^a_i be a nondegenerate 3x3 coframe and impose metric compatibility omega_i^{ab}=-omega_i^{ba}.  The discrete torsion-free Cartan equations

C^a_ij + omega_i^a_b E^b_j - omega_j^a_b E^b_i = 0

form a 9x9 linear system.  Exact symbolic elimination gives

det A(E) = -2 (det E)^3.

Hence every invertible frame determines a unique local metric-compatible torsion-free spin connection.

This closes the kinematic gap identified in Pass11559.  It does not supply an Einstein-Hilbert action or frame dynamics.

## 11564 — the arrow is not an affine-group character

For G=3^3:W(D3) with trivial F3 coefficients, the low-degree LHS data give

H^1(W(D3),F3)=0,
H^1(3^3,F3)^{W(D3)}=0,

so H^1(G,F3)=0.

For the four-direction permutation module M=F3^4, direct cocycle linear algebra gives

dim Z^1(W,M)=3,
dim B^1(W,M)=3,
H^1(W,M)=0,

while H^0(W,M) is the unique all-ones line.

Thus the temporal arrow is precisely a residual coefficient-module invariant, not a hidden group homomorphism of affine PSp.

## 11565 — the continuous Hesse stabilizer inside E6 is exactly F4

Using the executable Albert multiplication:

- all derivation commutators annihilate the trace covector;
- the traceless multiplication sector has dimension 26;
- the trace covector plus its 26 traceless-multiplication images span the full 27-dimensional dual space.

Therefore no nonzero traceless multiplication preserves the trace line.  Since E6(-26)=Der(J)+L(J_0),

Stab_E6(<T>) = Der(J) = F4,

dimension 52.

Because E6 already preserves the Albert norm N, the infinitesimal stabilizer of the Hesse polynomial-law two-plane <N,T^3> is also exactly F4.

## 11566 — canonical axis Dirac and Hamming-Wilson regulator

A stronger exact factorization exists than the null-edge B_i construction:

delta_i = T_i - T_i^{-1}

on the weight-1 Hamming orbital, and at N=3

D = delta_1 delta_2 delta_3.

Define

Q_axis = i sum_i Gamma_i delta_i.

Then exactly

Q_axis^2 = I4 tensor (6I-A1),

so the weight-1 Hamming Laplacian is simultaneously the square of the first-order Dirac operator.

This also corrects the global refinement picture.  The body-diagonal B_i used in Pass11557/11559 have Fourier symbol

8 i sin(k_i) product_{j!=i} cos(k_j),

which has 12 nodal lines in addition to corner zeros.  Pass11559 remains correct locally near k=0, but the globally clean refinement is the axis Dirac.

The same weight-1 Laplacian is a canonical Wilson operator:

H_W = Q_axis/(2a) + Gamma_0 r W1/a,
W1=6I-A1,

with exact cross-term cancellation.  For r != 0 and zero bare mass, its only Brillouin-zone zero is the origin.

## 11567 — heat trace counts 8 species before Wilson, 1 after

On the periodic 2pi three-torus with a=2pi/N, compare the lattice heat trace with the one-spinor continuum target

4 [sum_{n in Z} exp(-t n^2)]^3.

At N=81 and t=1:

- unregulated axis Dirac / continuum = 8.0339119054;
- Hamming-Wilson / continuum = 0.9825851383.

The first converges to eight lattice species; the second converges to one.

## 11568 — a genuine integer winding invariant

For

h(k)=(sin k1,sin k2,sin k3,m+2 sum_i(1-cos ki)),

the normalized map T^3 -> S^3 is defined away from m=0,-4,-8,-12.  Its exact corner formula is

nu = -1/2 sum_{n=0}^3 C(3,n)(-1)^n sign(m+4n).

The five phases are

0, +1, -2, +1, 0.

Direct Brillouin-zone integration reproduces those integers to better than 2e-6.

This upgrades the old 4+4 chirality language to a real lattice topological invariant, while not yet identifying it with an observed chiral index.

## 11569 — the present exceptional selectors stop at Spin(6), not the SM

The exact selector chain is now

E6 --Hesse trace--> F4 (52)
   --primitive idempotent--> Spin(9) (36)
   --event three-plane--> Spin(3)xSpin(6) (18)
   --commute with event Spin(3)--> Spin(6) ~= SU(4) (15).

The Standard Model gauge algebra has dimension 12.  Therefore the new Hesse+Peirce+event selectors do not yet derive it.

The repository already quotes the Todorov-Dubois-Violette target
S(U(2)xU(3)) = Spin(9) intersect (SU(3)xSU(3))/Z3 inside F4,
but that second SU(3)xSU(3) selector is not yet executable in this chain.  It remains the missing exact weld.

## Physical boundary

These passes derive finite Clifford/Jordan/cohomological structure, a kinematically unique discrete spin connection, a regulator-clean continuum Dirac candidate, heat-trace species reduction, and a lattice winding invariant.  They do not derive Einstein dynamics, physical Lorentzian signature, the Standard Model gauge group, observed chirality, masses/couplings, or vacuum selection.
