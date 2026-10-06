# Passes 11580–11589 — Spin(10), Pati–Salam, hypercharge and the external/internal split

This packet resolves the weak-isospin/spacetime conflation exposed by Pass11579 and corrects two older interpretations.

## 11580 — internal Spin(10) produces Pati–Salam and the full one-family charge spectrum

The old 12-dimensional F4∩Spin(9) algebra is exact, but its su(2) rotates the selected event three-plane. It therefore cannot simply be called weak isospin.

Lift instead to the exact Cl(10) system already reconstructed from the Albert Peirce-16. The old u(1) generator on that 16 has primitive spectrum

[
q_{BL}=3(B-L)=(-3)^2+(+3)^2+(-1)^6+(+1)^6.
]

Inside the 45-dimensional executable spin(10), solve the centralizer of su(3)+u(1)_{B-L}. Its dimension is seven and its centroid splits 3+3+1:

[
\boxed{su(2)_L\oplus su(2)_R\oplus u(1)_{B-L}}.
]

On one canonical Spin(10) Weyl 16, normalize the right Cartan by (r=2T_{3R}). The joint weights give

[
6Y=q_{BL}+3r,
]

with multiplicities

[
\boxed{(6Y)=1^6+(-3)^2+(-4)^3+2^3+6+0}.
]

These are exactly (Q,L,u^c,d^c,e^c,\nu^c). No hypercharges were fitted after the decomposition.

Equivalently,

[
Y=T_{3R}+\frac{B-L}{2}.
]

This is the key structural correction: the earlier ±1/±3 u(1) was (3(B-L)), not hypercharge.

## 11581 — the global group is the Z6 quotient and the family is anomaly-free

For (y=6Y), color triality (t\in\{-1,0,1\}), and weak parity (s\in\{0,1\}), every derived multiplet satisfies

[
y+2t+3s\equiv0\pmod 6.
]

The kernel of

[
(g_3,g_2,z)\mapsto(z^{-2}g_3,z^3g_2)
]

has order six, so the native compact group is

[
\boxed{[SU(3)\times SU(2)\times U(1)]/\mathbb Z_6
\cong S(U(3)\times U(2)).}
]

The five perturbative anomaly sums (SU(3)^3,SU(3)^2U(1),SU(2)^2U(1),\mathrm{grav}^2U(1),U(1)^3) all vanish exactly. There are four left weak doublets counting color, so the Witten SU(2) global anomaly is absent as well.

## 11582 — overlap index follows q² once the lattice resolves the charge

On the L=4 four-torus the unit-flux overlap index gives q=1 -> -1 and q=2 -> -4. For q=3 the Wilson gap nearly closes and the coarse lattice loses the expected index. Refining to L=5 restores

[
q=3\mapsto -9,\qquad q=4\mapsto -16,
]

with Ginsparg–Wilson residuals around (4\times10^{-13}).

Thus the index weld behaves like the continuum q² law in the resolved regime. High-charge sectors require an explicit admissibility/refinement check; the coarse result must not be extrapolated.

## 11583 — the clock cycles chiral gradings

The canonical Cl(10) volume element gives an exact chirality involution (Chi) with 16+16 eigenspaces. All 45 spin(10) bivectors commute with it.

The doubled order-eight clock behaves differently. One tick sends (Chi) to a new involution (Chi_1) with

[
{Chi,Chi_1}=0,
qquad
(ChiChi_1)^2=-1.
]

Two ticks send (Chi\to-Chi), and four ticks return it. The full clock therefore cycles a quaternionic-looking pair of chiral gradings rather than preserving a fixed Weyl sector.

This is a firewall: a physical chiral gauge theory must choose or dynamically stabilize a Weyl polarization; the full clock cannot automatically be identified with an internal Spin(10) gauge motion.

## 11584 — a unique single-family 10_H Yukawa channel

Direct Casimir decomposition on the executable Weyl 16 gives

[
\boxed{\mathrm{Sym}^2(16)=10\oplus126},
qquad
\boxed{\Lambda^2(16)=120}.
]

The vector 10 occurs with multiplicity one. Therefore the Spin(10)-invariant (16\,16\,10_H) Yukawa channel is unique up to one overall coefficient for one family.

This does not solve flavor. Once there are several families, a family-space coupling matrix is additional data unless another symmetry fixes it.

## 11585 — coarse variable-frame spectral action is not Einstein–Hilbert alone

Use the Cartan-selected torsion-free metric-compatible connection on the periodic Z3^3 variable frame and construct the spin-connected Wilson/Dirac heat trace at fixed total volume.

For the exact epsilon=0.1 frame,

[
sum_x R(x)=-\frac{360725}{209088},
qquad
\boxed{sum_x \det E(x)R(x)=0}.
]

Nevertheless the heat trace changes already at (O(\epsilon^2)). For example at t=0.1,

[
\Delta H(t)=8.58119\,\epsilon^2-18.74445\,\epsilon^4+\cdots.
]

So a single coarse Hamming cell contains spectral/frame invariants beyond the Einstein–Hilbert scalar term. This does not contradict continuum heat-kernel theory; it says the EH coefficient can only be claimed after a controlled refinement/renormalization limit.

## 11586 — primitive charge lattice

The derived (q_{BL}) and (6Y) lattices both have gcd one. Hence conventional hypercharge is quantized in units of (1/6). Combining (Y) with (T_{3L}) gives

[
Q_{em}=\left\{\frac23,-\frac13,0,-1,-\frac23,\frac13,1,0\right\},
]

for (u,d,\nu,e,u^c,d^c,e^c,\nu^c).

## 11587 — correction: trinification triality is not three generations

The same E6 fundamental 27 has two different branchings:

[
27=16+10+1
]

under Spin(10), and

[
27=9+9+9
]

under trinification. Therefore the three trinification nonets cannot simultaneously be three independent copies of the 16-dimensional family.

The historical repo slogan “the three nonets are the three generations” is superseded. The S3 permutes gauge factors/nonets inside one 27. Three independent E6 families require three 27-plets (81 states before breaking) or another separately derived family mechanism.

## 11588 — external and internal chirality belong on separate factors

The existing four-dimensional overlap background has index -1. Tensor it with the internal Spin(10) Weyl projector of rank 16. All 45 internal spin(10) generators preserve that internal projector and commute identically with external spacetime operators.

This gives the clean architecture

[
H_{physical}=H_{4D\ spacetime}\otimes H_{internal},
]

where external overlap/Cl(4) carries spacetime chirality and the internal Spin(10) Weyl factor carries gauge chirality. The resulting index multiplicity is -16 for one internal family.

This separation is what prevents weak (SU(2)_L) from being confused with spatial/event Spin(3).

## 11589 — synthesis

The strongest current internal chain is

[
\boxed{
Spin(10)
\to Spin(6)\times Spin(4)
\to SU(3)\times SU(2)_L\times SU(2)_R\times U(1)_{B-L}
\to [SU(3)\times SU(2)_L\times U(1)_Y]/\mathbb Z_6.
}
]

One Weyl 16 supplies one anomaly-free SM family plus (
u^c). The 10_H Yukawa channel is unique for one family.

The remaining TOE frontier is now narrower rather than solved: a real family-replication mechanism, flavor/Yukawa texture, vacuum selection, measured couplings, and a controlled variable-frame continuum limit yielding Einstein dynamics remain open.
