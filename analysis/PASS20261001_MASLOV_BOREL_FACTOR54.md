# 2026-10-01 — Maslov–Borel Sylow orbit-defect theorem

## Result

The remaining factor two in the 256-pair Maslov chirality is now explained
inside the exact finite symmetry already present in the W33 chamber geometry.

For either member of the time-reversed 256 pair,

\[
\boxed{\chi=\pm54=\pm(81-27).}
\]

The 81 is a free orbit of the chamber Sylow-3 group \(U_{81}\). The 27 is an
orbit with stabilizer a noncentral \(C_3\subset U_{81}\). Equivalently,

\[
\frac{|\chi|}{27}=3-1=2.
\]

This resolves the explicit open boundary left by Pass 11222.
## Effective relation stabilizer

Pass 11222 used a matrix stabilizer \(\Gamma\) of order 324. Its action is
projective on Lagrangian subspaces, so the scalar \(-I_6\) is ineffective.
The new verifier constructs the quotient explicitly:

\[
|\Gamma|=324,\qquad
|\bar\Gamma|=|\Gamma/\{\pm I\}|=162.
\]

A six-generator pure-Python certificate gives a bijective homomorphism

\[
\boxed{\bar\Gamma\cong N_{\mathrm{PSp}(4,3)}(U_{81})
      =U_{81}\rtimes C_2.}
\]

Independent GAP identification gives both groups SmallGroup [162,10] with
structure ((C3 x C3 x C3) : C3) : C2.
The raw matrix group is SmallGroup [324,68], and GAP independently verifies

\[
\Gamma\cong C_2\times\bar\Gamma.
\]

Thus the order 324 used in the earlier orbit computation contains one
projectively invisible scalar factor two.

## The false order-324 weld

A tempting first guess was

\[
\Gamma\stackrel{?}{\cong}N_{\mathrm{PGSp}(4,3)}(U_{81}),
\]

because both groups have order 324. It is false.

Their exact invariants separate them: \(\Gamma\) has center 6, derived
subgroup 27 and 36 elements of order 18, whereas
\(N_{\mathrm{PGSp}}(U_{81})\) has center 1, derived subgroup 81 and no
elements of order 18.
The correct identification appears only after quotienting the ineffective
\(\{\pm I\}\) from \(\Gamma\), and it lands in the symplectic flag Borel of
order 162, not the similitude normalizer of order 324.

## Restriction to the Sylow-3 core

For the \(Q=-1\) member, the Borel orbit imbalance from Pass 11222 is

\[
\Delta n_{27}=+1,\qquad
\Delta n_{81}=+1,\qquad
\Delta n_{162}=-1.
\]

Restricting those orbits to \(U_{81}\) gives only orbit sizes 27 and 81.
Every 54-orbit splits as \(27+27\); every 162-orbit splits as \(81+81\).
The imbalance simplifies to

\[
\boxed{\Delta n_{27}=+1,\qquad\Delta n_{81}=-1,}
\]

hence \(\chi=27-81=-54\). For the time-reversed \(Q=+1\) member the signs
reverse and \(\chi=+54\).
The 27-orbit stabilizers were checked objectwise. They all have order three,
but none is the center \(Z(U_{81})\); they are noncentral \(C_3\) subgroups.

Therefore the structural identity is

\[
\boxed{
|\chi|
=
|U_{81}|-|U_{81}/C_3^{\mathrm{noncentral}}|
=
81-27
=
54.
}
\]

The residual Borel \(C_2\) is not the chirality source. It merely fuses pairs
of \(U_{81}\)-orbits.

## Evidence and boundary

Machine certificate: data/w33_20261001_maslov_borel_factor54.json.

Executable producer: analysis/w33_20261001_maslov_borel_factor54.py.

Parents are Passes 11070, 11084, 11212 and 11222. The new increment is the
projective-Borel identification and the Sylow-3 orbit refinement that removes
the final unexplained factor two.

This is a theorem about one finite Maslov/Bargmann chirality count. It does
not make 54 a universal physical constant, elapsed-time rate, coupling,
particle mass, or continuum prediction.
