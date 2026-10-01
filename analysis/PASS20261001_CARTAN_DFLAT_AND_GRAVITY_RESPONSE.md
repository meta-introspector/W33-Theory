# Physical interfaces: D-flat Cartan balancing and a curvature-only gravity response

This extends the global G26 alignment result by checking actual field and
metric operators. It does not turn that constructed potential into a full
81-field action or derive spacetime.

## 1. The Cartan T coordinate is not the literal color qutrit

Producer: `w33_20261001_cartan_dflat_resource_bridge.py`.
Certificate: `data/w33_20261001_cartan_dflat_resource_bridge.json`.

Use the actual Pass 11261 signed 27x3 embedding and declare its canonical
positive kinetic norm Tr(V†V). At the Cartan ray [0:0:1], the column Gram is

    G = [[10,0,0], [0,13,5], [0,5,7]].

The SU3 moment G-Tr(G)I/3 has squared norm68. The eigenvalues are
10, 10-sqrt(34), 10+sqrt(34): no compact column rotation makes them equal.
The exact ambient Cartan Gram is [[59,2,0],[2,98,0],[0,0,30]]; the companion
reconstructed grade-pair Gram is ten times it. This directly connects that
metric candidate to the ambient kinetic choice. It differs from the standard
unitary G26 pullback metric: kinetic compatibility has not been supplied.

Interpreting normalized V as a 27 tensor 3 quantum amplitude gives a mixed
color marginal with purity **92/225**, rather than a pure T ray. More
generally canonical SU3 D-flatness requires V†V proportional to I, so that
interpretation gives I/3 on the color register. The abstract Cartan coordinate
and this external register are different objects.

There is an explicit partial repair:

    V -> V h,  h=det(G)^(1/6) G^(-1/2),  det(h)=1.

This noncompact SL3 transformation balances the column Gram to 660^(1/3)I,
reducing the norm squared from30 to **26.1197630735**. The actual signed cubic
operator transforms by congruence and retains rank78; the numerical covariance
residual is below1e-15. This solves the column moment only: E6 moments and the
full scalar vacuum equations remain open.

The optimal arbitrary one-sided rank-one measurement extracts a chosen pure
color target T with

    p_max=1/[Tr(G) T†G^(-1)T] = 0.265449158913.

Its explicit 27-component covector is stored. For the balanced state the
probability is1/3 for any target. This measurement contains target-dependent
ninth-root phases. Its cost is an input resource, not a free Clifford operation
or a proof of magic distillation. A classical field configuration is not a
prepared quantum state.

Prior ownership: Passes11260/11261, the framed magic contract, and
Pass10962's [Kempf–Ness/D-flat distinction](PASS10962_HIDDEN_SECTOR_DFLAT_TIERS.md).

## 2. A curved Dirac calculation that cannot confuse mass-volume with gravity

Producer: `w33_20261001_isovolume_dirac_gravity_response.py`.
Certificate: `data/w33_20261001_isovolume_dirac_gravity_response.json`.

BT1130 already owns the product heat identity

    C2(g)=N A2(g)-F2 A0(g).

The second term multiplies volume; it is not an Einstein–Hilbert curvature
term. A Ricci-flat background value alone cannot measure the Newton coupling.
Build an explicit off-shell test on an input four-torus of periods2pi:

    g_e=exp(2sigma) g_flat,
    sigma=e*cos(x1)-(1/4)log I0(4e).

The volume stays exactly fixed. The unitarily transported Dirac operator on
flat L2 is Dtilde=exp(-sigma/2)D0exp(-sigma/2), with explicitly stored Fourier
matrix coefficients. Its epsilon-squared heat response is

    K2(t)=-pi²/t + pi²*t/140 - pi²*t²/1260 + O(t³)

up to exponentially small Poisson images. The zero a0 response follows from
fixed volume. The zero a4 response also follows from conformal flatness and
Euler characteristic zero; the Gaussian calculation reproduces both.

Exact rational integration computes the asymptotic coefficients. A separate
finite-epsilon Dirac matrix check agrees to **4.97e-8**. Independent momentum
cutoffs42/46 agree to **9.1e-12** at t=.025. The numerically extracted normalized
curvature coefficient is **0.9999955481**, with the expected O(t²) correction.
These cutoff comparisons are numerical controls, not directed-interval bounds.

For constant Hermitian internal masses the finite heat factor multiplies this
response exactly. The vertex-Laplacian40-mode carrier and the separate
legacy480-mode fixture are explicitly distinct inputs. Their positive factors
retain the conventional Euclidean EH sign corresponding to positive inverse
Newton coupling. Their F2 volume correction cannot fill the curvature response.

This supplies a runnable curved spectral-action benchmark, rather than merely
an arithmetic identification of a finite moment with gravity. The torus,
metric family and spin structure are external inputs. It does not prove a
curved discrete refinement limit, emergent Lorentzian spacetime, the measured
Newton scale, or a resolution of the Euclidean conformal-factor problem.

External checks: [Chamseddine–Connes spectral action](https://arxiv.org/abs/hep-th/9606001),
[conformal Dirac covariance](https://arxiv.org/abs/1409.4983).
Prior ownership: BT1033's geometric route and BT1130's product bookkeeping.
