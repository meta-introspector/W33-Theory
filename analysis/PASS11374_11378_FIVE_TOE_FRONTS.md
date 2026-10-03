# Passes 11374–11378: five bounded advances on the physical TOE fronts

Producer: `analysis/w33_pass11374_11378_five_frontier.py`. Certificate: `data/w33_pass11374_11378_five_frontier.json`. The five sections below are exact within their stated ansatz or kinematic model. None selects observed masses, the Standard Model vacuum, gravity dynamics, or vacuum energy.

## 11374 — a CP-even selector for the existing three-ray flavor bridge

Use the prior rays $u=(1,0,0)$, $v=(1,1,0)/\sqrt2$, and the one-parameter continuation $w_\phi=(1,e^{i\phi},1)/\sqrt3$. Their gauge-invariant Bargmann loop and Gram determinant are

$$
\mathcal B(\phi)=\operatorname{Tr}(P_uP_vP_{w_\phi})=(1+e^{i\phi})/6,
\qquad \det(X^\dagger X)=1/6.
$$

With the positive weights $a=(1,2,4)$, $b=(2,3,5)$ of Pass 11361, the exact identity becomes

$$
\frac{\operatorname{Tr}[A,B]^3}{i}=\sin\phi,
\quad A=X\operatorname{diag}(a)X^\dagger,
\quad B=X\operatorname{diag}(b)X^\dagger.
$$

The **supplied** CP-even and ray-rephasing-invariant operator

$$
V_\phi=\lambda\bigl(\operatorname{Re}\mathcal B-1/6\bigr)^2
=\frac{\lambda}{36}\cos^2\phi,\qquad\lambda>0,
$$

selects exactly two minima modulo $2\pi$: $\phi=\pi/2,3\pi/2$. Both have positive phase curvature $\lambda/18$ and opposite nonzero commutator cubics $+1,-1$. The selected positive branch has exact Gram characteristic polynomials $t^3-7t^2+9t-4/3$ and $t^3-10t^2+59t/3-5$, with discriminants $2063/3$ and $142459/27$. Its mixing-invariant magnitude is exactly $1/(6\sqrt{(2063/3)(142459/27)})\approx8.74977\times10^{-5}$. These values follow from the supplied weights and are not observed mass inputs. Thus this particular ray phase can be selected **spontaneously** without an explicit CP-odd coefficient. CP cannot choose between the two branches in this action. The ansatz, weights and operator are supplied, and the full ray-field normal Hessian is not established. This is a bridge to the previously owned [Pass 11286 CP-even sextet toy selector](w33_pass11286_misaligned_family_vacua.py), not a claim that CP-even toy selection was absent. The physical interpretation of the commutator follows [Jarlskog's basis-independent analysis](https://cds.cern.ch/record/161795).

## 11375 — why any number of norm-only heavy thresholds cannot fill the hard block

Pass 11344 found 105 free entries in the symmetric $14\times14$ normal Hessian; Pass 11362 constructed one radial hard loop direction. This generalizes that calculation. Let the *total* renormalized zero-momentum heavy potential be $f(I)$, with the same invariant $I=\phi^\dagger\phi$ for every species. At a stationary point $I_0$, $f'(I_0)=0$, so the chain rule gives

$$
\partial_a\partial_b f(I)|_{I_0}=f''(I_0)(\partial_a I)(\partial_b I).
$$

This has rank at most one, regardless of the number of norm-only species, their masses or counterterms. This rank bound requires the norm-only sector's *own* total tadpole to vanish. If a different invariant sector cancels a nonzero $f'(I_0)$, the $f'(I_0)\,\partial_a\partial_b I$ term remains and the rank-one conclusion does not follow. Under the stated separate-tadpole hypothesis, at least 13 directions of the 14-dimensional quotient normal space remain untouched by this mechanism. It cannot determine the remaining 104 independent hard entries. For a rank-one-source construction to span *all* symmetric normal matrices, at least 105 independently controllable directions are needed: the abstract vectors $e_i$ and $e_i+e_j$ for $i<j$ give exactly 105 outer products spanning the block. No global E6-invariant realization of that abstract basis is supplied. A zero-momentum potential Hessian is also not the momentum-dependent pole matrix. See the [one-loop derivative framework](https://arxiv.org/abs/1606.07069) for the distinction.

## 11376 — exact kinematic count for native flat vertex-frame transports

The actual connected W33 Levi graph has 80 vertices, 160 edges and 81 independent cycles (replayed from Pass 11289's integral fundamental-cycle basis). Give each edge an independent $SO^+(1,3)$ transport: there are $160\cdot6=960$ local connection coordinates. If every transport is required to be $\Lambda_i^{-1}\Lambda_j$ for vertex frames, a spanning tree gauges away its 79 edges, leaving the 81 chords as independent fundamental holonomies. Requiring each to be the identity imposes $81\cdot6=486$ independent equations near the identity. The flat solution has $79\cdot6=474$ coordinates before vertex gauge and **zero after vertex gauge**.

This is an exact graph-connection statement, not a gravity degree-of-freedom count. It sharpens Pass 11367: adding independent transverse edge boosts without their compensating rotations creates holonomy, while an actual vertex-frame connection has none. Vielbeins and metric fields still carry variables; their lapse and secondary-constraint algebra remains open. The metric/vielbein distinction matters in published [multi-gravity cycle analyses](https://arxiv.org/abs/1410.7774).

## 11377 — a fixed-background axion sector of the wound wall

For the supplied derivative-only compact axions $\theta_x=x+\chi_x$, $\theta_y=y+\chi_y$, hold the Pass 11341 wall geometry fixed. Each Fourier mode with torus integers ($m_x,m_y$) has a nonnegative quadratic form proportional to

$$
h\int r b^2\left(|\partial_\xi\chi|^2+
\frac{m_x^2+m_y^2}{b^2}|\chi|^2\right)d\xi,
\quad 1\le b\le5/4,
$$

so each nonzero torus harmonic has Rayleigh quotient at least $16(m_x^2+m_y^2)/25$. The torus-zero radial block is nonnegative, with its regular zero given by a constant shift. The two constant shifts are exact *nonlinear* moduli of the stated derivative-only action: shifting either compact axion by a constant leaves its derivative and hence the full action unchanged. This resolves those two particular zeros. Nonconstant axion variations couple to metric, Maxwell and wall displacement; positivity of their principal block does not determine the Schur complement or total wall Morse index.

## 11378 — an exact scale-free relation, and the missing scale

At the zero-bare-$\rho$ point of the supplied quantized wall, $q=4/3$, $h=16/45$, $b_w=5/4$, the exact cap function and wall quantities are

$$
F(b)=-\frac{8(b-1)(2b-5)}{45b^2},\quad
F_w=\frac{16}{225},\quad
T^2=\frac{4096}{5625},\quad R_w=\frac{512}{1125}.
$$

Consequently $R_w/T^2=5/8$ in the stated $\kappa^2=1$ convention, or $R_{\rm phys}/(\kappa_{\rm phys}^4 T_{\rm phys}^2)=5/8$ with units restored. The relation survives any overall physical length $L$: $R_{\rm phys}=R_w/L^2$, $T_{\rm phys}=T/(\kappa_{\rm phys}^2L)$. Flux quantization and this exact ratio therefore do not determine $L$ or the gravitational coupling. This extends the prior scale-identifiability audit with a concrete invariant of the wall, not an observed cosmological constant or mass prediction.

## What remains necessary

The five executable sections pass, and each specifies a narrower open problem: a full-field UV-supported CP selector; genuinely independent E6 invariant gradients and momentum-dependent 1PI matching; a vertex-vielbein Hamiltonian constraint calculation; the coupled axion/metric/wall fluctuation operator; and a dynamical physical scale plus vacuum-state selection. Earlier toy selectors and native graph/cycle results are credited above. These calculations do not establish a theory of everything.
