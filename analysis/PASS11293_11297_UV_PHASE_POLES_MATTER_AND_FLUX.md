# Passes11293–11297: execute the five remaining physical-frontier targets

Reservation4b860f9ca follows the published11285–11292 packet067262c9b. The parallel inventory-only commitf042b70a5 was reviewed and integrated first. All five producers execute scoped PASS. Successful construction and successful obstruction tests are distinguished below. No observed parameter or completed TOE claim follows from this packet. These Higgs, flavor, condensate and gravity constructions are separately declared models, not one fully matched UV action.

## 11293 — explicit UV matching exposes the running obstruction

The desired operator is d_ABC ψ^{Ai}ψ^{Bj}v^C F†_ij/M, with ψ in(27,3), v in(27,1), and F† in(1,bar6). Introduce vectorlike Weyls X in(27,bar3) and Xbar in(bar27,3), with

M X^{A}_i Xbar_A^i + y1 d_ABC ψ^{Ai}v^B X^C_i + y2 Xbar_A^i ψ^{Aj} F†_ij.

Eliminating the pair gives the desired dimension-five interaction with coefficient −y1 y2/M (the Weyl mass matrix carries the corresponding identical-field factor). This is checked using the actual signed E6 cubic, an81-dimensional light sector and162-dimensional heavy sector: the heavy-block Schur complement equals −2(D_v⊗F†)/M. A heavy scalar H in(27,bar6) also generates the operator through its Yukawa and cubic H†vF† couplings.

However, T_E6(27)=3 fixes the added one-loop costs: six complex27 copies cost6, while two sets of three Weyl27 copies cost12. The existing composite inventory has bE6=34−6r and bSO(2r)=7r−34. Auxiliary AF needs integer r≥5; at r=5 the scalar channel gives bE6=−2 and the fermion channel gives−8. Neither single-propagator completion preserves joint AF anywhere in this inventory.

### An economical alignment escape

The obstruction motivated a different inventory: replace the third complex(27,10) frame W by one real symmetric traceless SO(10) tensor S in54. Keep K in45 and two frames U,V, and impose

[K,S]=0, S²+S−6I=0, AU=US,

alongside the old base, projection, mixed-factor and Gram constraints. Every constraint is still at most quadratic, so the scalar potential is quartic. At the actual reference, S=U†AU+its complex conjugate is real symmetric traceless, with eigenvalues2 (six times) and−3 (four times).

For the local kernel, first remove the20 K-orbit variations. Reality and commutation with K make S uniquely recoverable from U†AU; its differential is likewise unique. The intertwiner and S polynomial imply(A²+A−6I)U=0. Multiplication by V† maps the remaining kernel into the old exact11275 selector kernel. Subtract its66 E6 directions; the residual25 frame directions are U5 gauge transformations. The actual integer compact-orbit tangent has rank111 modulo101, providing the matching lower bound. Thus1365 alignment fields have111 gauge and1254 positive normal directions. This is analytic elimination using the prior rank theorem, not a floating full Hessian or a global uniqueness claim.

Removing W saves ten complex E6 fundamentals and27 complex SO10 vectors; adding real54 costs only T54/6=12/6=2 auxiliary beta units. The composite inventory becomes(bE6,bSO10)=(14,8). The explicit scalar mediator leaves(8,8), and the vectorlike fermions leave(2,8): **both now preserve positive E6/SO10 one-loop gauge coefficients**.

For the scalar channel, the added nonnegative renormalizable potential ||M H+α vF†||² makes the heavy-scalar zero locus a graph H=−αvF†/M over the alignment vacuum. Its stabilizing quartic is explicit, rather than an omitted negative threshold potential. The family-sextet dynamics remains a separate input. Positive gauge coefficients do not prove complete asymptotic freedom of all Yukawa/quartic couplings, and family SU3 still fails AF. The new inventory changes the architecture; it does not retract the old fixed-inventory obstruction. The earlier SO10 Higgs literature remains prior art: [SO10 symmetry breaking](https://doi.org/10.1103/PhysRevD.24.1005).

The representation argument covers the two single-propagator tree topologies with precisely the named external fields. It is not a no-go theorem for strong/composite, radiative, asymptotically safe, altered-inventory or more elaborate UV models. The matching object now exists; the small running budget does not survive it. Prior owners:11271,11285 and11291. Group conventions and indices use the standard [Slansky tables](https://doi.org/10.1016/0370-1573(81)90092-2).

## 11294 — hierarchical CP and the missing phase constraints

Replace the Fourier target by the declared one-parameter standard mixing matrix with s12=ε, s23=ε², s13=ε³ and δ=π/2. Four mixed-moment targets then fix its modulus matrix by the same distinct-spectrum Vandermonde argument as11286. The exact identity is

J²=ε^12(1−ε²)(1−ε^4)(1−ε^6)².

To remove the four relative sextet phases invisible to FF†, use Z=adj(Fu)Fd and its traces z1=Tr Z, z2=Tr Z². Under SU3 congruence, Z transforms by similarity, so both traces are invariant. Add positive terms |z_k²−2Re(t_k)z_k+|t_k|²|². Their coefficients are real and their roots are CP conjugates, so the potential remains CP even. The reference roots t_k are calculated from the declared vacuum.

At ε=0.18,0.22,0.30,0.40 the full24-field constraint Jacobian has rank16, leaving the eight gauge directions; the phase-invariant Jacobian on the old four flat directions is nonsingular. SU3 covariance, CP-conjugate zeros and phase displacements are independently controlled. At ε=0.22, J≈0.00011046. These are numerical local-isolation checks, not an interval theorem for every ε or a global classification.

The spectra, ε powers and invariant-root coefficients are inputs. This is an explicit phase-stabilized hierarchical toy architecture, not an observed CKM prediction. The highest polynomial degree is24. The adjugate is polynomial; numerical inverse evaluation is used only at nonsingular matrices. Prior owner:11286 and BT891. Classical flavor-potential context: [Alonso et al.](https://arxiv.org/abs/1103.2915).

## 11295 — conditional complex pole blocks with explicit matching freedom

Reuse the full22-field scalar/Weyl channel matrices from11287. Let P project onto the eight Goldstone directions. At zero momentum, take only the finite Goldstone-column limits Σ0P: massless scalar channels annihilate these columns, so the divergent normal-normal zero-momentum logarithm is not evaluated. A symmetric constant

C=−(Σ0P+PΣ0−PΣ0P)

enforces(Σ0+C)P=0. Arbitrary real changes of external basis preserve this construction. At the spacelike subtraction point s0=−0.0001, the normal-normal analytic mass and kinetic matching coefficients are set to zero **as a scheme input**.

For inverse propagator s−Mtree²+Σ(s), diagonalize Σ(Mtree²) inside each degenerate tree block. This gives fourteen leading complex massive squared poles plus eight massless Goldstone poles. Off-block mixing and momentum iteration enter at higher perturbative order. In the declared scheme, the lowest ten-mode squared-pole block starts around0.000330611; the radial squared pole is approximately0.009382185−i0.000033384, with leading width/mass0.003606259. These scheme-dependent shifts are not transferred from the old fixed-slice subtraction prescription.

The E6 vector mass matrix at trace-normalized gE=0.5 has70 massive vectors. Its smallest pair threshold squared is0.01733004, above the largest scalar tree value0.00925714. The eight unbroken generators act trivially on the quotient slice. Thus these massive-vector pair cuts are closed at these inputs. Heavy-vector analytic terms, gauginos and UV matching are still not computed.

Most importantly, local two-point Ward matching is weaker than an invariant tadpole/counterterm functional. Changing an allowed normal analytic coefficient shifts the massive poles at this same order. The outputs are conditional matched-EFT poles, not finished UV physical predictions. Prior owners:11282,11287,11290; primary checks: [Braathen–Goodsell](https://arxiv.org/abs/1609.06977), [general scalar self-energies and tadpoles](https://arxiv.org/abs/1910.02094).

## 11296 — matter coupling that preserves the metric redundancy

Declare the continuum EH action plus22 quotient scalar fields with positive regular target metric G_ab(φ). Every metric-source-field occurrence is through the existing g(source11). The canonical matter densities are

H_m=p_a G^ab p_b/(2√h)+√h h^ij G_ab ∂iφ^a∂jφ^b/2+√h V,
D_i=p_a∂iφ^a.

The lapse and shift multiply the total constraints and acquire no time derivatives. Since matter introduces no velocities for the11 metric-source fields, the four lapse/shift nulls and chart fiber survive. The exact chart fiber is checked again. The continuum matter bracket is

{H_m[N],H_m[M]}=D_m[h^ij(N∂jM−M∂jN)].

Potential and target-metric derivative terms cancel. An exact symbolic local identity and periodic Fourier controls at64 and128 points check the matter bracket. Total ADM closure follows from the declared covariant EH-plus-matter action; this is not an exact finite-lattice claim. Counting33 configuration variables and nine metric first-class constraints gives24 physical configurations: two gravitons plus22 matter scalars. Positive target metric holds only on the regular quotient chart.

The real base, action and Planck scale remain supplied. No separate source-field sigma kinetic is added:11288 already explained why that could spoil the constraints. Quantum anomalies, a global chart and a convergent W33 discrete gravity limit remain open. Prior owners:11283,11288 and11287; classical canonical gravity: [DeWitt](https://journals.aps.org/pr/abstract/10.1103/PhysRev.160.1113).

## 11297 — membrane dynamics selects a sector but does not solve vacuum energy

Declare a compact three-form, intensive flux f=*F4, charge q, offset f0 and membrane tension T. Homogeneous sectors have f_n=f0+nq and energy ρ_n=Λbare+f_n²/2. Use the known flat-space thin-wall downward transition with Δρ>0:

R=3T/Δρ, B=27π²T^4/(2Δρ³), Γ=A exp(−B).

Both stationary radius and bounce action are checked symbolically. A finite Markov generator implements these transitions and conserves probability. For q=1,f0=0.23,T=0.2,A=1,Λbare=−0.01, evolution from n=6 reaches the unique absorbing sector n=0, whose energy is0.01645. Varying Λbare leaves the rates and terminal sector unchanged but shifts the final energy directly. This is relaxation toward nearest lattice flux, not radiative vacuum-energy cancellation.

The crucial connection is the distinction between intensive f and extensive Q=∫F. The earlier uniform linear constraint Vol=Q/μ4 implies f=Q/Vol=μ4 and Maxwell density μ4²/(2e²), independent of Q. Moving Q along that branch is not the Brown–Teitelboim intensive ladder. A naive combination of the two homogeneous models therefore cannot be used to claim vacuum-energy selection.

This does not rule out inhomogeneous bubbles under a global constraint, nonlinear sigma, multiple fluxes or moduli-dependent membranes. Those require their actual junction and global equations. Gravitational/CDL corrections, upward transitions, cosmological measure, prefactors and a membrane spectrum from W33 remain missing. Prior owners:11284 and11289; imported mechanism: [Hirano](https://arxiv.org/abs/1804.09985), [Bandos et al.](https://arxiv.org/abs/2306.09412), [local sequestering](https://arxiv.org/html/1505.01492v2).

## Validation and reproduction

Run the five matching `analysis/w33_pass11293_*.py` through `w33_pass11297_*.py` producers. They write scoped JSON certificates. Run `python3 -m pytest -q tests/test_w33_pass11293_11297_frontiers.py` for seven independent controls, including exact symbolic hierarchy, signed-tensor Schur matching, phase covariance, arbitrary-basis Ward matching, Fourier constraint closure and flux dynamics. The existing11285–11292 packet is prior work; standard tree matching, Jarlskog identities, ADM algebra and membrane bounce formulas are applications, not claimed discoveries.
