# Passes 11194–11198 — five exact TOE frontier audits

Producer: `analysis/w33_pass11194_11198_five_toe_frontiers.py`  
Certificate: `data/w33_pass11194_11198_five_toe_frontiers.json`  
Regression: `tests/test_w33_pass11194_11198_five_toe_frontiers.py`

These passes execute the five fronts reserved in `af90f5404`. They were recomputed after reading the parallel
Passes 11199–11211. In particular, Pass 11207 now gives the exact all-register arrow law
\(A(S)=n-c(S)\), Pass 11209 identifies the phase-sensitive Maslov/Bargmann chirality, and Pass 11206 proves that
the Clifford layer has no T violation. Pass 11213, which arrived during this audit, then finds the minimal
one-qutrit escape from every substrate time reversal: a cubic phase together with a nonzero cyclic shift. That is a
gate-level CP/T mechanism.  The earlier finite charge-conjugation scaffold in
`analysis/2026-09-23_execute_all5_plus3_physics_frontier.md` already marked spatial parity, measured mixing phases,
EDMs, and spontaneous CP breaking as open; deriving any of those quantities or masses still requires the signed
27-coordinate dynamics.

## 11194 — reversible E₆-sector measures do not produce Yukawa weights

Relative to an E₆ frame, the 45 cubic monomials split as \(5+40\). The perfect-gate adjacency quotient is

\[
Q=\begin{pmatrix}4&8\\1&11\end{pmatrix},\qquad P=Q/12.
\]

Let every \(1\cdot10\cdot10\) monomial have weight \(a\) and every \(16\cdot16\cdot10\) monomial weight \(b\).
There are 40 cross-sector edges. Detailed balance on any one of them gives \(a/12=b/12\), hence \(a=b\).
The unique reversible measure on the connected symmetric 45-state walk is uniform. Its sector masses are
\((1/9,8/9)\) only because the sectors contain \((5,40)\) vertices. The ratio is a census, not a Yukawa hierarchy.

This closes the most conservative symmetry-compatible probability proposal. Nonreversible dynamics, a chosen
vacuum, or amplitudes carried by the cubic tensor remain open.

## 11195 — the optimal perfect-gate compiler is phase complete

The existing phase ABI gives an explicit 9×9 unitary to each of the 80 W33 transvections and tracks displacement
and the Weyl cocycle. A fresh breadth-first enumeration gives shortest transvection depths

\[
0^1,\;1^{80},\;2^{1980},\;3^{16650},\;4^{33128},\;5^1.
\]

After the exact basis conjugacy, the fixed perfect gate \(p\) and the middle local gate \(h_*\) both have canonical
depth four. The stored SUM target has canonical depth two. Expanding its optimal perfect-gate normal form
\(Lph_*pR\) gives 19 transvections before cancellation. The two 9×9 unitaries agree with scalar ratio \(1\) and
residual below \(1.5\times10^{-15}\). The compiler join therefore preserves the global phase convention as well as
the symplectic matrix.

This completes the algebraic Clifford compiler. Approximate universality still uses the separately typed Pass 10944
\(|0\rangle\)-controlled-X cubic resource with the qutrit Fourier gate. Pass 11213 independently shows that a cubic
phase plus a cyclic level shift generically breaks the substrate's anti-unitary Clifford time reversals. Physical
magic-state injection, the two-qutrit extension, and hardware calibration remain open.

## 11196 — a linear signed-incidence Yukawa ansatz is sign blind

Let \(B\) be the 45×27 tritangent/line incidence matrix and let \(D\) contain the canonical cubic signs. The signed
linear transport is

\[
B_s=DB,\qquad D^2=I,
\qquad B_s^TB_s=B^TB.
\]

Consequently its rank remains 21 and the centered Gram spectrum remains \(6^{20}+0^7\). Every singular value and
every right singular subspace is independent of the cubic signs. A mass or mixing ansatz using only this linear Gram
data cannot use the very sign information that distinguishes the E₆ cubic.

This is a route-selection result. The nonlinear cubic Jacobian is sign sensitive and already has exact tested ranks
20, 54, 54 and 78 on four certified backgrounds. A viable Yukawa construction should therefore use the cubic
gradient, Hessian, Jacobian, or products coupling different triads, followed by a vacuum-selection principle.

The current canonical artifact uses 23 negative and 22 positive triads. Reversing the overall cubic sign gives the
older prose convention of 23 positive and 22 negative; the invariant content is the 23-versus-22 imbalance.

## 11197 — the W(3,q) arithmetic tower has a collapsed heat spectrum

For \(W(3,q)\), the nonzero combinatorial Laplacian eigenvalues are

\[
q^2+1\quad\text{and}\quad(q+1)^2
\]

with multiplicities \(q(q+1)^2/2\) and \(q(q^2+1)/2\). After degree normalization, both eigenvalues tend to one.
The normalized spectral measure converges to \(\delta_1\), and for every fixed \(t>0\),

\[
\frac1{|V|}\operatorname{Tr}e^{-tL/k}\longrightarrow e^{-t}.
\]

The corresponding effective spectral dimension tends to \(2t\), rather than developing a scale-independent
dimension plateau. Every graph in the tower also has diameter two. This sharpens the repository's earlier
three-eigenvalue/no-Weyl-law result: arithmetic field extension creates more vertices but no manifold-like
low-frequency ladder. A metric spacetime limit still needs additional structure.

## 11198 — the Hecke trace zero is not vacuum cancellation

For the normalized perfect-gate walk,

\[
\operatorname{spec}(P)=\{1^1,(1/4)^{20},(-1/4)^{24}\},\qquad
\operatorname{Tr}(P^n)=1+\frac{20+24(-1)^n}{4^n}.
\]

Thus \(\operatorname{Tr}P=0\), but every tested higher moment \(n=2,\ldots,8\) is nonzero. More strongly, if a
\(\mathbb Z_2\) grading \(\Gamma\) commutes with \(P\), cancellation for every power would force the graded traces
on all three eigenspaces to vanish. That is impossible on the one-dimensional Perron eigenspace.

The isolated first-moment zero is unrelated to the Hodge indices \(-40\) and \(-80\), and it does not
cancel Pass 11099's independently computed positive ten-dimensional string vacuum energy. A finite Markov trace is
not a cosmological constant without a physical grading and Hamiltonian.

## Net effect on the TOE frontier

The packet removes three tempting shortcuts: orbit counts do not generate Yukawa amplitudes, cubic signs cannot act
through a linear incidence Gram matrix, and the Hecke first-moment cancellation is not supersymmetry. It also closes
the algebraic phase join for the universal-computation stack and strengthens the exact obstruction to obtaining
metric spacetime from the raw \(W(3,q)\) tower. Pass 11213 supplies the minimal gate-level T-violating mechanism. The
next positive flavor construction should be a nonlinear signed-cubic vacuum operator whose eigenspaces are tested
against the 27 matter coordinates and whose rephasing-invariant three-generation phase reduces to that cubic-plus-
mixing mechanism.
