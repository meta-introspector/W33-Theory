# Pass 11218 — the cubic rank 78+3 split is the Vinberg Cartan split

The signed $E_6$ cubic had repeatedly produced an $81\times81$ skew Jacobian of rank $78$.  Pass 11218 closes the gap between that computation and the classical order-three $E_8$ geometry: the residual three modes are an explicit Cartan subspace, not random numerical zero modes.

## Exact repository-gauge witness

In the canonical coordinates $n=3i+a$, the nine-coordinate background stored in the certificate has

\[
\operatorname{rank}D_v=78,\qquad \dim\ker D_v=3.
\]

Exact rational nullspace computation gives three integral kernel vectors.  All three pairwise brackets vanish, and the background is itself the second basis vector.  Therefore

\[
\ker D_v=\mathfrak c,\qquad [\mathfrak c,\mathfrak c]=0,
\qquad \dim\mathfrak c=3.
\]

The rank witness is compact.  Delete flat coordinates $1,2,16$, corresponding to $(i,a)=(0,1),(0,2),(5,1)$.  The complementary $78\times78$ principal minor has

\[
\det=2^{76}3^{24}5^2 7^2
 =\left(2^{38}3^{12}\cdot5\cdot7\right)^2\ne0.
\]

The square is required because the principal matrix is skew-symmetric; its absolute Pfaffian is $2^{38}3^{12}\cdot35$.

## Why this proves the generic theorem

The already-certified root grading is the inner order-three grading

\[
\mathfrak e_8=(\mathfrak e_6\oplus\mathfrak{sl}_3)
 \oplus(\mathbf{27}\otimes\mathbf3)
 \oplus(\overline{\mathbf{27}}\otimes\overline{\mathbf3}).
\]

Reeder, Levy, Yu, and Gross classify this grading in Table 21, $E_8$ row 3b, of *Gradings of positive rank on simple Lie algebras*.  Its theta rank is three, its little Weyl group is Shephard–Todd $G_{26}$, and its invariant degrees are

\[
6,\quad12,\quad18.
\]

The rank-three classification supplies the corank-at-least-three bound for the grade-one centralizer.  Our nonzero $78\times78$ minor attains that bound in the actual signed repository tensor.  Since maximal rank is a Zariski-open condition, $78$ is the generic Jacobian rank and (3) is the generic kernel dimension.

This also identifies the correct algebraic variables for vacuum selection: the generic quotient has three independent invariant coordinates of degrees $6,12,18$.  A physical potential should therefore be written on those three coordinates, then pulled back to the 81-dimensional matter space for Hessian and mixing calculations.

## The $G_{26}$ distinction

This is a genuine role for $G_{26}$, but it does not undo `data/w33_extended_clifford_g26_no_go.json`.  The action here is the little Weyl action on the three-dimensional Cartan slice.  The earlier no-go concerns an abstract identification of the full retained-phase extended qutrit Clifford group with $G_{26}$; their centers and projective quotients differ.  These are different objects acting on different spaces.

The corpus search also checked the earlier Shephard–Todd tower in `PASS1020_E8_TRANSITIVE_51840.md`, `analysis/w33_pass1020_e8_transitive_51840.g`, `analysis/w33_pass1039b_gaussian_base.g`, `analysis/w33_eisenstein_forcing.py`, `analysis/w33_eisenstein_grand_synthesis.py`, and `analysis/w33_pass1047_eisenstein_parabolic_ladder.g`, plus the Clifford and rank-four polytope treatments in `analysis/w33_pass8909_8924_the_centraliser_is_the_clifford_group.py` and `exploration/WITTING_W33_S12_SYNTHESIS.py`.  Those results concern the rank-four G32 action, its rank-three G25 parabolic, or Clifford centralizers.  None constructs the G26 Cartan slice in this signed 81-coordinate bracket.

## Boundary

The theorem supplies the generic vacuum quotient and an explicit Cartan slice.  It does not select a vacuum, prove stability or positivity, determine a potential's coefficients, or derive a mass or mixing observable.  The next constructive calculation is the most general low-degree $G_{26}$-invariant potential in the degree-$6,12,18$ generators and its exact Hessian on the certified slice.

Evidence:

- `analysis/w33_pass11218_e8_cubic_cartan_kernel.py`
- `data/w33_pass11218_e8_cubic_cartan_kernel.json`
- `tests/test_w33_pass11218_e8_cubic_cartan_kernel.py`
- [Reeder–Levy–Yu–Gross, Table 21](https://abel.math.harvard.edu/~gross/preprints/PosRank.pdf)
