# 2026-10-01 — E8/G26 mirror-separation theorem

Producer: `analysis/w33_20261001_g26_e8_mirror_separation.py`
Certificate: `data/w33_20261001_g26_e8_mirror_separation.json`

## Result

Pass 11261 and Pass 11269 put two exact relative invariants on the same
21-hyperplane arrangement in the new semisimple Cartan three-plane.

Let

- (A) be the product of the nine qutrit SIC mirror forms;
- (B) be the product of the twelve qutrit stabilizer/MUB mirror forms;
- (P) be the restricted signed-(E_8) complementary Pfaffian;
- (J) be the standard (G_{26}) reflection Jacobian, pulled back to the Cartan plane.

Their mirror valuations are

[
operatorname{div} P=(3,1),qquad
operatorname{div} J=(1,2)
]

in the ordered basis (SIC, stabilizer).
The valuation matrix therefore has determinant five:
[
detegin{pmatrix}3&1\\1&2end{pmatrix}=5.
]
More strongly,
[
2,operatorname{div}P-operatorname{div}J=(5,0),qquad
3,operatorname{div}J-operatorname{div}P=(0,5).
]

So the two mixed relative invariants separate the two qutrit resource strata
up to fifth power:
[
P^2/Jpropto A^5,qquad J^3/Ppropto B^5.
]

With the frozen normalisations, (delta=omega-omega^2), (delta^2=-3),
the cleared polynomial identities are
[
5,3^{20}P^2+2^{30}A^5J=0,
]
[
32J^3-3^{22}5^3delta,B^5P=0.
]
The degrees close independently: (78=45+33) and (99=60+39).
## Differential form

Away from the mirrors the same statement is
[
2,dlog P-dlog J=5,dlog A,
qquad
3,dlog J-dlog P=5,dlog B.
]

Thus the signed-(E_8) cubic Pfaffian and the (G_{26}) discriminant provide
two exact logarithmic probes whose integer combinations isolate the SIC walls
and the stabilizer/MUB walls separately. This is useful for any future
Cartan-plane potential because the two reflection strata can now be coupled
independently without inventing a basis-dependent classifier.

## Literature boundary

Shephard–Todd already identifies (G_{26}) as the order-1296 Hessian
reflection group with invariant degrees (6,12,18), and the Hessian case has
two reflection/homology classes. The (G_{26}) Jacobian factorisation is
therefore classical in character. The repository increment is the exact
comparison with the independently constructed signed-(E_8) Pfaffian on the
explicit Cartan plane, and the resulting index-five divisor separation.

## Firewall

The integer five is the determinant of a two-by-two valuation matrix. It is
not a family count, particle multiplicity, coupling, physical (mathbb Z_5)
symmetry, or prediction. No kinetic metric, Hamiltonian or vacuum selection is
derived here.
