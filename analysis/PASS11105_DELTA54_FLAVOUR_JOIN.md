# Pass 11105 — the flavour group of the W(3,3) string models is Δ(54) = H27 : ⟨−I⟩, the other track's physical qutrit Heisenberg group

Producer: `analysis/w33_pass11105_delta54_flavour_join.py` (parts A, B, D need no data; part C is frozen from the dumps)
Certificate: `data/w33_pass11105_delta54_flavour_join.json`
Regression: `tests/test_w33_pass11105_delta54_flavour_join.py`

This pass joins the two tracks. The other track built the qutrit Heisenberg group H27 = 3^(1+2), its Clifford
normaliser H27 : SL(2,3) of order 648 (= the W33 point stabiliser), and the two-qutrit Pauli group 3^(1+4) of order 243
as finite subgroups of E8. Those constructions live in `analysis/2026-09-21_physical_external_a2_h27.md`,
`..._physical_a2_clifford648_w33_bridge.md` and `..._e8_trinification_two_qutrit_pauli243.md`, which all state "no
heterotic vacuum inferred". This track built explicit heterotic vacua on T⁶/Z3, whose 27 fixed points form AG(3,3).
Here is how the two meet.

## A. The symmetry group of the Z3 selection rule is Δ(54) (rule level, no spectra)

On one torus, a twisted field θ^l carries a fixed-point class c ∈ F3. Two rules apply: the space-group rule Σc ≡ 0 and
the point-group rule Σl ≡ 0 (mod 3).

**Search.** I searched all permutation-phase maps of the classes in both twisted sectors: 6·6·27·27 candidates. A
candidate is kept if every allowed monomial up to degree 6 is invariant (311 monomials were checked) and every
forbidden one stays forbidden.

**Result:**

| | |
|---|---|
| order of the symmetry group | **54**, equal to Δ(54) |
| generators | X: c → c + l (translation), Z: phase ω^c (the space-group rule), R: c → −c |
| Heisenberg relation | ZX = P·XZ with **P = ω^l, the point-group twist** |
| ⟨X, Z⟩ | H27 = 3^(1+2), with centre = the orbifold twist |
| conjugation image on H27/Z = F3² | {+I, −I} only |
| index in Clifford-648 | **12** |
| quadratic phase ω^(c²) (the Clifford gate X → XZ) | **not** a symmetry: the collinear triple (0,1,2) has Σc² ≡ 2 |

So the physical vacuum realises the other track's H27 together with only the central −I of its SL(2,3). The rest,
A4 = PSL(2,3), would need quadratic phases, and the space-group rule forbids them. Δ(54) as the flavour group of a
Wilson-line-free Z3 torus is known (Kobayashi, Nilles, Plöger, Raby, Ratz, hep-ph/0611020). What is new is the
identification with the corpus objects.

## B. Geometry

A θ^l field with class label c sits at the geometric fixed point g = l·c; the θ² labels are inverted. In these terms the
rule reads Σ l_i g_i ≡ 0. Consequences:
* θθθ couplings are exactly the collinear triples of AG(3,3). There are 117 lines, and 12 in the two-torus Hesse
  plane AG(2,3).
* A θ × θ² × untwisted coupling requires the **same** point.

## C. One qutrit, never two (104 spectra of Pass 11095)

In 104 of 104 models, the spectrum is invariant under X_t and R_t for **exactly one torus**, the family-triplication
torus: torus 0 in 13 models, 1 in 38, 2 in 53. The Wilson lines of the other tori are random in the scan, so this is
forced by three generations, not built in.

The geometric flavour group is therefore one Δ(54). The two-qutrit Pauli group 3^(1+4) of order 243 would need two
Wilson-line-free tori, and it is realised in **0 / 104** models. Its projective commutation graph is recomputed here as
SRG(40,12,2,4) = W(3,3). **W(3,3) is therefore not the flavour geometry of these vacua; one torus's qutrit is.**

## D. The Higgs and Schur's lemma

A Higgs localized at a point p0 leaves the stabiliser ⟨Z·P^(−p0), R_p0, P⟩, of order 18 = Z3 × S3. On the other two
generations this stabiliser acts irreducibly: its commutant has dimension 1. While it is unbroken, the charm–up block of
Y†Y is therefore a multiple of the identity, so **m_c = m_u exactly**.
* Pass 11102 saw the order-level shadow of this.
* Pass 11103 shows that no alignment of the scalar vacuum lifts it parametrically in the up sector.

## Reading

* The other track's "physical H27" is realised, in explicit string vacua, as the family symmetry of the three
  generations. Its centre is the orbifold twist itself.
* The larger objects are not realised by the geometry: the Clifford normaliser, and the 243-element two-qutrit Pauli
  group whose geometry is W(3,3).
* The −I is the reflection that keeps charm and up degenerate.
