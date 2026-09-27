# Pass 10979 — the W(3,3) twist on the prime orbifolds: T6/Z3 and Z3×Z3

> **Framing corrected by Pass 11024.** Every number below stands and is independently certified by an exact
> lattice solver (0 differences). But the reading of section E was too strong. On these vacua W is a Z3-graded
> function of one or two composites, W = F(p) with F(ωp) = ωF(p), so the first balance at degree 12 is not a
> hurdle: grad F = 0 has generic solutions. Supersymmetric vacua **exist** in all 28, at ⟨φ⟩ ≈ 0.8–1.1 M_s. The
> numerical search here missed them because it started at the FI scale. The vacua are excluded instead by a
> hidden-representation index: six massless charge-1/3 states in each. See
> `analysis/PASS11024_Z3XZ3_SUSY_VACUA_FRACTIONAL_CHARGE_INDEX.md`.

Producer: `analysis/w33_pass10979_z3xz3_parity_vacua_f_obstruction.py`
Frozen input: `data/w33_pass10979_z3xz3_ledger.json.gz` (exact-rational left-chiral ledger, space-group and
R discrete charges, the scan inputs and statistics)
Certificate: `data/w33_pass10979_z3xz3_parity_vacua_f_obstruction.json`
Regression: `tests/test_w33_pass10979_z3xz3_parity_vacua_f_obstruction.py`

Passes 10974 and 10978 closed every family built on the W(3,3) twist so far, but several of those closures
had to settle which R-selection rule is valid first: the per-plane rule holds only on prime planes
(Bizet et al., arXiv:1301.2322). T6/Z3 and T6/Z3×Z3 on SU(3)³ have **only prime (order-3) planes**, so
every rule orbifolder reports is established. This is the first place the matter-parity question can be
decided with **no rule in doubt**. The twist is the A8 Kac class; the 27 fixed points of T6/Z3 match the
27 points far from a given point of W(3,3).

## A. Scans

* **T6/Z3.** The A8 Kac pair that is modular invariant for Z3 is unique. With it, 64,000 randomised
  Wilson-line draws gave **0 Standard Models**. This corroborates the order-three law already in the paper
  (Holotrade ed1eaca/a6e1c69/a779053): for a cyclic Z_N, u^c and e^c survive together only for even N.
* **T6/Z3×Z3** (V1 = the A8 pair, V2 a modular-invariant second shift). 40,000 draws gave **13 Standard
  Models, 12 of them distinct spectra.** The order-three law is a statement about *cyclic* groups. A second
  independent order-3 twist evades it, consistent with the known Z3×Z3 MSSM-like landscape (e.g.
  Olguín-Trejo, Pérez-Martínez, Ramos-Sánchez, arXiv:1808.06622; Carballo-Pérez, Peinado, Ramos-Sánchez,
  JHEP 12 (2016) 131). What is new here is the W(3,3)-twisted sample and its vacuum analysis. The existence
  of Z3×Z3 Standard Models is not new.

## B. Parity census (torus × space group × every R-combination; all established)

| models | parity exists | closed by Farkas | FI-cancelling rays realizable |
|---|---|---|---|
| 13 | 4 | 1 | **3** |

This is the first family of the programme with parity-preserving FI-cancelling D-flat directions under
established rules.

## C. The minimal vacua are flat-or-massive

Every realizable FI support was run through the exact MSSM-viable matcher (each X̄ paired with an X of
opposite phase, three odd families, even H_u and H_d):
* c1: no viable assignment;
* c3: 24 of 24 viable;
* c4: 12 of 12 viable.

On all 36, the single exotic d-triplet pair is **massless to all orders** (non-negative monomials, ILP),
and W restricted to the support vanishes identically. This is the flat-or-massive dichotomy again. But
outside singlets enter W *linearly* at order 3–4 (n_63 and n_65 in c3, n_62 in c4), so those vacua are
not F-flat and the singlets are forced on.

## D. The vacua the linear terms force

Adding the forced singlets to each support gives 36 distinct extended vacua.
* **8** keep the exotic triplet massless, so they are dead.
* **28** survive, with these properties:
  * the joint parity is still MSSM-viable (one parity element for the whole vacuum);
  * the exotic d is massive;
  * the vacuum is FI-cancelling D-flat with **every** field nonzero (LP certificate).

Their Higgs sector:
* In 22 of the 28 the Higgs pair is massless to all orders, so μ must come from supersymmetry breaking.
* In the other 6 it gets a mass at order 11–12.

## E. F-flatness, exactly

With every vacuum field nonzero,

    F_i = φ_i⁻¹ Σ_k c_k e_{k,i} m_k(φ),

so F = 0 forces Σ_k (c_k m_k) e_k = 0 with every c_k m_k ≠ 0. **The exponent vectors of the monomials
of W on the vacuum must be linearly dependent.** The producer enumerated every W-allowed monomial up to
degree 14 in the vacuum fields (exact integer charge test under every selection rule). Results:

| monomials by degree | vacua | first dependent degree |
|---|---|---|
| {4: 2} or {4: 1} | 12 | none through 14 |
| {3: 2, 12: 5} | 6 | 12 |
| {3: 1, 4: 1, 12: 1, 13: 1, 14: 1} | 6 | 12 |
| {3: 1, 12: 1} | 4 | 12 |

**No supersymmetric vacuum exists from the superpotential through degree 11, for any nonzero couplings, in
any of the 28.** In 12 of them that holds through degree 14. The first possible balance sets cubic or
quartic terms against degree-12 terms, |⟨φ⟩/M_s|⁹ ~ c_low/c_12, so the VEV is fixed by a *coupling ratio*
rather than by the FI term. If ⟨φ⟩ is a few tenths of M_s, as the FI term sets it elsewhere in the
programme, the degree-12 coupling must exceed the cubic by 10³–10⁶. A generic numerical F = D = 0 search
(random O(1) couplings, 6 seeds per vacuum, in the scratch record) found no root with every field
nonzero, in agreement.

## Verdict and the open question it names

T6/Z3×Z3 is the first W(3,3) family where matter parity survives the Fayet–Iliopoulos term under
established rules. Its vacua cannot be supersymmetric through degree 11. **Whether heterotic couplings can
supply the inverted hierarchy c_12 ≫ c_3 is a worldsheet-instanton question.** Cubic couplings between
twisted fields at different fixed points are exponentially suppressed, e^{−A R²}, so a large compact
radius *could* invert it. This is the one precisely named open door of the string route: compute the
instanton suppression of the cubic terms in the c3/c4 vacua.

Scope: the census covers the 13 Standard Models found, not all of Z3×Z3. The obstruction is exact
through the stated degree and says nothing beyond degree 14. The numerical search is corroboration only.
