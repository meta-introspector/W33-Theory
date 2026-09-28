# Pass 11096 — the non-supersymmetric twins of Z12-I, Z2×Z6-I, Z3×Z6 and Z6×Z6

Producer: `analysis/w33_pass11096_twins_of_the_other_w33_families.py`
Frozen evidence: `data/w33_pass11096_family_twin_evidence.json`, `data/w33_pass11096_z2xz6i_twin_chiral_content.json`
Certificate: `data/w33_pass11096_twins_of_the_other_w33_families.json`
Regression: `tests/test_w33_pass11096_twins_of_the_other_w33_families.py`

Passes 11092–11094 settled the Z6-I twins. Pass 11089 listed further W(3,3) families whose Standard Models
admit non-supersymmetric twins with the same shifts and Wilson lines. This pass computes all of them.

## A. A correction to Pass 11089

For a Z_M × Z_N twin, modular invariance with unchanged shifts needs **gcd(M,N)(v₁′·v₂′ − v₁·v₂) to be even**.
Pass 11089 tested only that it is an integer. Its listed twists for **Z2×Z6-I, Z3×Z6 and Z6×Z6 give −1** and
fail; orbifolder rejects them ("2 (V_N x V_M − v_N x v_M) = 1 = 0 mod 2 failed"). The existence claims survive,
because correct twins exist: 26, 26 and 52 single-generator insertions of (−1)^F respectively (34 for
Z2×Z6-II, 62 for Z12-I), and still none for Z3×Z3. The Pass 11089 producer now tests evenness and has been
regenerated. Representatives used here:

| family | twin twist (the other generator unchanged) |
|---|---|
| Z12-I | (1/12, −5/12, −2/3) |
| Z2×Z6-I, Z3×Z6 | v₂ = (0, 7/6, −1/6) |
| Z6×Z6 | v₁ = (7/6, 0, −1/6) |

None changes the order of a twist on fermions.

## B. Where tachyons can live (exact)

The tachyonic right-mover levels, computed sector by sector:
* **Z12-I:** −1/12 in θ, θ⁵, θ⁷, θ¹¹; −¼ in θ³, θ⁹.
* **Z2×Z6-I:** −⅓ and −⅙ in (0,1) and (0,5).
* **Z3×Z6:** −⅓ and −⅙ in (0,1) and (0,5); −⅙ in (1,5) and (2,1).
* **Z6×Z6:** −⅓ and −⅙ in (1,0), (1,5), (5,0) and (5,1); −⅙ in (1,1), (1,4), (5,2) and (5,5).

## C. Results (patched orbifolder; tachyons from the mass-level engine of Pass 11094 at every level of B)

| family | SM twins | N = 0 | anomaly-free (one GS U(1)) | tachyonic | tachyon-free |
|---|---|---|---|---|---|
| Z12-I | 289 | 289 | 289 | **289** (12–96 states) | 0 |
| Z2×Z6-I | 29 | 29 | 29 | 21 | **8** |
| Z3×Z6 | 5 | 5 | 5 | **5** | 0 |
| Z6×Z6 | 10 | 10 | 10 | **10** | 0 |

The supersymmetric parents of all four families have no twisted state at the tachyonic levels (control).

## D. The 8 tachyon-free Z2×Z6-I twins are not Standard Models

Applying each supersymmetric parent's hypercharge (Y·Y = 5/6) to its twin's fermions: none of the 29 Z2×Z6-I
twins, and in particular none of the 8 tachyon-free ones, has the Standard-Model chiral content. Examples:
* twin 8 has 4 net quark doublets, 4 lepton-type doublets and a net chiral (3,1)_{1/6};
* several twins keep net chiral fermions (1,1)_{±½} of half-integer charge.

orbifolder's Standard-Model test rejects all 29 twins as well. (That test is supersymmetric in spirit, so it is
not decisive on its own; the hypercharge count is.) Unlike Z6-I, where all 87 twins kept three generations, the
(−1)^F twin of Z2×Z6-I changes the chiral content.

## Reading

Across every W(3,3) family with a non-supersymmetric twin (Z6-I, Z12-I, Z2×Z6-I, Z3×Z6, Z6×Z6; Z2×Z6-II has no
Standard Models, and Z3 and Z3×Z3 have no twins), **no twin is both tachyon-free and a three-generation
Standard Model**. The twin route (same shifts, (−1)^F in the twist) is closed for the W(3,3) twist. The
separate SO(16)×SO(16) route is not; see Pass 11095.
