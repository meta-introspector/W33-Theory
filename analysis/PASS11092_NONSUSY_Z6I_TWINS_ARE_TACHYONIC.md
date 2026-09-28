# Pass 11092 — the 87 non-supersymmetric Z6-I W(3,3) twins: consistent, three generations, and all tachyonic

Producer: `analysis/w33_pass11092_nonsusy_z6i_twins_are_tachyonic.py`
orbifolder patches: `analysis/w33_pass11092_orbifolder_n0_patches.py`
Inputs: `data/w33_pass11092_z6i_w33_models.json` (the 87 shifts and Wilson lines),
`data/w33_pass11092_orbifolder_twin_evidence.json` (frozen output of the patched orbifolder)
Certificate: `data/w33_pass11092_nonsusy_z6i_twins_are_tachyonic.json`
Regression: `tests/test_w33_pass11092_nonsusy_z6i_twins_are_tachyonic.py`

Pass 11089 left one door open on the string route. Every one of the 87 Z6-I W(3,3) Standard Models has a
**non-supersymmetric twin**:
* the same shift and Wilson lines;
* the twist (0, ⅙, ⅙, −⅓) replaced by v′ = (0, ⅙, ⅙, ⅔), a 2π rotation attached, i.e. (−1)^F.

Without superpartners there are no dimension-4 or dimension-5 proton-decay operators, and those operators are
what closed the supersymmetric versions. The twins needed a spectrum engine that orbifolder 1.2.1 did not
have. This pass builds one and computes all 87.

## Prior art, stated first

* **The twist is not new.** Font and Hernández, *Non-supersymmetric orbifolds*, Nucl. Phys. B634 (2002) 51,
  [hep-th/0202057]. Their Table 1 lists (⅙, ⅙, ⅔) among the non-supersymmetric Z6 actions. Their §4.3 reports
  tachyon-free models with group SO16 × SO10 × SU2 × U1² in both heterotic strings. They also give the
  level-matching mechanism by which putative tachyons disappear.
* **Non-supersymmetric orbifold software exists.** The *non-SUSY orbifolder* (Escalante-Notario,
  Pérez-Martínez, Ramos-Sánchez, Vaudrevange, [arXiv:2504.20137], 2025) builds orbifolds of the
  **SO(16)×SO(16)** string. There the Witten twist is a separate freely acting Z₂ with gauge shift
  V₀ = (1,0⁷;1,0⁷), and each point-group twist is itself supersymmetric. It builds on two earlier works:
  * Blaszczyk, Groot Nibbelink, Loukas and Ramos-Sánchez, JHEP 10 (2014) 119, [arXiv:1407.6362];
  * the landscape scan of Pérez-Martínez, Ramos-Sánchez and Vaudrevange, PRD 104 (2021) 046026,
    [arXiv:2105.03460], over all 138 Abelian geometries, Z3×Z3 included.
* **The twins here are a different construction.** (−1)^F sits inside the point-group generator, with **no**
  gauge shift of its own. So they are Z6 orbifolds of the E8×E8 string in the Font–Hernández sense, not
  SO(16)×SO(16) models. The five orbifolder patches below were written for this construction. They were not
  compared with the non-SUSY orbifolder's code.

## 1. The engine: five faults in orbifolder at N = 0

Each fault was localised by gdb or a debug print. Each fix keeps the supersymmetric control **byte-identical**:
all 87 parent spectra agree before and after. A clean rebuild from the committed patch module reproduces the
working build byte for byte on 6 sample models.

| # | where | fault | fix |
|---|---|---|---|
| 1 | `CState::FindSUSYMultiplets` | no N = 0 branch | each right-mover is its own multiplet, classified by helicity (−½ / +½ Weyl fermion, ∓1 vector, ∓2 graviton, 0 real scalar pair "H") |
| 2 | `CFixedBrane::FindSUSYMultiplets` | counts over an empty vector (segfault) | N = 0 has one, empty, supercharge combination |
| 3 | `CSector::SortByEigenvalue` | **a right-moving oscillator with number operator N contributes −N mod 1; level matching needs +N** | flip the sign |
| 4 | `CState::CreateRepresentations` | a neutral untwisted fermion is relabelled as a modulus and its unset q_sh is read (segfault) | no moduli at N = 0 |
| 5 | `COrbifold::Create` | anomalous-U(1) generator built only for N = 1, so the anomaly check read an unrotated basis | build it for N = 0 too |

**Fault 3 is the one Pass 11089 left open, and 11089's diagnosis ("missing vacuum phase") was wrong.**
**Prior art (found in Pass 11095):** the non-SUSY orbifolder (arXiv:2504.20137) already carries exactly this
correction in `CSector::SortByEigenvalue` ("wrong trafo sign for R-moving oscillator excitations in original
Orbifolder"). Our patch 3 is an independent rediscovery, and the fix is theirs. Pass 11095 also adds patch 6 (the
order of a twist on spinors) and patch 7 (a mass-level engine for tachyons). The
debug print showed a mismatch of exactly **1 + 2N_R = ⅓ (mod 1)** in every failing state. All of them were the
massless scalar with q_sh = (0, ⅙, ⅙, −⅓) and N_R = ⅙. In supersymmetric models no massless right-mover ever
carries an oscillator, so the sign had never been exercised.

## 2. The twins are consistent three-generation models

All from the patched orbifolder (frozen evidence):
* **Anomaly-free up to Green–Schwarz, 87/87.** SU(3)³ = 0. Every SU(2) has an even number of doublets. Every
  U(1) anomaly is universal, so there is exactly one Green–Schwarz U(1). This is a strong check on the
  engine: a wrong projection phase would break universality.
* **Standard-Model test, 87/87.** orbifolder's own `CAnalyseModel` passes on the left-handed fermions. Net
  generations are **3 q, 3 u, 3 d, 3 l, 3 e in every twin**, as in the parents.
* **Higgs doublet scalars.** (1,2)_{±½}, hidden-singlet: 10 in 55 twins, 30 in 24, and 26 in 8.
* **θ-sector fermions flip helicity, 87/87.** The θ-sector fermions are right-handed in the parent and
  left-handed in the twin. This happens because (−1)^F exchanges the two spinor classes of the right mover in
  odd sectors. The net Standard-Model chirality is unchanged.
* **Gauge group, 87/87.** It equals the parent's: identical vector content and number of U(1)s. Gauge bosons
  have q·n = 0, so (−1)^F does not project them.

## 3. Every twin is tachyonic

**Right movers** (pure Python, exact; `right_movers` in the certificate). In the twin the only tachyonic
right-mover is the scalar q_sh = (0, ⅙, ⅙, −⅓) with N_R = 0. It occurs in the θ and θ⁵ sectors, at
**M²/8 = −⅙**. The sectors k = 0, 2, 3, 4 have none, and the supersymmetric twist has none at all.

A tachyon therefore needs a θ-sector left mover with **½p_sh² + N_L = 7/12**. At each θ fixed point the
centraliser is generated by the constructing element, so every level-matched state survives the projection.

**Left movers.** This is an independent engine (exact enumeration of E8×E8 near V + mW₅, with m = 0, 1, 2),
separate from orbifolder.
* **Control.** Its count of massless θ-sector states (½p_sh² + N_L = ¾) equals orbifolder's count in
  **87/87** models.
* **Known-positive control.** For the Font–Hernández shift (⅙, ⅙, ⅔, 0⁵; 0⁷, 1) it finds gauge group
  SO10 × SU2 × U1² × SO16 (42 + 112 roots) and **no tachyon**: the minimum of ½p_sh² is exactly ¾.

**Result: all 87 twins are tachyonic.**
* Each has **2–20 tachyonic left-moving states**, at 1, 2 or 3 of its 3 fixed points (in 9, 24 and 54 twins).
* Two mechanisms occur: ½p_sh² = 7/12 with no oscillator, and ½p_sh² = 5/12 with one α_{−1/6}.
* An explicit momentum is certified for each model. In model 17 (SM_20260917_55), at the fixed point without
  Wilson line, p = (0⁶, −1, −1; 0⁸) gives ½p_sh² = 7/12.

**Why every one of them.** For this twist the tachyonic level lies exactly ⅙ below the massless level, and ⅙
is the smallest θ-sector left-moving mode. Removing or adding one α_{−1/6} therefore maps twin tachyons to
massless **oscillator-excited** states of the supersymmetric parent, and back. So:

> **A fixed point of the twin is tachyonic iff the supersymmetric parent has a massless oscillator-excited
> left mover there.**

The certificate checks this on all **261** fixed points. It is a restatement of level matching for this twist,
in the Font–Hernández spirit, not a new mechanism. Its content is that tachyon-freedom requires every massless
θ-sector state of the parent to be an oscillator-free ground state with ½p_sh² = ¾. None of the 87 W(3,3)
Standard Models meets this, while the Font–Hernández model does.

## Reading

The last open door of the string route with the W(3,3) twist is closed at the orbifold point. The twins are
anomaly-free, three-generation and proton-safe in the supersymmetric sense. But every one carries charged
tachyons at its θ fixed points, so none is a stable vacuum.

Scope:
* "Tachyonic" is a statement at the orbifold point. What the tachyons condense to, for example a blow-up of
  the θ fixed points, is not computed here.
* The twins keep the parents' shifts and Wilson lines. Other non-supersymmetric completions of the W(3,3)
  twist are not classified: different Wilson lines, or an extra freely acting Z₂ (the SO(16)×SO(16) route).

[hep-th/0202057]: https://arxiv.org/abs/hep-th/0202057
[arXiv:2504.20137]: https://arxiv.org/abs/2504.20137
[arXiv:1407.6362]: https://arxiv.org/abs/1407.6362
[arXiv:2105.03460]: https://arxiv.org/abs/2105.03460
