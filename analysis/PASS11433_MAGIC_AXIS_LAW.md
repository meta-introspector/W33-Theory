# Pass 11433 — the magic-axis law: a proved grading structure, and an exhaustive four-qutrit check of the fixed-axis case

Producer: `analysis/w33_pass11433_magic_axis_law.py`
GAP: `analysis/gap/w33_pass11433_stab_z1_classes.g` (generators in `analysis/gap/w33_pass11433_gens4.g`) →
`data/w33_pass11433_stab_z1_classes_n4.txt`
Certificate: `data/w33_pass11433_magic_axis_law.json`
Regression: `tests/test_w33_pass11433_11437.py`

**The law (Pass 11420, F′).** If M²z₁ = z₁, then v = z₁ + Mz₁ is fixed by M. Every frame a with ω(a, v) ≠ 0 violates.

## Grading lemma (proved)

* **Statement.** Let M²z₁ = z₁ with v ≠ 0. Then Π = W(v) commutes with V_M T₁. So
  Π U Π⁻¹ = ω^c U for U = W(a)V_M T₁, with c = ω(v, a) up to sign.
* **Proof.**
  * V_M W(v) V_M† = W(Mv) = W(v).
  * Mz₁ is collinear with z₁, so its x₁-component vanishes and W(v) has no X on qutrit 1. It therefore commutes
    with T₁, which is diagonal on qutrit 1.
  * Finally, W(v)W(a)W(v)⁻¹ = ω^{ω(v,a)}W(a).
* **Consequences when c ≠ 0** (exactly the frames F′ calls violating):
  * U maps the Π-eigenspace E_j to E_{j+c};
  * its spectrum is invariant under multiplication by ω;
  * tr(Π^m U^r) = 0 unless 3 divides r;
  * any reversal must carry Π to an element grading U⁻¹ with the same degree.
* **Checked** on all 1296 two-qutrit classes of both "bad" cells. Π commutes with V_M T₁ in all of them, and all
  69 984 frames with c ≠ 0 are graded, with ω-invariant spectrum.
* **What the lemma does not do.** The grading alone does not force violation, so F′ stays a conjecture. The cubic phase
  of T₁ cancels around each Π-cycle (Pass 11420, L2), so the missing step is not a scalar holonomy.

## Exhaustive check of the fixed-axis case at four qutrits

* **The fixed-axis cell is a group.** It is exactly H = Stab(z₁) ⊂ Sp(8,3), with |H| = 20 056 328 248 320 =
  |Sp(8,3)|/6560. Its H-conjugation orbits are the conjugacy classes of H.
* **GAP finds only 549 of them.**
* Every class representative was decided by the sparse decider of Pass 11421. F′ requires the violating frames to
  contain {a : a_{x₁} ≠ 0}.

| verdict | classes | share of H (by class size) |
|---|---|---|
| violating frames are exactly {a_{x₁} ≠ 0} | 224 | 95.37% |
| violating frames contain {a_{x₁} ≠ 0} | 8 | 2.23% |
| beyond the decider's cap (first pass, cap 3¹¹) | 317 | 2.40% |
| **F′ fails** | **0** | **0** |

* A second pass with cap 3¹⁴ on the 317 undecided classes ran more than 2.5 CPU-hours per worker without finishing and
  was stopped. These classes need the frame-orbit method of Pass 11373 with the Weyl criterion at n = 4; that is not
  done here. The producer keeps the second pass as an option, off by default, so the committed certificate is the
  default run.

**Reading.**
* At four qutrits the fixed-axis half of the law holds on every class the decider can reach, about 97.6% of H by
  class size, with no exception.
* The remaining classes are small: high-symmetry elements whose equation (S) has a very large solution space.
