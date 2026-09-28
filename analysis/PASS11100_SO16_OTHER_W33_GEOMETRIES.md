# Pass 11100 — SO(16)×SO(16) with the W(3,3) shifts on the other geometries: only T⁶/Z3 works

Producer: `analysis/w33_pass11100_so16_other_w33_geometries.py`
Frozen evidence: `data/w33_pass11100_so16_other_geometries_evidence.json` (every base with its Witten shift V₀, and
the per-base scan counts)
Certificate: `data/w33_pass11100_so16_other_w33_geometries.json`
Regression: `tests/test_w33_pass11100_so16_other_w33_geometries.py`

Pass 11095 built 104 tachyon-free three-generation SO(16)×SO(16) models on the W(3,3) A8 shift of T⁶/Z3. This pass
runs the same construction on the other geometries that carry W(3,3) shifts. The W(3,3) shifts are held fixed.
For each base a Witten shift V₀ = (a; b) is chosen, with a, b ∈ ½·{norm-4 E8 vectors} and V₀·V_k ∈ Z for every
shift: standard-like candidates are tried first, and the first one the non-SUSY orbifolder loads is used. Every base
got a loadable V₀. The Wilson lines are then drawn at random, with the non-SUSY orbifolder's tachyon and
Standard-Model tests.

| geometry | W(3,3) bases | draws | tachyon-free | SM-like draws | tachyon-free **and** SM-like |
|---|---|---|---|---|---|
| T⁶/Z3 (Pass 11095) | A8 pair | 80,000 | **100%** | 153 | **153** |
| T⁶/(Z3×Z3) | 2 | 13,500 | 0.7% | 12 | 0 |
| T⁶/Z6-I | 58 | 29,000 | 38% | 3 | 0 |
| T⁶/Z12-I | 117 | 35,100 | 4.2% | 0 | 0 |

**Only the prime geometry T⁶/Z3 gives tachyon-free Standard-Model-like models.** There, every one of 80,000 draws
was tachyon-free. In the other geometries tachyons are common, and the few Standard-Model-like draws (15 in all)
are all tachyonic.

Bookkeeping. The Z3×Z3 input directory held four identical copies of each of its two bases (c1 to c4). The first
pass therefore counted every distinct draw four times. The table uses the 2 distinct bases: 3,000 draws each with
seed 20260928, plus a rescan of 3,500 and 4,000 draws with seed 20260929.

Scope. This is a sample, with one V₀ per base and random Wilson lines, and it proves nothing. It does not say the
other geometries cannot host a tachyon-free W(3,3) Standard Model. It says that on this sample the T⁶/Z3
construction of Pass 11095 is the one that works.
