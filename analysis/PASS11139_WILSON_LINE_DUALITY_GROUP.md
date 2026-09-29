# Pass 11139 — the duality group of the Wilson-line modulus is Γ0(3)

Producer: `analysis/w33_pass11139_wilson_line_duality_group.py`
Scan: `analysis/w33_pass11139_scan_duality_group.py`
Frozen: `data/w33_pass11139_duality_group_tests.json`
Regression: `tests/test_w33_pass11139_wilson_line_duality_group.py`

Z_β(τ; γT) was compared with Z_β(τ; T): the full Witten-sector partition function, with both Wilson-line tori at T. An
exact symmetry shows a deviation that shrinks as the coset cutoff K grows; a broken one does not.

| map | max deviation K=3 → K=4 | verdict |
|---|---|---|
| T → T + 1 | 0 → 0 | exact |
| T → T/(3T + 1) | ~5·10⁻⁶ → ~10⁻⁹ | **exact** (Γ0(3) generator) |
| T → (−T + 1)/(−3T + 2) | ~3·10⁻⁴ → ~3·10⁻⁷ | **exact** (order-3 elliptic element) |
| T → −1/T (S) | 3% | broken |
| T → T/(T + 1) (Γ⁰(3)) | 18% | broken |
| T → −1/(3T) (Fricke) | 0.34% (√3), 0.043% (golden), K-independent | approximate only |

The group is **Γ0(3)**, the level-3 congruence subgroup expected for Z3 Wilson lines. Its only elliptic point, of order 3,
is (3 + i√3)/6: **exactly the centre of the charged tachyon disk** of Pass 11136. The golden neutral disk is centred at
the Fricke point i/√3, which is a symmetry point only approximately. This matches Pass 11115's approximate Fricke duality,
now quantified.
