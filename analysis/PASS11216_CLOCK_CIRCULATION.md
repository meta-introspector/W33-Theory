# Pass 11216 — the one-way clock circulates: current spectra and oriented triangles

Producer: `analysis/w33_pass11216_clock_circulation.py` (reads the intersection matrices frozen in Pass 11205)
Certificate: `data/w33_pass11216_clock_circulation.json`
Regression: `tests/test_w33_pass11213_11216.py`

**Result.**
* **Current operator.** For the two pairs where A_R and A_{Rᵀ} commute, J_R = A_R − A_{Rᵀ} has eigenvalues exactly
  2i·Im λ:
  * **±24√3·i and ±36√3·i** on the 256 pair;
  * **±72√3·i and ±324√3·i** on the 6912 pair.

  These are integer multiples of √3·i, the qutrit phase scale. The 2304 pair does not commute, and its current
  eigenvalue (≈ 185.9) is not of that form.
* **Oriented triangles.** Through each split pass 7168 oriented 3-cycles of the 256 relation (x → z → y → x), against
  2304 transitive triangles. For the 6912 relation the counts are 3 462 912 against 2 903 040. Cyclic, one-way flow
  dominates.
* **Relation to chirality.** The Kashiwara–Maslov chirality of Pass 11209 (merged from the parallel session) is ∓54
  on the 256 pair and ±18 on the 6912 pair. No simple relation between these numbers and the circulation data is
  found here, and none is claimed.

**Scope.** Spectral and counting data of a finite relation algebra. No physical time scale is assigned.
