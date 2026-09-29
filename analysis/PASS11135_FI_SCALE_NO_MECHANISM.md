# Pass 11135 — the FI term would put the condensate in the right range, but no FI mechanism exists without supersymmetry

Producer: `analysis/w33_pass11135_fi_scale_no_mechanism.py`
Scan: `analysis/w33_pass11135_scan_fi.py`
Frozen: `data/w33_pass11135_fi_data_21.json`
Regression: `tests/test_w33_pass11135_fi_scale_no_mechanism.py`

**Numbers (21 neutral-exit survivors):**
* Each model has exactly one anomalous U(1).
* The tachyon's anomalous charge is |q_T,A| = 3 or 4.
* Tr Q_A = 12|u_A|² in all 21. This is a tautology of the orbifolder's definition t_A = (1/12)·Σ_f p_f, recorded as a
  consistency check.
* The supersymmetric formula ⟨φ⟩² = g²·Tr Q_A·M_P²/(192π²|q|) gives ⟨T⟩/M_P = 0.146–0.203 at g² = ½. That lies close to
  the ε ≈ 0.12–0.15 that Pass 11130 needs for m_b/m_t.

**Why this is not a mechanism.** The FI mass g²·q·ξ|φ|² belongs to the scalar of a *chiral multiplet*. For a complex
scalar without supersymmetry:
* nothing distinguishes φ from φ*; the dumps list both, e.g. Hd_j = conj(Hu_k) exactly (Pass 11130);
* CPT gives particle and antiparticle the same mass;
* a mass term linear in q is not covariant under φ ↔ φ*.

So no D-flatness condition fixes ⟨T⟩ in this vacuum. The scale coincidence is recorded and explicitly **not** claimed as
a prediction. The supersymmetric value also changes by √2 between the two standard normalisations of the FI term.
