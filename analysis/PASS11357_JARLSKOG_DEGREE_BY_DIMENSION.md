# Pass 11357 — the degree of the time-reversal witness depends on the dimension: 10 for qubits, 6 for qutrits, 4 for d = 5

Producer: `analysis/w33_pass11357_jarlskog_degree_by_dimension.py`
Certificate: `data/w33_pass11357_jarlskog_degree_by_dimension.json`
Regression: `tests/test_w33_pass11355_11360.py`

**Setup.**
* Single qudits of dimension d = 2, 3, 5, each with its Clifford group modulo phase (orders 24, 216, 3000; generated
  by the Fourier, phase and shift gates).
* A diagonal magic gate: the qubit T = diag(1, e^{iπ/4}), or the cubic phase diag(ζ^{k³}).
* Exact verdicts by searching every anti-unitary Clifford reversal.

**Which witnesses vanish identically** (polynomial identity test on 40 Haar-random unitaries).

| d | J₂ | J₄ | J₆ | J₈ | J₁₀ | first non-zero degree |
|---|---|---|---|---|---|---|
| 2 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | 1.73 | **10** |
| 3 | ≈ 0 | ≈ 0 | 1.72 | 23.9 | 255 | **6** |
| 5 | ≈ 0 | 0.14 | 4.82 | 128 | 3243 | **4** |

The table shows the maximum |J_{2t}| over the 40 Haar samples. "≈ 0" means below 10⁻¹³, at rounding level.

**Completeness on Clifford+magic words** (720 words per d, 1–6 magic gates):

| d | violators | first complete witness |
|---|---|---|
| 2 | 170 | **J₁₀** (J₂–J₈ detect none) |
| 3 | 418 | **J₆** |
| 5 | 533 | **J₆** (J₄ detects 517 of 533) |

**Correction to Pass 11355.**
* J₂ ≡ 0 follows from the 2-design property.
* J₄ does **not** follow from it. The degree-(t,t) invariants come from the commutant of C^{⊗t} ⊗ C̄^{⊗t}, whose
  dimension is the frame potential at 2t. d = 5 is a 2-design with J₄ not identically zero.
* For qutrits, J₄ ≡ 0 holds as a verified identity.

**Reading.**
* Polynomial visibility of magic-induced time-reversal violation needs degree 10 for qubits, 6 for qutrits, and 4 (to
  appear) or 6 (to be complete) for d = 5.
* A sketch for qubits: for SU(2), Uᵀ = σ_y U† σ_y with σ_y a Pauli, so qubit reversibility asks whether U⁻¹ is
  Clifford-conjugate to U. The low-degree real invariants cannot see that. This is a mechanism sketch, not a proof.
* The qutrit witness is the cubic one, matching the degree of the quark-sector Jarlskog invariant.
