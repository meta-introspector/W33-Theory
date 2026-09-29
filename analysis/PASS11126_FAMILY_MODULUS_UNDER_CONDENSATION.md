# Pass 11126 — the neutral condensate cannot move the family modulus, so it cannot touch the family hierarchy

Producer: `analysis/w33_pass11126_family_modulus_under_condensation.py`
Data: `data/w33_pass11126_family_torus_quantum_numbers.json`
Certificate: `data/w33_pass11126_family_modulus_under_condensation.json`
Regression: `tests/test_w33_pass11126_family_modulus_under_condensation.py`

The tachyons of the six neutral-exit survivors were checked in four places: just below the diagonal onset, at
Im T_WL = 1, and on each Wilson-line torus shrunk alone. Every one of them has **zero momentum and winding on the family
torus**, the one without Wilson lines. Each model has exactly one complex pair (4 real states), living on one Wilson-line
torus.

Consequences:
* **No tree-level coupling to T\*.** The condensate's vertex operators carry no family-torus lattice momentum. T\* enters
  the potential only through the Narain factor Γ₂,₂(T\*, ρ) of the Witten sectors (Pass 11106). The one-loop potential
  therefore stays exactly SL(2,Z)-invariant in T\*, with critical points at i and ρ.
* **The flavour structure is untouched.** The condensate is a Δ(54) singlet (Pass 11117), so m_c = m_u and the tree-level
  ratio m_{c,u}/m_t = |Y_dist/Y_same|(T\*) (Pass 11109) do not change. Even a jump of T\* between the two fixed points
  gives ½ at ρ or 0.366 at i, never the 0.0036 the data need (which requires Im T\* ≈ 3.2).

The condensate lives in the Wilson-line sector; the family problem must be solved elsewhere.
