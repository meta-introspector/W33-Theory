# Pass 11117 — the neutral condensate breaks exactly one U(1): the Standard Model and the family symmetry survive

Producer: `analysis/w33_pass11117_neutral_condensate.py`
Data: `data/w33_pass11117_neutral_condensate_charges.json` (explicit enumerator of Pass 11107; the field-dump U(1) basis)
Certificate: `data/w33_pass11117_neutral_condensate.json`
Regression: `tests/test_w33_pass11117_neutral_condensate.py`

Pass 11112 found that in models 53 and 57 the one-loop potential reaches a Standard-Model-neutral winding tachyon
first. This pass determines what its condensate breaks. The setup is Wilson-line torus 1 at 1.4i, torus 2 at 3i and the
family torus at ρ.

| | model 57 | model 53 |
|---|---|---|
| U(1)s in the model | 8 | 9 |
| hypercharge of the tachyons | 0 | 0 |
| SU(3)_c, SU(2)_L, hidden non-abelian | singlets | singlets |
| distinct U(1) charge vectors | ±one vector | ±one vector |
| U(1) directions broken by the condensate | **1** | **1** |
| charged under the anomalous U(1) | yes (±4) | yes (±4) |
| cubic self-coupling allowed | no | no |

## Reading

* **Exactly one U(1) combination is broken.** It has a component along the anomalous U(1), which is already massive by
  the Green–Schwarz mechanism, and components along non-anomalous U(1)s. Hypercharge, colour, SU(2) and every hidden
  non-abelian factor remain unbroken, so electromagnetism survives.
* **The flavour symmetry is untouched.** The states wind in a Wilson-line torus and carry no fixed-point label. The
  physical state is the Z3-invariant combination of winding images, a singlet of the family Δ(54). The condensate
  therefore cannot lift m_c = m_u.
* **No cubic term.** The only gauge-neutral combinations are even in the charge. The condensation starts as a
  second-order, Higgs-like instability, not a first-order jump.

Not computed: the quartic coefficient, which needs string four-point amplitudes of the winding states, and the
back-reaction on the moduli. Together these decide where the condensation ends.
