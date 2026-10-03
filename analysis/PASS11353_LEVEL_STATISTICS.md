# Pass 11353 — the substrate arrow is visible in level statistics: reversible magic is COE, violating magic is CUE

Producer: `analysis/w33_pass11353_level_statistics.py`
Certificate: `data/w33_pass11353_level_statistics.json`
Regression: `tests/test_w33_pass11350_11354.py`

**Setup.** Floquet ticks on n qutrits (dimension 3ⁿ):
* W = a Clifford+T circuit: 12 (n = 4) or 15 (n = 5) layers, each a random Clifford word followed by a cubic gate on a
  random qutrit.
* Ensembles:
  * **violating:** U = W. Long magic circuits are reversible with probability → 0 (Pass 11332).
  * **reversible:** U = W·Wᵀ. Then the complex conjugation K, an anti-unitary Clifford, reverses it: Ū = U⁻¹.
  * **Clifford only:** U a random Clifford.
* Statistic: the mean spacing ratio ⟨r⟩ of eigenphases. Random-matrix reference values (Atas et al., PRL 110, 084101,
  2013): Poisson 0.3863, COE 0.5307, CUE 0.5996, CSE 0.6744.

| ensemble | ⟨r⟩, n = 4 (81 levels × 40) | ⟨r⟩, n = 5 (243 levels × 40) | class |
|---|---|---|---|
| violating magic tick | 0.603 ± 0.005 | 0.600 ± 0.002 | **CUE** |
| reversible magic tick | 0.529 ± 0.005 | 0.536 ± 0.003 | **COE** |
| Clifford only | 42% of spacings degenerate | 68% degenerate | no level repulsion (finite order) |

**Why the symplectic class cannot occur.**
* An anti-unitary Θ on an odd-dimensional space cannot square to −1, because Kramers pairs need even dimension.
* Every substrate time reversal therefore has Θ² = +1. Checked: all 40 anti-unitary Clifford reversals recovered for
  two-qutrit magic ticks have V·V̄ = +I.
* So reversible magic dynamics on qutrits is **orthogonal-class only**, never symplectic-class.

**Reading.**
* The combinatorial arrow of Passes 11252–11352 (does an anti-unitary Clifford reverse the tick?) has a spectral
  fingerprint: COE level repulsion when the arrow is absent, CUE when it is broken.
* Clifford dynamics shows neither, because its spectrum is degenerate.
* The random-matrix mechanism is Dyson's threefold way and is standard. What is specific here is the realisation on
  the substrate's gate set and the exclusion of the symplectic class by odd dimension.
* Scope: ensemble averages. A single tick's ⟨r⟩ fluctuates; this is not a per-tick reversibility test.
