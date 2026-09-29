# Pass 11131 — across all 21 neutral-exit survivors the condensate never unlocks an up-type Yukawa or a μ-term

Producer: `analysis/w33_pass11131_unlocked_across_21.py`
Scan: `analysis/w33_pass11131_scan_unlock_21.py` (hidden-gauge check enforced; neutral state selected independently)
Frozen: `data/w33_pass11131_unlock_scan_21.json`
Regression: `tests/test_w33_pass11131_unlocked_across_21.py`

| sector | models with unlocked entries (of 21) |
|---|---|
| up-type Yukawa | **0** |
| μ-term | **0** |
| down-quark | 14 |
| charged lepton | 17 |
| neutrino Dirac | 20 |

* A8SM_20260982_47048 unlocks nothing at all.
* The hidden-gauge check removes no entry in any model.
* The positive control passes in all 21: a 6·q_T probe is forbidden without T and allowed at order 6 with it.
* The six of Pass 11124 reproduce its counts exactly, with the neutral state chosen independently. This is an independent
  re-derivation of 11124's counts.

Reading:
* "Up and μ untouched" holds across the whole class.
* The single-Higgs lowering of m_b/m_t from ε³ to ε² (Pass 11130) needs unlocked down entries, which 14 models have. It
  is verified only in four of the six.
