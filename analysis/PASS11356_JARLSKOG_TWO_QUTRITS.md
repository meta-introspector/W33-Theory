# Pass 11356 — the substrate Jarlskog invariant on two qutrits: complete on every tick checked

Producer: `analysis/w33_pass11356_jarlskog_two_qutrits.py`
Certificate: `data/w33_pass11356_jarlskog_two_qutrits.json`
Regression: `tests/test_w33_pass11355_11360.py`

**Setup.** J₆(U) = avg_C |tr(U†·CUC†)|⁶ − avg_C |tr(U†·CUᵀC†)|⁶, averaged over the full two-qutrit Clifford group
(51 840 symplectic classes × 81 Pauli frames). See Pass 11355 for the definition and Theorem 1. It is tested against
exact verdicts:
* one-gate ticks W(a)·V_M·(T⊗I), with verdicts from the F₃-linear decider (Pass 11350), half drawn from bad classes;
* words with 2–4 cubic gates on random qutrits, with verdicts from the Weyl criterion (Pass 11252).

| | |
|---|---|
| ticks checked | 240 (173 violators) |
| **mismatches (J₆ > 0 ⇔ violating)** | **0** |
| largest \|J₆\| on reversible ticks | 7×10⁻¹⁵ |
| smallest J₆ on sampled violators | 0.0246 |

**Named ticks.**

| tick | J₆ | J₄, J₂ | verdict |
|---|---|---|---|
| (T⊗T)·SUM (minimal violator) | **0.0030** | 0, 0 | violating |
| (I⊗T)·SUM | 0 | 0, 0 | reversible |
| (I⊗T)·SUM·(I⊗T²) (the F_min² tick) | **2.26875 = 363/160** (numerical) | 0, 0 | violating |

**Reading.**
* The degree-six substrate Jarlskog invariant also separates reversible from violating ticks on two qutrits.
* The minimal two-qutrit violator has a small J₆. The deeper single-sector violator (I⊗T)·SUM·(I⊗T²) has a large
  one, even though its reversal fidelity is lower (Pass 11251).
* Completeness on two qutrits is verified on this sample, not proved.
