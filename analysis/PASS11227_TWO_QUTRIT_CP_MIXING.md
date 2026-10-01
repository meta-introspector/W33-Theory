# Pass 11227 — two qutrits: entanglement supplies the mixing, both sectors need the phase, and the quantum is 2π/9

Producer: `analysis/w33_pass11227_two_qutrit_cp_mixing.py`
Certificate: `data/w33_pass11227_two_qutrit_cp_mixing.json`
Regression: `tests/test_w33_pass11227_11228.py`

**Question (from Pass 11213).** On one qutrit, the minimal ticks that no substrate time reversal inverts are the cubic
phase T = ζ₉^{x³} composed with a cyclic shift X^a of the three levels. Can the entangling gate SUM, a shift of one
qutrit controlled by the other, supply that mixing instead?

**Method.** The test is exact. U is substrate-T-reversible iff some anti-unitary two-qutrit Clifford Θ = V K inverts it,
i.e. V U* V† U ∝ I. The search runs over all 51 840 symplectic classes × 81 Paulis, with representatives generated from
H, S and SUM. The best time-reversal fidelity is F_T(U) = max_V |tr(V U* V† U)|/9 (Pass 11228); F_T = 1 iff reversible.

**Result.**

| circuit | substrate time reversal? |
|---|---|
| T ⊗ T (cubic phases only) | reversible |
| SUM (Clifford only) | reversible |
| (I⊗T)·SUM, (T⊗I)·SUM | reversible |
| SUM·(I⊗T)·SUM†, SUM·(I⊗T)·SUM | reversible |
| (I⊗T)·SUM·(I⊗T), with two cubic gates on the same qutrit | reversible |
| **(T⊗T)·SUM** | **T-violating** |
| **(T⊗T)·SUM·(T⊗T)·SUM** | **T-violating** |
| T·X ⊗ I (control: one-qutrit violator) | T-violating |

* **Entanglement replaces the explicit shift, but only when both qutrits carry a cubic phase.** One cubic phase with an
  entangler, in any placement tried, is reversible. So are two cubic phases on the same qutrit.
* **The violation is quantised and universal.** Every violating circuit tested has the same best time-reversal fidelity,
  F_T = (1 + 2cos 2π/9)/3 = 0.8440296287…. That is the same value as the one-qutrit minimal violators. The circuits
  tested are (T⊗T)·SUM, (T⊗T†)·SUM, (T†⊗T)·SUM, (T²⊗T)·SUM, SUM·(T⊗T), (T⊗T)·SUM†, the doubled circuit and T·X⊗I.
  * Opposite phases on the two qutrits do **not** cancel.
  * Two layers do **not** add.
  * The irreducible T-odd phase is 2π/9.

**Reading (analogy, not derivation).** In the Standard Model, CP violation needs phases in *both* quark sectors that
mixing cannot simultaneously remove: the Jarlskog invariant contains the mass differences of the up *and* down sectors.
On the substrate, T-violation needs cubic phases on *both* qutrits, linked by an entangler. One sector alone, however
dressed, can always be undone by a substrate time reversal. No CKM matrix, mass or measured phase is derived.

**Scope.** Ten circuits checked by exact search, plus eight fidelity evaluations. Random two-qutrit Clifford+T words
and larger registers are not sampled here.
