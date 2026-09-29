# Pass 11124 — the neutral condensate unlocks the missing lepton and down-quark Yukawas, and nothing in the up sector

Producer: `analysis/w33_pass11124_condensate_unlocked_couplings.py`
Scan: `analysis/w33_pass11124_scan_dumps.py` (needs the Pass 11108 orbifolder dumps)
Frozen data: `data/w33_pass11124_{unlock_spacegroup,unlock_separate_winding,target_classes,rank_and_control}.json`
Regression: `tests/test_w33_pass11124_condensate_unlocked_couplings.py`

Method: the selection-rule engine (FastOrders MILP) was run for every up, down, lepton, neutrino-Dirac and μ coupling of
the six neutral-exit survivors. It was run once with the hidden-singlet SM-neutral VEV fields, and once with the tachyon T
(and T̄) added. T's quantum numbers:
* U(1) charges U·ℓ, exact;
* Witten sector k = 1, Z3 sector l = 0;
* R = 0 (NS ground state);
* its winding class N mod 3 in the space-group class of its own torus.

| model | lepton | down | ν Dirac | up | μ |
|---|---|---|---|---|---|
| 36621 | 27 | 0 | 54 | 0 | 0 |
| 46043 | 27 | 0 | 54 | 0 | 0 |
| 24165 | 81 | 81 | 162 | 0 | 0 |
| 40521 | 108 | 81 | 324 | 0 | 0 |
| 5904 | 81 | 81 | 270 | 0 | 0 |
| 17224 | 81 | 81 | 270 | 0 | 0 |

These entries are forbidden at every order without the condensate and allowed at orders 2–6 with it.

* **Mechanism.**
  * q_T lies outside the span of the singlet-VEV charges; the rank rises by exactly 1 in all six.
  * Every coupling target in span(S, q_T) \ span(S) needs exactly ±1·q_T (102 targets). So each unlocked Yukawa carries
    a single condensate insertion and is suppressed by ⟨T⟩/M_s.
* **Control.** A probe with charge 6·q_T is forbidden without T and allowed at order 6 with it, in all six models.
* **Up quarks and μ are untouched.** This matches the Pass 11117 argument: the condensate is a Δ(54) singlet, so
  m_c = m_u stays.
* **Validity of the winding rule.** 1854 of the 1863 entries contain a θ-twisted SM field, where winding counts only
  modulo (1−θ)Λ. The other 9 (40521, ν Dirac) are conditional.

**Correction made within this pass.** The first encoding treated T's winding as a separate, exactly conserved Z3 and found
0 unlocked couplings in all six. That was an artefact: winding modulo (1−θ)Λ *is* the fixed-point class, which twisted
fields compensate. The strict count (0) is kept as the bound for couplings with no twisted field. This matches the repo's
standing rule that a clean zero needs a named mechanism and a positive control. Here the control passed, and the
mechanism turned out to be the encoding.

Scope:
* Selection rules only: no coefficients and no mass eigenvalues.
* Whether the new charged-lepton entries lift the Δ(54) degeneracy m_μ = m_e is not decided here.
* The condensate VEV itself is not computed (Pass 11123: tree level does not fix it).
