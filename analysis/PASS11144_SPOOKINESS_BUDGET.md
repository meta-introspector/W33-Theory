# Pass 11144 — the spookiness budget: W(3,3) is classical in space and in time; spooky action needs magic

Producer: `analysis/w33_pass11144_spookiness_budget.py`
Frozen optimum: `data/w33_pass11144_temporal_optimum_params.json`
Certificate: `data/w33_pass11144_spookiness_budget.json`
Regression: `tests/test_w33_pass11144_spookiness_budget.py`

The yardstick is the qutrit Bell inequality of Collins–Gisin–Linden–Massar–Popescu (CGLMP; two settings, three
outcomes). Local hidden variables give at most **2**.

| scenario | max I | status |
|---|---|---|
| local hidden variables (all 81 strategies) | 2 | exact |
| **space, W(3,3) resources only** (stabilizer states, Pauli bases, all labellings) | **2** | exact: no violation |
| **time, W(3,3) resources only** (stabilizer start, Pauli bases, Lüders updates) | **2** | exact: no violation |
| space, maximally entangled pair, CGLMP bases | 2.872934 | Collins et al. |
| **time, one qutrit, maximally mixed start** | **2.872934** | identical to the Bell pair |
| space, optimal state (Acín–Durt–Gisin–Latorre) | 2.914854 | spatial quantum optimum |
| **time, one qutrit, optimal start and bases** | **3.1628065** | 40 restarts; no closed form (not √10) |
| time, classical system with ≥ 2 memory states | 4 | algebraic maximum |

Findings:
1. **Inside W(3,3) there is no spooky action, in space or in time.** In odd dimension, stabilizer states and Pauli
   measurements have a non-negative discrete Wigner function (Gross). Its phase space F₃² × F₃² is the point space of
   W(3,3), and that is a local hidden-variable model.
2. **W(3,3) sits exactly on the edge.** Along the path from the Pauli bases (t = 0) to the optimal CGLMP bases (t = 1),
   I(t) = 2 + (π²/6)·t² + O(t³). Any non-Clifford rotation in this direction violates, starting at second order. The
   resource that turns W(3,3) correlations into Bell violations is magic.
3. **One qutrit in time is a Bell pair in space.** Started maximally mixed, a qutrit measured twice gives exactly the
   Bell-pair statistics (2.872934). This is Pass 11143's state ↔ gate bijection, now at the level of measured data.
4. **For qutrits the space–time equality breaks.** For qubits the temporal CHSH optimum equals the spatial one, 2√2
   (Fritz; reproduced here as a method check). For qutrit CGLMP the temporal optimum, 3.1628, *exceeds* the spatial
   optimum, 2.9149. The optimal start is strongly biased: its first-outcome probabilities are 0.758, 0.121, 0.121.
5. **The orderings invert.**
   * In space: classical 2 < quantum 2.915.
   * In time: quantum 3.163 < classical-with-memory 4. Measurement disturbance (non-orthogonal post-measurement states)
     caps how much a quantum system can "tell its future self".

**What this says about spooky action at a distance.**
* The correlation geometry of a Bell pair is exactly the geometry of one qutrit correlated with itself across time,
  related by a partial transpose (11143).
* The temporal version needs nothing but persistence, and its larger values measure invasiveness.
* The spatial version is spooky precisely because it cannot signal, and it needs magic, a resource outside W(3,3).
* So the "action at a distance" is the no-signalling shadow of a self-correlation in time. That is a statement about the
  structure of the correlations, not a mechanism: nothing travels backwards in time.

Prior art:
* Fritz, NJP 12, 083055 (temporal CHSH);
* Budroni–Emary, arXiv:1309.3678 (temporal correlations grow with dimension);
* qudit CHSH and Wigner negativity, arXiv:2405.14367;
* Collins et al. and Acín et al. (CGLMP values).

The temporal CGLMP optimum for one qutrit (3.1628), the exact W(3,3) = 2 statements, and the π²/6 onset were not found in
these sources.
