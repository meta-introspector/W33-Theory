# Pass 11165 — a perfect three-qutrit tick secret-shares every past: forgetting ANY one future erases it

Producer: `analysis/w33_pass11165_three_qutrit_arrow.py`
Regression: `tests/test_w33_pass11165_three_qutrit_arrow.py`

**Formula.** For an n-qutrit Clifford gate S we use the stabilizer formula of Passes 11156/11161:
S(X) = |X| − 2n + rank(G restricted to the legs outside X), with G = [I; S].
* Pass 11166 checks this formula against brute-force Choi entropies.
* For each input leg and each set Y of output legs we compute I(A_in : Y).

**Results.** The sample is 12 perfect gates (the F₉ gate K = P² + i𝟙, the Pass 11158 example, and random ones) and
12 non-perfect controls.
* **Perfect gates:** I(A_in : Y) = 0 for every proper subset Y of the three outputs, and 2 for all three together.
  Each input's past is shared among the three futures with threshold (3,3).
* **Controls:** all 12 non-perfect gates leak, meaning some proper subset of outputs keeps information.

**Reading.**
* Three qutrits answer the question "does forgetting two partners cost 4 trits?": **no**.
* The two-qutrit arrow (Pass 11161) needed the partner discarded. In a perfect three-qutrit tick, forgetting **any
  one** future already erases the past completely.
* The erased record is still exactly two trits, since that is the maximum mutual information of one qutrit with
  anything. The arrow's price is intrinsic to one qutrit, not to the number of partners.

**Prior art (this property is known).** Hiding an input from every proper subset of outputs is the secret-sharing
property of AME states (Helwig, Cui, Latorre, Riera, Lo, PRA 86, 052335 (2012)). Pahari (arXiv:2607.00210) calls it
exact local masking and studies it for two-qudit brickwork circuits. New here is only the reading in the paper's
arrow of time, and the contrast with the paper's own clock (Pass 11166).
