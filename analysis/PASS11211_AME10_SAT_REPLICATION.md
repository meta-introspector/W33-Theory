# Pass 11211 — independent replication of Pass 11190 (AME(10,3) sign structure) with two further SAT solvers

Producer: `analysis/w33_pass11211_ame10_sat_replication.py`
Frozen: `data/w33_pass11211_ame10_sat_replication.json`
Regression: `tests/test_w33_pass11211_ame10_sat_replication.py`

**Numbering and scope.** This work was done independently on branch `claude/gallant-sagan-i33c7e` as "Pass 11190"
(19:25 UTC on 2026-09-30), in parallel with master's Pass 11190 (`analysis/PASS11190_AME10_SIGN_STRUCTURE.md`, 20:49 UTC).
Master's version is canonical and more complete. It adds:
* the uniqueness of the sign structure, with automorphism group PGL(2,9);
* a direct combinatorial exclusion of C4+C6 (a prism has no 4-cycle of the required kind);
* the explicit S(3,4,10) ⇔ block correspondence.

What remains here is a **replication with a different encoding and different solvers**.

**Same mathematics, independently derived.** The derivations match master's:
* sign partitions of the 6-sets (uniform or 3 + 3);
* Klein-quadric realisability of the local sign graphs (7K1 and K33+K1 elliptic; prism+K1 and K4+3K1 hyperbolic);
* the isotropy duality;
* the dictionary det S[o, i] = −1 ⇔ i and o lie in one sign class of {i} ∪ O.

**Replication.**
* The model is a boolean CNF with one-hot sign states on the 210 six-sets, 455 allowed labelled local graphs per 7-set,
  and 1575 duality pairs: 56 910 variables and 488 580 clauses.
* It was solved with **CaDiCaL 1.5.3** and **Glucose 4** (PySAT); master used OR-tools CP-SAT. Both solvers give:
  * K4+3K1 on {0, …, 6}: UNSAT;
  * any pattern other than C10 or cross+perm at inputs {0, …, 4}: UNSAT;
  * controls (C10 and cross+perm each realisable): SAT.
* So three independent solver families agree on master's two infeasibility results.
* Positive control: all 71 tabu graphs satisfy every constraint; the dictionary holds on 12 600 of 12 600 blocks; and the
  census is 5112 C10 / 12 780 cross+perm.

No new mathematical claim is made here; the theorem and its credit belong to master's Pass 11190.
