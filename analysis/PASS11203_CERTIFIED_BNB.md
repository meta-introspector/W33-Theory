# Pass 11203 — the torus-class theorem 4/√15 is now rounding-certified

Producer: `analysis/w33_pass11203_certified_bnb.py`
Certificate: `data/w33_pass11203_certified_bnb.json`
Regression: `tests/test_w33_pass11203_certified_bnb.py`

**Why.** Pass 11187's branch and bound used float64 with a 10⁻⁹ margin, and the Pass 11188–11192 prior-art sweep
flagged it as not using directed rounding.

**Method.**
* IEEE-754 +, −, ×, ÷ and √ are correctly rounded. Nudging every intermediate result one ulp outward with `nextafter`
  therefore gives a guaranteed enclosure.
* Two closing rules are used, both fully outward-rounded:
  1. the monotone bound (A increasing in p and c, decreasing in q) is below 4/√15, itself rounded down;
  2. the Pass 11187 interval-Hessian concavity test, on 6³ cells whose endpoints are pinned so that they provably cover
     the hull, after which the Pass 11167 slice theorem applies unchanged.
* The non-rigorous centred form is dropped.
* An infeasible box is discarded only if its rounded-down Σlo exceeds 1. Split midpoints are shared exactly.

**Result.**
* The search completes: **9 219 253 boxes**, of which 4 384 961 are closed by the rigorous bound, 224 627 by the
  rigorous concavity core, and 39 are infeasible. It took 1406 s.
* Controls:
  * A_up encloses the exact value at 300 random points (mpmath, 50 digits);
  * the rigorous Hessian test still rejects a non-concave cell and accepts the optimum's neighbourhood;
  * the target is provably below 4/√15.

**Status.** "4/√15 is the maximum of N_AB + N_AC on the torus-symmetric class" is now a computer-assisted theorem with
certified rounding. The only remaining assumption is the proof in the Pass 11167 slice theorem. The global problem is
the Allen–Meyer (2017) conjecture specialised to the sum, and remains open.
