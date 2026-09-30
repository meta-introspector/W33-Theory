# Pass 11199 — the local geometry of AME(10,3) is a regular spread, but its reguli are not the Steiner blocks

Producer: `analysis/w33_pass11199_ame10_local_spreads.py`
Certificate: `data/w33_pass11199_ame10_local_spreads.json`
Regression: `tests/test_w33_pass11199_ame10_local_spreads.py`

**Question (open in Pass 11190).** The uniform 4-sets of an AME(10,3) stabilizer state form S(3,4,10), the inversive
plane of order 3, which is the plane of reguli of a regular spread of PG(3,3) (Thas 1997; the repo's regular-spread
packets are BT2088–2092). Are the ten parties the ten lines of such a spread, with the blocks as its reguli?

**Test.** For every 3-set A (71 graphs × 120 = 8520 cases):
* W_A (the stabilisers trivial on A) is 4-dim.
* The seven parties outside A give seven kernel lines in PG(W_A) = PG(3,3).

**Result (all 8520 cases).**
* The seven lines complete to a spread in exactly one way, and that spread is **regular**. Every 3-set of parties
  therefore carries a canonical regular spread of PG(3,3), i.e. a local elliptic quadric on the Klein quadric.
* **No labelling** of the three missing lines by the parties of A makes all 30 Steiner blocks reguli.
* Among the seven kernel lines there are exactly 5 reguli. Exactly 2 of them are Steiner blocks. The other 3 are 4-sets
  whose complementary 6-set splits 3|3 in a way that divides A between the two classes. No regulus has A as a whole class.
* Of the 5 Steiner blocks inside the seven, the other 3 are not reguli.

**Answer.** Negative in the strong form. Each local geometry is a regular spread, so each carries its own Miquelian
inversive plane, but that plane is not the state's Steiner system: the two meet in exactly 2 of 5 local circles.

**Scope.** Stabilizer AME(10,3) states. Whether a global spread, in some other space, realises the blocks remains open.
