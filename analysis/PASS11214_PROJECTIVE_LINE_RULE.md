# Pass 11214 — the AME sign structures are cross-ratio geometry: harmonic pairs on PG(1,5), chiral cross-ratios on PG(1,9)

Producer: `analysis/w33_pass11214_projective_line_rule.py`
Certificate: `data/w33_pass11214_projective_line_rule.json`
Regression: `tests/test_w33_pass11213_11216.py`

**Background.** Passes 11190 and 11200 found that the sign structures of AME(6,3) and AME(10,3) stabilizer states are
unique, with symmetry PGL(2,5) and PGL(2,9); Pass 11211 (merged) replicates 11190. This pass writes the structures down
on the projective lines.

**Six qutrits: harmonic conjugacy on PG(1,5).** Label the parties by PG(1,5).
* Every 4-set of parties splits 2|2, and **the split is its harmonic pairing**: the unique pairing {a,b}|{c,d} with
  cross-ratio (a,b;c,d) = −1.
* In each of 12 samples, exactly |PGL(2,5)| = 120 of the 720 labellings realise this.

**Ten qutrits: chiral cross-ratio geometry on PG(1,9)**, with F9 = F3[i] and i² = −1.
* The uniform 4-sets (the Steiner system S(3,4,10)) are exactly the **Baer sublines**, the harmonic quadruples with
  cross-ratio in F3. Every state has exactly 2 such labellings once three points are fixed, a Frobenius pair.
* Take a 4-set K that is not a subline. Its three pairings have cross-ratios in the three classes {±i}, {1+i, −1+i}
  and {1−i, −1−i} of F9 ∖ F3.
* **Rule.** Two points x, y outside K are in the same sign class of the 6-set K^c iff both hold:
  * {x, y} is harmonic to exactly one pair {a, b} of K;
  * the pairing {a,b}|{c,d} of K is not in the excluded class.
* The excluded class is {1+i, −1+i} for one labelling of the Frobenius pair and its conjugate {1−i, −1−i} for the
  other. This holds for **all 71 states, with no exception** (44 of one handedness relative to the normalisation, 27 of
  the other).
* Before the rule was found, the refined statistics had already shown it to be decisive: 720 / 720 / 360 pairs, no
  exceptions.

**Reading.** The sign structure tells a cross-ratio from its Galois conjugate, so it is chiral. That is exactly why its
symmetry is PGL(2,9) and not PΓL(2,9), and why a Steiner system carries two sign structures (Pass 11190), the Frobenius
mirrors of each other. Six qutrits are harmonic conjugacy over F5, which is achiral since F5 has no Frobenius. Ten
qutrits are cross-ratio geometry over F9 with a handedness.

**Correction during the pass.** A first exploration encoded F9 with index 1 = i instead of 1. The rule it suggested was
structured but mislabelled. The producer uses a checked field class with ONE = 1.
