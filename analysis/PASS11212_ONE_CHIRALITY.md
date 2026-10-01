# Pass 11212 — one chirality: the Gauss-sum twist and the Maslov chirality agree, and a label-free Choi phase exists

Producer: `analysis/w33_pass11212_one_chirality.py`
Certificate: `data/w33_pass11212_one_chirality.json`
Regression: `tests/test_w33_pass11212_one_chirality.py`

**Background.** Two τ-odd ("chiral") invariants of three-qutrit split relations were found in parallel. The two
sessions work in parallel, so this pass compares their invariants directly rather than re-deriving either.

* **Master's Pass 11201** uses the third-moment invariant I(σ) = tr[ρ^{⊗3} ∏_p P_{σ_p}] of the gate's Choi state,
  an exact quadratic Gauss sum.
  * For one labelled pattern (CHIRAL) it equals Q·i√3/243 on the 256 pair, where Q is the twist of master's Pass 11189.
  * For the 6912 pair a second pattern gave ±i√3/729 on the representatives.
  * Its stated caveat: the sign needs registers labelled compatibly with the relation.
* **Pass 11209** counts class-1 minus class-3 triples of product Lagrangians (Kashiwara–Maslov index, oriented by a
  gauge-invariant profile).
  * The result is χ = ∓54 on the 256 pair and ±18 on the 6912 pair.
  * It sums over all product Lagrangians of both splits, so it needs no labelling.

## Results

**1. One chirality on the 256 pair: χ = 54·Q.**
* On **all 512 members** of the two 256 orbitals, the Maslov chirality, the Gauss-sum twist and the Choi-state phase
  agree: χ = 54Q and Im I(CHIRAL) = Q·√3/243.
  * Orbital with Q = −1: χ = −54, Im = −√3/243, on 256 of 256 members.
  * Orbital with Q = +1: χ = +54, Im = +√3/243, on 256 of 256 members.
* Hence χ = 54 · (243/√3) · Im I(CHIRAL) exactly on this pair.
* Master's certificate and Pass 11209 therefore describe the same orientation, reached two ways:
  * from the Weil representation, through the Choi state's third moment;
  * from the Maslov index, through the Bargmann phases of product stabilizer states.

**2. On the 6912 pair, Pass 11201's pattern has the right sign but is not an invariant.**
* Pass 11201 evaluated its 6912 pattern on one representative per orbital. Members of the same orbital carry different
  labellings of the second split's factors.
* Sampling 64 members per orbital gives these values of the pattern (in units of i√3/729):

  | 6912 orbital | χ | values (count) |
  |---|---|---|
  | chiral | +18 | −1 (14), 0 (50) |
  | chiral | −18 | +1 (22), 0 (42) |
  | achiral | 0 | −1 (13), 0 (42), +1 (9) |
  | achiral | 0 | −1 (22), 0 (24), +1 (18) |
  | achiral | 0 | −1 (34), +1 (30) |
  | achiral | 0 | 0 (64) |

* So on the chiral pair the value **vanishes on most labellings**, but where it is nonzero its sign is set by the
  relation, −sign χ.
* On achiral orbitals both signs occur. A nonzero value there carries no orientation.
* The sign reported in Pass 11201 is therefore correct for the relation. The value itself is a property of the
  labelling, and a single member's zero or nonzero value proves nothing.
* The relation-level invariants of that pair are χ = ±18 (Pass 11209) and master's SL(2,3) holonomy (N vs −N,
  Pass 11189), plus the label-free tagged phase below.
* This is a scope refinement to one line of Pass 11201. The 256-pair statement there is confirmed on every member.

**3. A label-free Choi-state phase.**
* Tag each of the 36 relabellings (π_in, π_out) of the six registers by the block-rank code matrix of the relation,
  read in the pattern's slots.
  * The codes are: 0 for a zero block, 1 for rank one, 2 or 3 for an invertible block with det 1 or 2.
  * They are local-Clifford invariant and τ-blind.
* The multiset of (tag, I) over all 36 relabellings is independent of how the relation's factors are labelled. It is
  constant on orbitals and complex-conjugated by τ.
* It is **chiral exactly on the four τ-swapped relations** (the 256 and 6912 pairs) and achiral on the other 16. This
  holds for both patterns (CHIRAL and Pass 11201's 6912 pattern).
* It separates all 20 relations.
* So the Choi-state phase does measure the direction of time on *unlabelled* registers, once the relabellings are sorted
  by block type. This is the same device, a gauge-invariant orientation, that makes the Maslov count label-free.
* This refines Pass 11201's caveat. Plain averaging over relabellings erases the sign, as 11201 found; sorting by block
  type keeps it.

**Correction during the pass.** A first version read the code matrix as [input, output]. Its rows are actually the
outputs (the factors of F₀) and its columns the inputs (the factors of F). That tag was not covariant, and one 6912
orbital came out non-constant. Direct tests confirmed that the Gauss-sum evaluator is local-Clifford invariant and
matches direct contraction, which isolated the error. After the fix the invariant is constant on every orbital.

## Scope

* The comparison is exact on the members checked:
  * all 512 for the 256 pair;
  * 64 per 6912 orbital for the labelling test, with χ recomputed on 4 of them;
  * 4 per orbital, including the stored representative, for the tagged multiset.
* The regression recomputes on stored representatives without the 110 565-split enumeration:
  * χ = 54Q and Im I = Q on both 256 representatives;
  * the tagged multisets of the 6912 pair are conjugate and chiral;
  * χ = ±18 on that pair.
* The run took 46 minutes.
* Why the proportionality constant is exactly 54 = 2·27 is not derived here. Plausibly it is the ratio of the triple
  count to the Gauss-sum normalisation, but that is not proved.
