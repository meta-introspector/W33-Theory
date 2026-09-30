# Pass 11201 — the time-reversal twist is a measurable chiral phase: a quadratic Gauss sum in the Choi state

Producer: `analysis/w33_pass11201_twist_chiral_invariant.py`
Certificate: `data/w33_pass11201_twist_chiral_invariant.json`
Regression: `tests/test_w33_pass11199_to_11206_frozen.py`

**Question.** Pass 11189 separated the time-reversal pairs of three-qutrit split relations by finite-field signs: an
SL(2,3) holonomy class for the 6912 pair and a symplectic twist Q = ±1 for the 256 pair. Are these physical, i.e.
visible in local measurements of the gate?

**Invariants.** For the Choi state ρ of the gate, I(σ) = tr[ρ^{⊗3} ∏_p P_{σ_p}], with σ_p ∈ S3 permuting the three
copies on party p.
* These are local-unitary invariant third moments, accessible by randomized measurements.
* Complex conjugation, the time reversal, gives I(σ)^* = I(σ⁻¹).

**Exact evaluation.** For a stabilizer state with Lagrangian L, the value of I(σ) is
3⁻¹⁸ ∏_p t_p Σ_{(v₁,v₂)∈C} ζ^{Σ_{cyclic p} ±ω_p(v₁,v₂)/2}.
* t_p = 27, 9, 3 for an identity, a transposition and a 3-cycle.
* C is cut out by linear constraints: v₁,p = v₂,p = 0 for an identity; the fixed copy vanishes for a transposition.
* The stabiliser's linear phases cancel, so I depends only on L and σ. It is a Gauss sum.
* The evaluator reproduces direct tensor contraction.

**What follows.**
* If all 3-cycles turn the same way, isotropy kills the phase and I is real and positive (30 of 30 random checks).
* I is real whenever an odd relabelling of the copies preserves the pattern of transpositions. Chirality needs 3-cycles
  of both orientations *and* two different transpositions.
* **For such a pattern, on the 256 relation, I = Q · i√3/243 exactly**: purely imaginary, a quadratic Gauss-sum value,
  with the sign of the twist. It is constant on all 256 members of each orbital (−1 on one, +1 on its reverse, in units
  of i√3/243) and changes sign under time reversal.
* The 6912 pair has its own chiral pattern, with value ±i√3/729, which also flips.
* **Caveat (scope).** The sign is carried by labelled registers. Of the 36 relabellings of inputs and outputs, 6 give
  +1, 6 give −1 and 24 give 0, so the fully symmetrised relation invariant is real. The arrow is locally measurable
  only when the registers are labelled compatibly with the pairing the relation itself defines (its invertible
  blocks).

**Reading.** The finite-field twist that tells a relation from its time reverse is the sign of a Gauss sum: an
imaginary third moment ±i√3/243 of the gate's Choi state. The √3 is the qutrit's i√3 = ζ − ζ̄. Time reversal
conjugates it.

**Correction during the pass.** A first version of the exact evaluator had the 3-cycle orientation reversed, which
gave the conjugate value. It was caught by the built-in comparison with direct contraction, which the real-valued
calibration patterns could not detect. The sign relation is **+Q**, as checked.
