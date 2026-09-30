# Pass 11202 — where an E6 reflection is local: a double-six and GQ(2,2)

Producer: `analysis/w33_pass11202_reflection_double_six.py`
Certificate: `data/w33_pass11202_reflection_double_six.json`
Regression: `tests/test_w33_pass11202_reflection_double_six.py`

**Background.** Pass 11192 found that among the time reversals of the two-qutrit substrate (the anti-symplectic half
of W(E6)), exactly the 36 reflections have the maximal arrow A = 4. Each reflection is nevertheless block-monomial in
15 splits.

**Result (all 36 reflections).** The dictionary used is 45 splits = tritangent planes, 27 frames = lines, and
line ∈ plane ⇔ split ∈ frame.
* The anti-symplectic involutions fall into two classes: 36 move 12 frames (the reflections) and 540 move 24.
* A reflection's **12 moved frames form a double-six**. They are two sixes of pairwise skew lines, each meeting all
  lines of the other six except its partner, and the reflection swaps partners.
* The 36 reflections give the **36 distinct double-sixes**, matching the classical roots ↔ double-sixes of E6.
* Its **15 fixed splits are exactly the tritangent planes made of three fixed lines**. The 15 fixed frames and 15
  fixed splits form the **generalized quadrangle GQ(2,2) = W(2)**; the GQ axioms were checked.
* In all 15 fixed splits the reflection **swaps the two qutrits** (both diagonal blocks vanish).

**Reading.** Each E6 root is a time reversal that acts as a pure SWAP of the two qutrits in the 15 splits of the W(2)
complementary to its double-six, and exports everything in every split. The line/plane geometry is classical
(Schläfli; Coxeter). The time-reversal and swap reading on the qutrit substrate is the new content.

**Cross-track (added with Passes 11213-11216).** The parallel session reached the same 36-reflection / 15-split
statement independently (Pass 11210, merged), adds that the reflections are exactly the Kramers time reversals
(theta^2 = -1, possible only for an even number of qutrits), and points to the classical prior art already in the
repo: `analysis/PASS7217_7232_double_six_doily_spread_code.md` (the all-c_ij syntheme slice of a double-six).
