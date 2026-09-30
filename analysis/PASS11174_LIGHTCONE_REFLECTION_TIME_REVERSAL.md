# Pass 11174 — the perfect tick's forced pattern is the time-reversing swap of GL(2,3); one oblique Hesse clock as a potential suffices

Producer: `analysis/w33_pass11174_lightcone_reflection_time_reversal.py`
Regression: `tests/test_w33_pass11174_lightcone_reflection_time_reversal.py`

**The link.**
* The 27 events are Sym²(F₃²), via v = (x², xy, y²) (Pass 10946).
* g ∈ GL(2,3) acts by Sym²(g), which preserves q and has det Sym²(g) = det g.
* The 48 elements give 24 Lorentz maps (PGL(2,3) = S₄): 12 of det +1 and 12 of det −1. The other 24 elements of O(q) are
  their negatives.
* **The only coordinate permutations among them are the identity and the light-cone reflection v₀ ↔ v₂**, and the
  reflection is Sym² of the swap [[0,1],[1,0]], whose determinant is −1.

**So the forced pattern of Pass 11168 is time reversal.** All 108 perfect kick–tick–kick steps share π = the light-cone
reflection. That is the Sym² image of an element of the det −1 coset, the coset the repository's determinant character
marks as time-reversing. Pass 10955 identifies that character as three things at once: the extended-Clifford
unitary/antiunitary character, the clock CP grading, and the D₄ half-spin sheet swap.

**Hesse-clock potentials.** The four Hesse clocks are the four projective null directions:
* the light-cone axes n₁ = (1,0,0) and n₂ = (0,0,1), from rays (1,0) and (0,1);
* the oblique directions n₃ = (1,1,1) and n₄ = (1,2,1), from rays (1,1) and (2,1).

Among the 81 potentials N = Σ cᵢ nᵢnᵢᵀ, **exactly 36 make V(N)·K·V(N) perfect**:
* c₁ and c₂ are irrelevant: they are local phase gates on the light-cone qutrits, sitting at both ends of the step.
* The condition is on the oblique clocks alone: c₃ ≠ c₄ and c₃ + c₄ ≠ 2. This follows from the Pass 11168
  determinants, with N₀₁ = N₁₂ = c₃ − c₄ and N₀₂ = c₃ + c₄.
* **A single oblique clock is enough**, e.g. N = (x₀ + x₁ + x₂)². A single light-cone clock is not.

**Reading.** The paper's own clocks supply the interaction that makes the tick a perfect scrambler. The only one that
works is a clock reading along an oblique null direction. The scrambler then time-reverses along exactly the element
of GL(2,3) that the clock group already calls time reversal.

**Related corpus files (read).** Passes 10951 and 10954 record that the det −1 elements of the clock group flip
orientation parity and swap the bipartition classes (the order-8 generator [[0,1],[1,1]] has det −1). That is consistent
with the reading here. None of them identifies the light-cone reflection or the perfect-tick pattern.
