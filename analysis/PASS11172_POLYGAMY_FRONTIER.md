# Pass 11172 — spatial polygamy is a thin bulge, and one closed-form family carries its whole frontier

Producer: `analysis/w33_pass11172_polygamy_frontier.py`
Regression: `tests/test_w33_pass11172_polygamy_frontier.py`

**Question.** What is the whole achievable region of (N_AB, N_AC) over three-qutrit pure states, not just its
symmetric point 4/√15?

**Method.** For each weight λ, compare two maxima of N_AB + λ·N_AC:
* the maximum over the closed-form family of Pass 11167 (α|000⟩ + β Σ|aa0⟩ + γ Σ|a0a⟩);
* the global maximum from a weighted dual see-saw (monotone), with 150 random starts per λ.

**Results.**
* Family = global to < 10⁻¹² at every λ tested (0, 0.5, 0.9, 0.93, …, 1).
* **For λ ≤ λ_c = 0.92143, the maximum is exactly 1**, at the corner (1, 0): a maximally entangled AB pair with C
  decoupled. There is no polygamy at all there.
* For λ_c < λ ≤ 1 the maximiser moves along the family arc. It starts at (0.7206, 0.3033) at λ_c and reaches the
  symmetric point (2/√15, 2/√15) at λ = 1, then mirrors beyond.

**Shape.** The upper boundary of the convex hull of the polygamy region has three parts:
* the segment from (1, 0) to the arc endpoint;
* the family arc;
* the mirror segment to (0, 1).

Space allows a bulge beyond the monogamous line N_AB + N_AC = 1 of at most 0.0328, and only for nearly balanced
sharing. In time the same sum is 2 (Pass 11148).

**Scope.** Exact on the family. The claim that the family carries the global boundary is numerical: every weight
checked, with both dual and primal searches.
