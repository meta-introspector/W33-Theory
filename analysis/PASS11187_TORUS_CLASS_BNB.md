# Pass 11187 — 4/√15 is the maximum of N_AB + N_AC on the torus-symmetric class: a certified branch and bound

Producer: `analysis/w33_pass11187_torus_class_bnb.py`
Regression: `tests/test_w33_pass11187_torus_class_bnb.py`

**The class.** The states fixed by the maximal torus of the optimum's U(2) symmetry:
ψ = α|000⟩ + Σ_{a=1,2} (βₐ|aa0⟩ + γₐ|a0a⟩), with phases removable. With sₐ = βₐ², tₐ = γₐ², c = α²:

    N_AB = Σₐ (√(tₐ² + 4c·sₐ) − tₐ)/2 + √(s₁s₂),     N_AC = Σₐ (√(sₐ² + 4c·tₐ) − sₐ)/2 + √(t₁t₂).

This is checked against direct negativity. Pass 11167 settled the 2-parameter slice s₁ = s₂, t₁ = t₂; this pass settles
the full 4-parameter class.

**Proof.**
1. **AM–GM.** f ≤ Ψ = ½ Σₐ φ(sₐ, tₐ; c), with φ(s,t;c) = √(t² + 4cs) + √(s² + 4ct). Equality holds on the symmetric
   slice.
2. **Core.** Suppose φ(·,·;c) is concave, for every c in the box's range, on the rectangle hull of both (sₐ, tₐ). The test
   splits that hull into 6×6×6 cells, each checked with interval Hessian bounds. Then Jensen gives
   Ψ ≤ Ψ(sym) = f(sym) ≤ 4/√15, the last step by Pass 11167.
3. **Elsewhere.** A box is closed when an upper bound falls below 4/√15. The bound is the minimum of the monotone corner
   bound and the first-order centred form with an interval gradient. The centred form is used only on boxes strictly
   inside the simplex.

**Result.** The branch and bound **terminates**: 265 725 boxes processed, 131 374 closed by bounds, 1450 by the
concavity core, 39 empty, in about 16 s.

**Controls.**
* The analytic Hessian matches finite differences (4·10⁻⁶).
* The interval test rejects a non-concave cell (s = 0.4, t = 0.05, c = 0.3).
* It accepts the optimum's neighbourhood.
* 400 unconstrained optimisations on the class all converge to 4/√15 at sₐ = tₐ = 2/15.

**Scope.**
* Floating-point bounds with a safety margin of 10⁻⁹ (not directed-rounding interval arithmetic).
* A restricted class (torus-invariant states), larger than Pass 11167's.
* The global problem remains open. The evidence for it is 3000/3000 dual see-saws (Pass 11167) and the frontier match
  (Pass 11172).
