# Pass 11168 — no Lorentz-invariant dynamics of the clock is perfect; a kick–tick–kick step can be, and then it always time-reverses along the light-cone reflection

Producer: `analysis/w33_pass11168_clock_dual_pair.py`
Regression: `tests/test_w33_pass11168_clock_dual_pair.py`

**Setup.** The three qutrits are the coordinates (v0, v1, v2) of F₃³, with q = v0v2 − v1². The polar Gram matrix is
M = [[0,0,1],[0,1,0],[1,0,0]], and M = M⁻¹.
* K is the clock tick of Theorem 4.3 (Pass 11166, the kinetic term e^{−iq(p)}): Z_b → Z_b X_{Mb}.
* V(N) is a potential kick e^{−iN(x)}: X_a → X_a Z_{Na}.
* V = V(M) is the Fourier dual of K.
* O(q) is the Lorentz group of the finite light cone: 48 linear maps.

**Results.**
* **⟨K, V⟩ has order 24**, a copy of SL(2,3) built from M, and **no perfect element**. Because M₀₁ = M₁₂ = 0, the four
  blocks between the transverse qutrit and the light-cone pair vanish in every element.
* **⟨K, V, O(q)⟩ has order 576**, since the Lorentz maps commute with K and V (LMLᵀ = M). **It has no perfect element
  either.** The Lorentz-invariant potentials are exactly N = cM, so they stay inside this group.
* **One tick is never enough.** K·V(N) has off-diagonal determinants −M_ij N_ij, which vanish at (0,1) for all 729 N.
* **Tick–kick–tick is never perfect.** The (0,1) determinant of K·V(N)·K is N₂₁N₀₁ − N₂₁N₀₁ = 0 identically.
* **Kick–tick–kick V(N)·K·V(N) is perfect for exactly 108 of the 729 potentials**, all of them Lorentz-breaking. The
  simplest is N(x) = x₁(x₀ + x₂), which couples the transverse coordinate symmetrically to both light-cone coordinates.
* **Every one of the 108 has π = (v0 ↦ v2, v1 ↦ v1, v2 ↦ v0)** (π as in Pass 11169). This is forced:
  * the light-cone diagonal blocks have determinant (1 + N₀₂)², a square, so they can never reverse orientation;
  * the transverse block has determinant 1 + N₀₁N₁₂, which perfection forces to −1.
* Also: V(N₁)·K·V(N₂) is perfect for 23 328 = 2⁵·3⁶ ordered pairs, and K·V(N)·K·V(N) for 30 potentials.

**Reading.**
* The paper's clock and everything built from its own symmetry never scramble perfectly.
* A perfect tick needs an interaction that breaks the finite Lorentz invariance.
* When it scrambles perfectly, each light-cone coordinate's past reappears time-reversed in the *other* light-cone
  coordinate's future, and the transverse coordinate's past reappears time-reversed in its own.
* So π is the light-cone reflection v0 ↔ v2, the only nontrivial coordinate permutation preserving q.

**Scope.** "Perfect" refers to the paper's coordinate split into three qutrits. The obstruction for Lorentz-invariant
dynamics uses M₀₁ = M₁₂ = 0, which holds in these coordinates. The Pass 11166 rank-≤1 theorem for translation-invariant
ticks is basis-free.
