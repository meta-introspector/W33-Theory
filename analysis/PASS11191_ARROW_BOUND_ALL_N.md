# Pass 11191 — the arrow bound A ≤ n for every number of qutrits whenever S + S⁻¹ is semisimple (two thirds of every group), and at five and six qutrits by search

Producer: `analysis/w33_pass11191_arrow_five_six_qutrits.py`
Frozen: `data/w33_pass11191_arrow_five_six_qutrits.json`
Regression: `tests/test_w33_pass11191_arrow_five_six_qutrits.py`

**Re-scoping.** Pass 11191 was reserved (4065cb3b) for "polygamy beyond the torus class". A first attempt paired the
one-cut bounds of Pass 11167 across parties: N_AB ≤ min(h(p_C), g(λ_B)) and N_AC ≤ min(h(p_B), g(λ_C)). It only lowers
the rigorous global bound from 25/18 = 1.3889 to about 1.3827. Closing the gap to 4/√15 needs a joint (DPS-type)
relaxation, so that item stays open. The pass number is used instead for the arrow-of-time bound of Pass 11188.

**Setting (Pass 11188).** A(S) = 2n − max Σ dim(P ∩ SP), taken over orthogonal families of nondegenerate planes.
* A ≤ n − 1 ⇒ S fixes a nondegenerate plane.
* A fixed plane P gives A(S) ≤ A(S|P^⊥).
* Exact for n ≤ 4: A ≤ n, with equality iff S fixes no nondegenerate plane.

**Lemma (bi-Lagrangian criterion).** Let T = S + S⁻¹, which is ω-self-adjoint, so ω_T(x, y) = ω(x, Ty) is alternating.
Suppose L is a Lagrangian for both ω and ω_T, with L ∩ SL = 0. Then S has a *half-moving split*, and so A(S) ≤ n.
*Proof.*
* On L, b(x, y) = ω(x, Sy) satisfies b(x, y) − b(y, x) = ω_T(x, y) = 0, so b is symmetric.
* b is nondegenerate: b(x, L) = 0 means x ∈ (SL)^⊥ = SL, and L ∩ SL = 0.
* Over F₃, b therefore has an orthogonal basis eᵢ with b(eᵢ, eᵢ) ≠ 0.
* The planes Pᵢ = ⟨eᵢ, Seᵢ⟩ are nondegenerate, since ω(eᵢ, Seᵢ) = b(eᵢ, eᵢ).
* They are mutually orthogonal, since ω(eᵢ, eⱼ), ω(eᵢ, Seⱼ), ω(Seᵢ, eⱼ) and ω(Seᵢ, Seⱼ) all vanish for i ≠ j.
* Each Pᵢ meets its image in Seᵢ. ∎

Every half-moving split found below by search is of this form (checked: `all_bilagrangian`).

**Theorem (all n).** If T = S + S⁻¹ is semisimple, then A(S) ≤ n, with equality iff S fixes no nondegenerate plane.
In particular this holds for every tick of order prime to 3, since semisimple S gives semisimple T.

*Proof.* Induct on n. If S fixes a nondegenerate plane P, restrict to P^⊥; T stays semisimple there. Otherwise, T is
self-adjoint, so its primary components V_f (one for each irreducible factor f of its minimal polynomial) are
ω-orthogonal, S-invariant and nondegenerate. It is enough to build a bi-Lagrangian L_f transverse to S L_f in each.

On V_f, T acts as a scalar t in F = F₃[x]/(f). Write ω = Tr_{F/F₃} ∘ h with h F-bilinear and alternating; S preserves h,
and S + S⁻¹ = t. Any h-Lagrangian F-subspace is Lagrangian for both ω and ω_T = Tr(t·h). There are three cases:
* **x² − tx + 1 irreducible over F.** V_f is a space over E = F[S] ≅ F_{q²}, and S = μ with μ⁻¹ = μ^q. Then
  h(ax, y) = h(x, a^q y). Define the sesquilinear form H by Tr_{E/F}(c·H(x, y)) = h(cx, y). It is skew-Hermitian, and κH
  is Hermitian for any κ with κ^q = −κ. Over a finite field a Hermitian form has an orthonormal basis eᵢ. Take L = F-span
  of the eᵢ: it is h-isotropic because Tr(κ⁻¹) = 0, and it spans V_f over E, so L ∩ μL = 0.
* **Distinct roots μ ≠ μ⁻¹ in F.** V_f = V_μ ⊕ V_{μ⁻¹}, both h-isotropic and paired by h. Take L = {x + φx}, with φ
  sending a basis of V_μ to its h-dual basis. L is isotropic because h(x, φy) is symmetric. And L ∩ SL = 0 because SL is
  the graph of μ⁻²φ and μ² ≠ 1.
* **Double root μ = ±1.** Let N = S − μ, so N² = 0 and N is ω-skew. N ≠ 0, since V_f contains nondegenerate planes. If
  rank N < ½ dim V_f, then ker N is not isotropic and contains an S-invariant nondegenerate plane, which is excluded.
  So ker N = im N is Lagrangian. Any Lagrangian L transverse to it has L ∩ SL = 0, because Nx ∈ L forces Nx ∈ L ∩ ker N = 0.
  And ω_T = tω. ∎

**Coverage (exact, from the GAP class lists of Pass 11188).**

| n | fraction of PSp(2n, 3) with T semisimple | with S semisimple |
|---|---|---|
| 2 | 18000/25920 = 25/36 | 9280/25920 |
| 3 | 3137097600/4585351680 ≈ 0.684 | ≈ 0.406 |
| 4 | ≈ 0.679 | ≈ 0.367 |

So the theorem settles A ≤ n, and A = n exactly when there is no invariant plane, for about two thirds of every group,
at every n.

**Five and six qutrits (search).** Planes span(x, Sx) with ω(x, Sx) ≠ 0 meet their image. A depth-first search for n
mutually orthogonal such planes, in seeded random order with restarts, gives:
* **n = 5:** 301 ticks, namely 300 random and the regular unipotent (one Jordan block). 136 of them fix no nondegenerate
  plane. **Every one of those 136 has a verified half-moving split**, including 27 with irreducible characteristic
  polynomial and the regular unipotent. The other 165 fix a plane, so A ≤ 4 by the exact n = 4 theorem.
* **n = 6:** 41 ticks. All 22 with no invariant plane have a verified half-moving split, including the regular unipotent.
  Aside: a fixed point order is pathological for this tick, and random order finds a split at once.

**Status of the conjecture A ≤ n.**
* Proved for all n when S + S⁻¹ is semisimple.
* Proved exactly for n ≤ 4.
* Supported by every sampled tick at n = 5, 6.
* Open: ticks whose S + S⁻¹ is not semisimple (for example, regular unipotents), for n ≥ 5. The lemma reduces it to the
  existence of a bi-Lagrangian of the pencil (ω, ω_T) transverse to its S-image, a pure linear-algebra question about
  self-adjoint operators with nilpotent part.
