# Pass 11024 — the Z3×Z3 parity vacua are supersymmetric, and carry massless fractional charges

Producer: `analysis/w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.py`
Exact solver (new, shared): `analysis/w33_exact_monomial_orders.py`
Certificate: `data/w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.json`
Regression: `tests/test_w33_pass11024_z3xz3_susy_vacua_fractional_charge_index.py`

Pass 10979 left one door open on the string route. There were 28 T6/Z3×Z3 vacua with these properties:
* matter parity survives the Fayet–Iliopoulos (FI) term, under established selection rules;
* the exotic triplets are massive;
* D-flatness holds with every field nonzero.

Their superpotential has linearly independent exponent vectors through degree 11. The paper said a
supersymmetric vacuum "would need a twelfth-order coupling to outweigh a cubic one". This pass settles those
vacua. It also sharpens that sentence.

## A. The superpotential on each vacuum is Z3-graded, to all orders

On every one of the 28 vacua S:
* **The invariants form a free monoid.** The monomials invariant under all continuous charges form a free
  monoid on k composites p_i = φ^{r_i}. k = 1 in 9 vacua and k = 2 in 19. The composites are products of
  3 or 4 fields with disjoint supports, and the lattice index is 1.
* **W sees exactly one residue class.** A monomial p^a satisfies the space-group and R selection rules iff
  **a₁ + … + a_k ≡ 1 (mod 3)**. The rule is checked on a full period, so it holds to all orders.

Hence **W|_S = F(p) with F(ωp) = ω F(p)**, ω³ = 1: a Z3 R-symmetry acting on the composites. With every
field nonzero, F_φ = 0 ⟺ ∇_p F = 0.

All exponent vectors lie in a k-dimensional span, so the first linear dependence is the first monomial
beyond the k linear ones, at degree 4·min(deg p). That is Pass 10979's "degree 12" for cubic composites.
For the quartic composites the prediction is degree 16, beyond 10979's enumeration to degree 14. The
predicted first dependent degree matches 10979 in **28 of 28** vacua.

## B. Supersymmetric vacua exist (correcting a framing in Pass 10979)

∇F = 0 is a system in one or two variables: linear plus cubic terms. It has the Bézout number of roots:
3 for k = 1 and about 9 for k = 2. The producer draws random complex O(1) couplings, solves ∇F = 0, and
lifts every root to fields. The lift uses log|φ| from the composite values and the D-terms with the FI
term, ξ = g² Tr Q₀ / 192π² (Holotrade cfdc1f2 convention, g² = 1/2). Every lifted point is verified to
F = D = 0 with every field nonzero. **Verified SUSY roots exist in 28 of 28 vacua.**

The VEV is set by the coupling ratio, |p|³ ~ |c₁/c₄|, not by ξ. With O(1) couplings the largest field VEV
is **0.80–1.07 M_s** (median over roots, per vacuum).

So the Pass 10979 sentence should read: *a supersymmetric vacuum exists for generic couplings, and a
twelfth-order coupling larger than the cubic one only lowers its scale.* Its exact statement (no balance
below degree 12) stands. Its numerical search found nothing because it started at the FI scale and never
reached the region near M_s.

## C. Massless fractional charges: a hidden-representation index

**A mass term needs a conjugate pair in SM × hidden representations.** So the excess
#X − #conj(X) in each class stays massless in any vacuum that leaves the hidden gauge group unbroken.
Charge counting shows the whole spectrum is vector-like with respect to the SM alone, 87 against 87. Not
with respect to SM × hidden:

| model | chiral fractional states | of which hidden singlets | b(SU(4)) | b(SU(3)') |
|---|---|---|---|---|
| c3 (Z3Z3_0001_c3_739) | 48 (charge ±1/3) | **6** | −7 | −16 |
| c4 (Z3Z3_0001_c4_3729) | 48 | **6** | −4 | −13 |

The hidden-charged part consists of SU(4) 4, 4̄ and 6 multiplets. All 28 parity vacua contain only
hidden-singlet fields, so **every one carries six exactly massless colour-singlet fermions of electric
charge 1/3.**

Two ways out fail:
* **Confinement.** Both hidden groups are infrared free (b < 0), so they do not confine.
* **Hidden hadrons.** Even if SU(4) confined, every SU(4) invariant carries charge in Z/3, so the lightest
  fractional state would be stable.

Context: fractionally charged states are the generic price of a hypercharge that is not embedded in a
GUT at level one (Schellekens, "Electric charge quantization in string theory", Phys. Lett. B 237 (1990)
363). What this pass adds is the refined, hidden-representation index for the W(3,3)-twisted vacua. It
decides exactly when such states must stay massless, and which hidden factor must break to avoid that.

Control: the same index vanishes for c2_4158 and for most Z6 models (see the certificate). It is a
property of these spectra, not of the method.

## D. Hidden condensates: the index can be cured, the spectrum cannot

**Which condensates could possibly help.** Break each hidden factor in turn and recompute the index:
* **c3, c4:** the fractional excess vanishes **iff the hidden SU(4) is broken**. Breaking SU(3)′ alone
  leaves all 48 states.
* **c1:** it vanishes iff the first hidden SU(2) is broken.

So a vacuum can avoid the index only if it contains a condensate of that *curing* factor.

**The census.** By Luty–Taylor, D-flat directions correspond to holomorphic invariants. Every quadratic
hidden invariant is added to the exact parity × FI census as a composite singlet:
* SU(N) mesons f·f̄, SU(2) doublet pairs, and 6·6 of SU(4);
* each carries summed U(1), space-group and R charges;
* the parity may include the broken Cartan, so only the composite need be even.

Every FI-cancelling extreme ray **containing a curing condensate** is searched for realizability, depth
first; realizability is monotone. That makes the search complete for this class. Baryonic invariants
(585 in c3) are left as named residue: they make the discrete-choice search explode.

| model | FI rays | rays with a curing condensate | realizable | MSSM-viable parity | parity aligned with the broken Cartan |
|---|---|---|---|---|---|
| c1_2822 | 68 | 4 | 4 | 1 | **0** |
| c3_739 | 105 | 39 | 39 | 39 | **39** |
| c4_3729 | 117 | 105 | 15 | 15 | **15** |
| c2_4158 | Farkas-closed with all composites | | | | |

**The candidates.** For each viable vacuum, every SM-charged class is analysed after the breaking:
* components are expanded: SU(4) → SU(3) gives 4 → 3+1 and 6 → 3+3̄; SU(2) doublets split;
* centre charges (N-ality) of every hidden SU(N) enter as extra discrete charges;
* mass ranks come from exact orders on the vacuum fields.

Treating any N-ality-allowed contraction as a coupling can only over-count masses, so the light counts
are lower bounds.

* **c1 (a near miss).** Breaking the first hidden SU(2) zeroes the index *exactly*. The broken factor's
  doublet classes mirror the singlet excess count for count: 12 = 6×2 and 6 = 3×2. But two condensates of
  one SU(2) cannot share a parity: P(n₂₀) ≡ P(n₃₆) or P(n₂₀) ≡ P(n₃₉) admits no MSSM-viable parity. Even
  ignoring that, 158 fractional multiplets stay massless at all orders. The vacuum also leaves 5 unbroken
  U(1)s, i.e. 4 extra massless U(1)′ gauge bosons.
* **c3 (39 genuine parity vacua with SU(4) broken).** Every one is supersymmetric (W ≡ 0 on the vacuum) and
  has the right MSSM content: 3 quark doublets and one light Higgs pair. In the best one, a single meson
  n₂·n₄₀ plus five singlets, the hidden-singlet charge-1/3 class is now balanced, 42 against 42, so the
  index is cured. But the couplings give rank 19, not 42. **Every one of the 39 keeps 128–166 fractionally
  charged multiplets massless at all orders.**
* **c4 (15 more).** The same picture: all supersymmetric, parity aligned, and 120–132 fractional
  multiplets massless in each.
* **Why: an unbroken discrete symmetry.** Classify each absent mass term. Either no integer exponents
  exist even allowing negative ones, so an exact symmetry of the vacuum forbids it, or only non-negative
  ones fail, a holomorphy obstruction. In the best c3 vacuum, 1475 of 1512 hidden-singlet charge-1/3
  entries are forbidden at the lattice level. The vacuum preserves U(1)_Y × **Z₃⁴ × Z₉**, the torsion of
  the charge lattice modulo the vacuum's, of order 3⁶ = 729. That symmetry protects **120 of the 128**
  massless fractional states: to all orders, perturbative or not, and against Kähler (Giudice–Masiero)
  terms. In c1 it protects 146 of 158. Across all 54 SU(4)-condensed vacua only U(1)_Y survives, and the
  unbroken discrete group is always a **3-group**:
  * Z₃²×Z₉², of order 3⁶, in 36 vacua;
  * Z₃⁴×Z₉, of order 3⁶, in 12;
  * Z₃²×Z₉, of order 3⁴, in 6.

  It protects 106–154 of the massless fractional states in c3 and 120–124 in c4.
* **Single-meson rescue of the 28 singlet vacua** (adding one of the 70 SU(4) mesons to each): no meson
  pairs the 42 states; the best rank is 33. Without a meson, at least 14 (c3) and 16 (c4) states stay
  massless.

## E. An independent certificate for the integer programs of Passes 10978–10980

Those passes computed monomial orders with mixed-integer programs, putting the Z_N congruences in
integer slack variables. `w33_exact_monomial_orders.py` replaces tolerance with exact arithmetic:
* it solves the charge system by Smith normal form;
* it LLL-reduces the kernel and recentres the particular solution in exact integers;
* it finishes with a congruence-free integer program in the 1–3 kernel dimensions;
* it re-verifies every solution in rational arithmetic.

Results:
* **Pass 10979:** 72 checks, 0 differences.
* **Pass 10978:** 57 checks, 0 differences, including Z12I_1063's linear n_36 at order 4.
* **Pass 10980:** the rerun with exact orders reproduces its summary identically.

Process note. A first version of this audit disagreed with the integer programs and looked like a
correction. The cause was a sign convention: a coupling X·Y·S^e needs Σe·v = w − v(X) − v(Y). The new
solver and its brute-force check shared the flipped sign. Only reproducing a *known positive* of the old
tool exposed it. The API now encodes the convention (`coupling(fields, w)`), and the regression test
includes a known positive.

## Verdict

T6/Z3×Z3 with the W(3,3) twist has matter-parity-preserving, supersymmetric, FI-cancelling vacua. Pass
10979's door was real, and more of them exist than it knew. Whether the hidden sector is left intact or
broken, every one carries massless fractionally charged matter:
* **Hidden sector intact:** six exactly massless charge-1/3 fermions, forced by a hidden-representation
  index. The hidden groups are infrared free and cannot confine them.
* **Hidden SU(4) or SU(2) broken:** there are 55 such parity vacua (1 in c1, 39 in c3, 15 in c4). The
  index is cured, but at least 120 fractional multiplets stay massless in each, most of them protected by
  an exact unbroken discrete symmetry of the vacuum, a 3-group of order 3⁴ or 3⁶.

The string route therefore closes in all nine families built on the W(3,3) twist:
* in eight, by matter parity or massless exotics;
* in Z3×Z3, by fractional charge.

**Scope.** "Closes" means under the established selection rules. It covers:
* every realizable FI direction built from singlets;
* every realizable FI direction containing a quadratic condensate of the factor whose breaking can cure the
  index;
* the single SU(4)-meson additions to the 28 singlet vacua.

Not covered, and so the precise residue of open problem 1: baryonic condensates, invariants mixing two
hidden factors, and non-aligned breaking patterns (e.g. SU(4) → SU(2)).
