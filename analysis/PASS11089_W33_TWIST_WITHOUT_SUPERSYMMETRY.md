# Pass 11089 — the W(3,3) twist without supersymmetry

Producer: `analysis/w33_pass11089_w33_twist_without_supersymmetry.py`
Certificate: `data/w33_pass11089_w33_twist_without_supersymmetry.json`
Regression: `tests/test_w33_pass11089_w33_twist_without_supersymmetry.py`

The heterotic route with the W(3,3) twist closed on supersymmetric obstructions. In Z6, matter parity dies
at the Fayet–Iliopoulos term, and dimension-4 proton decay operators return. In Z3×Z3, fractional charges
appear. This pass asks what survives, and what opens, without supersymmetry.

Prior art (added in Pass 11092):
* Font and Hernández, *Non-supersymmetric orbifolds*, hep-th/0202057, list the twin twist (⅙, ⅙, ⅔).
* Non-supersymmetric heterotic orbifolds of the SO(16)×SO(16) string have a dedicated tool,
  arXiv:2504.20137, and landscape scans, arXiv:1407.6362 and arXiv:2105.03460.

## 1. What does not depend on supersymmetry

The two Z3×Z3 obstructions of Passes 11024 and 11090 involve only massless **fermions** and **left-moving**
quantum numbers.
* **The index.** A Dirac mass needs a conjugate pair in SM × hidden representations. So the chiral excess
  of the hidden-refined index stays massless in every vacuum that preserves the hidden gauge group.
  Supersymmetry breaking by F-terms, by gaugino condensation, or by anything else that keeps the gauge
  group does not change this.
* **The charge character.** The charge class c = 3Q + t is fixed by the gauge-lattice shift of each fixed
  point. It is therefore the same for every right-moving completion of the same shifts.

**Breaking supersymmetry does not rescue the Z3×Z3 vacua.**

## 2. The prime W(3,3) orbifolds are intrinsically supersymmetric

A non-supersymmetric twin attaches (−1)^F to a generator: v → v′ = v + n, where n is an integer vector
of odd sum (a 2π rotation, which is −1 on spinors). Two parity facts settle the prime orders:
* **Every order-3 shift is even.** For every order-3 shift of E8 (3V in the lattice), **9V² is even**.
  In the integral class 9V² = Σuᵢ² with Σuᵢ even. In the half-integral class
  9V² = Σuᵢ(uᵢ+1) + 2.
* **Every order-3 twin twist is odd.** For every order-3 twist with (−1)^F attached,
  **9v′² = 2 + 6(n₁ − n₃) + 9n² is odd**.

Modular invariance, 3(V² − v′²) ∈ 2Z, needs the two parities to agree, and they never do. So **no
non-supersymmetric T⁶/Z3 or T⁶/Z3×Z3 of the E8×E8 string exists with an order-3 shift**, the W(3,3) (A8 Kac)
shift included.

> **Correction (Pass 11092).** An earlier version said "for any shifts". That is an over-read. The proof
> assumes 3V ∈ Λ.
> * The twisted generator has order 6 on fermions. With a shift of order 6, the Z3 geometry does carry
>   non-supersymmetric models: Font–Hernández (hep-th/0202057) treat the Z6 action v = (⅓, ⅓, ⅓) on the
>   T⁶/Z3 lattice.
> * The SO(16)×SO(16) string has been orbifolded on all 138 Abelian geometries, Z3×Z3 included
>   (arXiv:2105.03460).
>
> What stands is the statement about **twins**: the same order-3 shift, with (−1)^F attached.

This is checked three ways:
* by the parity proof above;
* by enumerating order-3 shifts (29.9 million during development; a smaller bound in the producer);
* by orbifolder. It rejects the twin of the W(3,3) Z3 model: "3 (V_M² − v_M²) = 2.33 = 0 mod 2 failed". The
  supersymmetric control loads with N = 1 and gauge group SU(9)×SO(14)×U(1).

## 3. Which families have non-supersymmetric twins

Adding (−1)^F changes every sign combination ±v₁±v₂±v₃ by an odd integer. A supercharge survives iff
some combination is an even integer. So a twin can break supersymmetry only if **no sign combination of
the original twist is an odd integer**. Computed from the orbifolder geometry twists:

| family | non-SUSY twin | same shifts and Wilson lines |
|---|---|---|
| Z3, Z3×Z3 | as a twist only | **never** (§2: no shifts at all) |
| **Z6-II** | **none**: ⅙ + ⅓ − ½ = 0, supersymmetric again | — |
| **Z6-I** | yes, (⅙, ⅙, ⅔) | **yes** |
| Z12-I, Z2×Z6-I, Z2×Z6-II, Z3×Z6, Z6×Z6 | yes | yes |

**All 87 Z6-I W(3,3) Standard Models have non-supersymmetric twins.** Each keeps the same shifts and
Wilson lines. orbifolder loads all 87 twins with **N = 0**, passing modular invariance and the Wilson-line
conditions; the originals load with N = 1.

## 4. The open door, and what blocks it

orbifolder 1.2 cannot compute an N = 0 spectrum: its spectrum code groups states into supermultiplets,
and it segfaults on the twins. So the massless and tachyonic spectra of the 87 twins need a
**non-supersymmetric spectrum engine**.

The failure has been localized, and two of the three faults are fixed in a patched build (`~/orb/buildns`).
With both fixes the supersymmetric control output is **bit-identical**, so the N = 1 path is untouched.

| fault | location | status |
|---|---|---|
| no N = 0 branch in the multiplet classification | `CState::FindSUSYMultiplets` | **fixed**: each right-moving state classified by helicity q₀ (−½ left-handed Weyl fermion, +½ right-handed, ∓1 vector, ∓2 graviton, 0 scalar) |
| `RecursiveCounting(..., MaxDigits[0], ...)` on an empty vector | `CFixedBrane::FindSUSYMultiplets` | **fixed**: N = 0 has one, empty, supercharge combination |
| "State is not invariant under constructing element" | `CFixedBrane::CreateStates` | **fixed in Pass 11092**. The diagnosis given here ("missing vacuum phase") was wrong: the right-mover eigenvalue counts an oscillator's number operator with the wrong sign (−N instead of +N), which only matters when a massless right-mover carries an oscillator, never in SUSY models |

Pass 11092 fixes this fault and two more (a modulus-labelling segfault, and the anomalous-U(1) generator being
built only for N = 1). It then computes all 87 twins: they are anomaly-free three-generation models, and
**every one is tachyonic** in its θ sector.

This was the most promising route left on the string side. Pass 11092 closes it at the orbifold point. With no superpartners there are no
dimension-4 or dimension-5 proton-decay operators, which is the obstruction that closed the
supersymmetric Z6-I route. The questions to settle, in order:
1. tachyons in the odd twisted sectors (**answered in Pass 11092: present in all 87**, θ and θ⁵ sectors only);
2. whether the chiral content is still three Standard-Model generations (**answered in Pass 11092: yes, 87/87**);
3. the one-loop tadpole and cosmological constant (moot at a tachyonic point).

Scope. "Twin" means (−1)^F attached to one generator with the same shifts. Other non-supersymmetric
constructions, for example a different lattice from the start, are not classified here. The
SO(16)×SO(16) string cannot host the twist as an E8 automorphism: SU(9) ⊄ SO(16).
