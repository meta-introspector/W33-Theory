# Pass 11310 — the U81 cocycle orientation is not the cubic gate's orientation (a tempting identification ruled out)

Producer: `analysis/w33_pass11310_magic_orientation_cocycle.py`; GAP: `analysis/gap/w33_pass11310_magic_sylow.g`,
`analysis/gap/w33_pass11310_u81_sylow.g`
Certificate: `data/w33_pass11310_magic_orientation_cocycle.json`
Regression: `tests/test_w33_pass11308_11312.py`

**Background (Codex).**
* Passes 11067–11069 and `w33_20261001_u81_antiunitary_orientation.py` write the chamber group as
  (h,d)·_s(h',D) = (hh', d + D + s·κ(h,h')), with κ = A²b − 2Ac.
* s = ±1 gives U81, and qutrit complex conjugation exchanges s = +1 and s = −1.
* "An energetic coupling selecting one sign remains open." The corpus already identifies U81 as the chamber Sylow-3
  subgroup of Sp(4,3) (Passes 11061–11075).

**The tempting identification.** Complex conjugation also exchanges the cubic gate T with T⁻¹ = T*. So is the cocycle
orientation the magic gate's choice of cube root of Z (T³ = Z versus (T⁻¹)³ = Z⁻¹)?

**Answer: no.** The one-qutrit "magic Sylow" group ⟨X, S, T⟩ modulo scalars was computed exactly, with elements
X^b·diag(ζ^{f(x)}) and f taken modulo constants:

| | ⟨X,S,T⟩/scalars | U81 (s = ±1) |
|---|---|---|
| order, centre, derived, class | 81, 3, 9, 3 | 81, 3, 9, 3 |
| elements of order 9 | **18** | **36** |
| GAP SmallGroup id | **(81, 9)** | **(81, 7) = C₃≀C₃ = Syl₃ Sp(4,3) = Syl₃ PSp(4,3)** |

* The groups are not isomorphic.
* No assignment e₀ ↦ g₀, e₁ ↦ X extends to a homomorphism, for either s.
* T and T⁻¹ generate the same magic group.

**Consequence.** The cocycle orientation lives at the Clifford (symplectic) level of the substrate, not in the
one-qutrit magic layer. Selecting it is not the same as choosing T over T†. Whatever selects s must be found among
symplectic or chamber data; this pass removes one candidate.
