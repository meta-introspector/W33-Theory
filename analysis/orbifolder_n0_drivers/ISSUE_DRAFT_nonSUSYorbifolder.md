# DRAFT (not posted) — issue for github.com/StringsIFUNAM/nonSUSYorbifolder

**Title:** Order of a SUSY-breaking twist is wrong when Σvᵢ < 0 (empty spectrum for an equivalent twist vector)

**Body:**

Thank you for the non-SUSY orbifolder: it made an independent cross-check possible for us. With a patched
orbifolder 1.2.1 our N = 0 engine reproduces the massless spectra of all 19 Z₂W × Z_N sample models, gauge group
and field by field. We also found that your `CSector::SortByEigenvalue` already carries the right-moving
oscillator sign correction that we had to add independently.

**The issue.** `CTwistVector::UpdateData()` (ctwist.cpp) treats a twist as SUSY-breaking only when
`sum > prec` ("assuming sum>0"). A twist with Σvᵢ < 0 falls through to the geometric-order branch. Such a
twist is the same group element written with a different lattice representative, for example the Witten twist
(0, −1, −1, −1). In that branch its order comes out as the geometric order n (n·v integer) instead of the order on
spinors (2n when Σ n vᵢ is odd).

**Reproducer.**
1. Copy `Geometry/Geometry_Z3_1_1.txt` and change the first twist line ` 0/1 1 1 1` to ` 0/1 -1 -1 -1`.
2. Point `Models/ZN_models/modelZ3_1_1.txt` at the copy.
3. Result:
   * original: gauge group SO(10)×SU(3)×SO(16)×U(1), 438 fields;
   * modified: "local condition 1 (V_loc² − v_loc²) = −1.00 = 0 mod 2 failed", **empty spectrum**.

**Suggested fix.** Make the test sign-independent: take the smallest n with n·v integral, then set the order to
2n if Σᵢ n vᵢ is odd, and n otherwise. With this rule a patched orbifolder 1.2.1 gives byte-identical spectra
for (0,1,1,1) and (0,−1,−1,−1), and is unchanged on all SUSY twists (Σvᵢ = 0).

**A related limitation, not a bug.** `chalfstate.cpp` hard-codes the right-mover lattice class
(`From_SO8S_Lattice`, "to input only the cospinor lattice ..."). E8×E8 orbifolds whose point-group generator
itself breaks supersymmetry therefore cannot be loaded. An example is the Z6-I twist (1/6, 1/6, 2/3) of Font and
Hernández, hep-th/0202057. The load fails with "State is not invariant under constructing element". An option to
select the right-mover lattice from the twist would cover that class.

Environment: Ubuntu 24.04 (WSL2), g++ 13.3, compiled without readline and driven through the C++ API.
