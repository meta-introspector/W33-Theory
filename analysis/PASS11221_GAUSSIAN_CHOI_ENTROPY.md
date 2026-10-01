# Pass 11221: the Gaussian arrow is an entanglement rate

Producer: `analysis/w33_pass11221_gaussian_choi_entropy.py`
Certificate: `data/w33_pass11221_gaussian_choi_entropy.json`
Regression: `tests/test_w33_pass11221_gaussian_choi_entropy.py`

**Background.**
* Pass 11219 proved that, for Gaussian (free bosonic) dynamics, A(S) = n − c(S) is the least total number of quadrature
  channels, Σ_k rank S[P_k^⊥ ← P_k], through which modes must leak into one another.
* Its first entropic reading, the vacuum's entanglement, failed. A loxodromic quartet can leave the vacuum a product
  state.
* The right object is the **operator** entanglement of the Gaussian unitary U_S. It is measured on its Choi state with
  squeezing r:
  * each mode is maximally correlated with a reference through a two-mode squeezed vacuum (TMSV);
  * the system modes are then evolved by S;
  * the cut is (mode k + its reference) | (all the rest).

## Result

For every mode of every mode decomposition,

    S_k(r) = 2 · rank S[P_k^⊥ ← P_k] · r + O(1)   (nats),

so

    **A(S) = ½ · min over mode decompositions of lim_{r→∞} d/dr Σ_k S_k(r).**

The arrow of a Gaussian dynamics is half its minimal total operator-entanglement rate per unit squeezing.
* It is 0 exactly for normal-mode-decomposable dynamics (elliptic and hyperbolic).
* Each loxodromic quartet or Jordan defect adds 2: its best split puts one channel on each of its two modes, i.e.
  4 nats per unit squeezing.
* A generic split has rank 2 on every mode, i.e. 4 nats per mode per unit squeezing.

**Why the slope is 2 per channel.**
* Given its reference, mode k's input is known up to e^{−r} noise. The remaining uncertainty of its output is
  V · S_{kR} S_{kR}ᵀ, with V = cosh 2r / 2 the input variance of the other modes and S_{kR} = S[P_k ← P_k^⊥].
  * This block has the same rank as S[P_k^⊥ ← P_k], by symplecticity.
* So det σ_A ≍ V² · det(V S_{kR} S_{kR}ᵀ + O(e^{−2r})). This gives:
  * rank 0: det σ_A ≍ 1, so ν₁, ν₂ → ½ and the slope is 0;
  * rank 1: det σ_A ≍ V², so ν₁ ∝ V and ν₂ = O(1), and the slope is 2;
  * rank 2: det σ_A ≍ V⁴, so ν₁, ν₂ ∝ V, and the slope is 4.
* Each symplectic eigenvalue ν ∝ V = e^{2r}/4 contributes log ν ≈ 2r, the entropy of one TMSV link.

## Checks

* Setup: 18 dynamics from Pass 11219's families (elliptic, hyperbolic and loxodromic mixtures up to six modes), each
  with its constructed optimal split and two random splits. Also the elliptic Jordan pair at a Krein collision, with its
  numerically found split.
* Method: Choi-state entropies of every mode at r = 5 and r = 6. Larger r amplifies float64 rounding of the symplectic
  matrices by e^{4r}; a first run at r = 8–9 showed exactly that, with nonzero slopes on rank-0 modes.
* Results:
  * **slope_k = 2·rank_k on every mode of every split** (to 10⁻²);
  * half the optimal total slope equals A (= 2 × the number of quartets) in every case;
  * random splits have slope 4 on every mode.

## Scope

* This is an operator (Choi-state) statement, so it does not depend on the initial state.
* Rates of state-entanglement growth for generic initial states (Bianchi–Hackl–Yokomizo, JHEP 03 (2018) 025) are a
  different quantity and are not claimed.
* The finite-field analogue is the export count of Pass 11207. A stabilizer Choi state's entanglement across a qutrit
  cut is rank-valued, which is how the arrow was defined.
