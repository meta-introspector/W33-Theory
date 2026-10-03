# Pass 11354 — the Clifford+T walk on PU(3) representations: an exact 3/8 on low degrees, slower decay for reversibility

Producer: `analysis/w33_pass11354_decay_rate.py`
Certificate: `data/w33_pass11354_decay_rate.json`
Regression: `tests/test_w33_pass11350_11354.py`

**Setup.**
* The one-qutrit random word U_k = C_k T ⋯ C₁ T has law μ^{*k}.
* On a representation π of PU(3), its Fourier transform is (P_π π(T))^k, where P_π projects onto Clifford-invariant
  vectors.
* The eigenvalue 1 occurs exactly on the SU(3)-invariants. The next modulus ρ is the convergence rate in that
  representation.
* Computed on the tensor spaces V^{⊗p} ⊗ V̄^{⊗q} with p ≡ q mod 3 (so the centre acts trivially), up to dimension 2187.

| (p, q) | dim | Clifford invariants | SU(3) invariants | ρ |
|---|---|---|---|---|
| (1,1) | 9 | 1 | 1 | 0 (only the trivial summand; the adjoint 8 has no Clifford invariant: 2-design) |
| (3,0) | 27 | 1 | 1 | 0 |
| (2,2) | 81 | 2 | 2 | 0 |
| (4,1) | 243 | 3 | 3 | 0 |
| (6,0) | 729 | 5 | 5 | 0 |
| (3,3) | 729 | **7** | 6 | **3/8** |
| (5,2) | 2187 | **15** | 11 | **3/8** |

**Results.**
* The first Clifford-anisotropic invariants appear in degree (3,3), where the Clifford group stops being a design.
  There the walk contracts by **exactly 3/8 per cubic gate** (0.375000000000000), and (5,2) gives the same value.
* Pass 11312's reversible fraction decays at **about 0.72 per gate** (successive ratios 0.65 → 0.76 → 0.67 → 0.73 →
  0.70 → 0.71 → 0.72), much **slower** than 3/8.
* **Reading.** Low-degree observables equidistribute fast. The reversible fraction is the measure of a neighbourhood of
  a closed Haar-null set, a singular observable whose decay is controlled by high-degree representations. So its rate
  is not given by the first gap.
* **Scope.** Finitely many representations. Neither a uniform spectral gap nor the 0.72 is derived here.
