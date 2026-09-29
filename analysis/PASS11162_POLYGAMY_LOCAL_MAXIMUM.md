# Pass 11162 — 4/√15 is a certified local maximum of spatial polygamy; 96 restarts find nothing higher

Producer: `analysis/w33_pass11162_polygamy_local_maximum.py`
Scan: `analysis/w33_pass11162_scan_restarts.py`, frozen as `data/w33_pass11162_restarts.json`
Regression: `tests/test_w33_pass11162_polygamy_local_maximum.py`

The test point is Pass 11159's state ψ* = √(7/15)|000⟩ + √(2/15) Σ_{a=1,2} |a⟩(|0a⟩ + |a0⟩), with
N_AB + N_AC = 4/√15.

1. **Smoothness.** The partial transposes of ρ_AB and ρ_AC at ψ* have spectrum
   {−0.1915 (×2), −2/15, 2/15 (×3), 0.3249 (×2), 7/15}. No eigenvalue lies within 2/15 of zero, so the negativity is
   real-analytic near ψ* and a second-order test is legitimate.
2. **Hessian.** The test is on the 52-dimensional real tangent space of the unit sphere modulo phase.
   * The gradient is ~3·10⁻⁹.
   * The Hessian has **32 strictly negative eigenvalues** (−4.58 … −0.72), **0 positive**, and **exactly 20 zero**
     eigenvalues.
   * 20 is the dimension of the local-unitary orbit through ψ*, computed independently from the orbit's tangent vectors.
   * So ψ* is a Morse–Bott maximum: strict modulo local unitaries. The margin, 0.72, is far above the finite-difference
     error (~10⁻⁷).
3. **Global evidence.** 96 random restarts all converge to 4/√15 to 10⁻¹⁵, and none exceeds it. Each restart runs
   Nelder–Mead then Powell from a random scale in 0.2–3. These add to Pass 11148's 16 restarts and the Pass 11154 scan.

**Scope.** This certifies a local maximum and gives strong numerical evidence for a global one. It is **not** a proof
of global optimality, which stays open.
