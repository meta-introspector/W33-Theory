# Pass 11134 — the zero-radius tachyon is one complex species that exists only through winding

Producer: `analysis/w33_pass11134_dual_tachyon_species.py`
Scan: `analysis/w33_pass11134_scan_dual_tachyon.py`
Frozen: `data/w33_pass11134_dual_tachyon_multiplet_21.json`
Regression: `tests/test_w33_pass11134_dual_tachyon_species.py`

As the neutral torus shrinks (Pass 11128), every winding class becomes light. The Δ → −½ ground states are the vectors
l = π + V0 + N·W with l² = 1 and β-phase −1. Here W is the neutral torus's Wilson line, and N is taken mod 3 because 3W
lies in the lattice.

In all 21 neutral-exit survivors:
* **One complex species.** There are exactly two such vectors: N = 1 and its conjugate N = 2. None has N = 0.
* **It exists only through winding.** It is a momentum state of the T-dual circle; its KK tower is the 4 → 6 → 8 growth
  seen in Pass 11128.
* **Split across the two E8s:** l₁² = 2/9 (13 models), 5/9 (5), 4/9 (2), 7/9 (1). These are ninths, reflecting the Z3
  Wilson line.
* **Not a vector-class multiplet.** It is neither the 32 of the SO(32) tachyonic string nor the 16 of the SO(16)×E8 one.
  The Wilson lines keep one component of whatever 10D multiplet it descends from.

**Correction to Pass 11128.** l² = 1 holds, but it does not identify a 10D string. The identification stays open.
