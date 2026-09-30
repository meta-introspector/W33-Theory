# Pass 11178 — three-qutrit factorisations: a rank-20 geometry; the Choi signature sees 15 of its 20 relations; perfect = one orbital

Producer: `analysis/w33_pass11178_three_qutrit_orbitals.py`
GAP: `analysis/gap/w33_pass11178_orbital_signatures.g`, frozen as `data/w33_pass11178_gap_orbital_signatures.txt`
Regression: `tests/test_w33_pass11178_three_qutrit_orbitals.py`

**GAP result.** Sp(6,3) acts on the 110 565 three-qutrit factorisations F₃⁶ = P₁ + P₂ + P₃ with **rank 20**. The
subdegrees are 1, 36, 96, 256, 256, 768, 864, 2304 (×3), 3456, 4608, 6912 (×6), 10 368, 41 472.

**Signature of each orbital.** For a representative F' = S(F₀), block (i,j) is the gate's party block S_ij, coded as
zero / rank one / invertible det +1 / invertible det −1. Up to relabelling rows and columns, this is the gate's
Choi-entanglement signature.

**Results.**
* **The 20 orbitals give 15 signature classes.**
  * 13 orbitals are determined by their signature.
  * One class merges two orbitals (6912 + 6912).
  * One class merges five (256 + 256 + 2304 + 6912 + 6912 = 16 640): "one orientation-preserving invertible block in
    each row and column, every other block rank one".
* **The perfect relation is the single orbital of size 3456.** Every block is invertible, with one orientation
  reversal per row and column (Pass 11169). And 3456/110 565 = 128/4095, as in Pass 11163.
* A fresh sample of 300 000 elements of Sp(6,3) matches all 15 class frequencies to within 2.6σ, with no missing or
  extra classes.

**Comparison.** For two qutrits (Pass 11177) the action has rank 3 and the signature is a complete invariant. For three
qutrits the action has rank 20 and the signature misses 5 of the 20 relations. The finer invariant that separates them
is open.
