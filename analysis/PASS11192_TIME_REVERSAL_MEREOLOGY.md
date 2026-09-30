# Pass 11192 — time reversal on the substrate: the 36 reflections of W(E6) are the maximal arrows

Producer: `analysis/w33_pass11192_time_reversal_mereology.py`
Certificate: `data/w33_pass11192_time_reversal_mereology.json`
Regression: `tests/test_w33_pass11192_time_reversal_mereology.py`

**Setting.**
* Aut W(3,3) = PSp(4,3).2 ≅ W(E6). Its outer coset is the anti-symplectic maps (TᵗJT = −J), realised by
  anti-unitaries: Wigner's time reversals.
* The identification of anti-unitary Cliffords with anti-symplectic maps is Appleby (arXiv:quant-ph/0412001), and
  W(E6) as the symmetry of the two-qutrit Pauli graph is Planat 2011 (J. Phys. A 44, 045301). The repo already owns
  the combined statement (Forty Points Cor. 9.7).
* New here: the mereology and arrow of these time reversals on the 45 two-qutrit splits.

**Definitions.**
* local(T) is the number of splits in which T is block-monomial.
* A(T) is the minimum over the 45 splits of the rank-export. For an anti-unitary this equals the export of the unitary
  τ_F T, where τ_F is the local complex conjugation of the split F.

**Result (all 25 920 projective anti-symplectic elements).**
* A ∈ {0, 2, 4}, with 9180 / 16 704 / **36** elements. Unitary ticks only ever reach A = 2 (Pass 11183).
* **A = 4 occurs for exactly one conjugacy class, of size 36: the 36 reflections of W(E6).** In every one of the 45
  splits a reflection exports everything (both qutrits, both trits), although it is block-monomial, i.e. swaps the two
  factors, in 15 of them.
* The other outer involution class (size 540, local in 7 splits) has A = 0.
* **10 944 anti-unitaries are local in no split.** That is the same fraction, 19/45, as for unitary ticks (Pass 11180).

**Reading.** The time reversals that no choice of subsystems can tame are exactly the reflections of E6. This ties the
maximal arrow to the E6 root system: 36 reflections ↔ 36 positive roots. It is a counting statement on the substrate's
own symmetry group, and no dynamics is claimed.
