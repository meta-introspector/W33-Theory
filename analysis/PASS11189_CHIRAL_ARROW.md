# Pass 11189 — the direction of time between two subsystem splits is a chirality, and only a Berry phase sees it

Producer: `analysis/w33_pass11189_chiral_arrow.py`
Frozen: `data/w33_pass11189_chiral_arrow.json`
Regression: `tests/test_w33_pass11189_chiral_arrow.py`

**The question (Pass 11184).** The 20 relations between two three-qutrit splits are classified by their entanglement data,
except for two reversal pairs, {256, 256} and {6912, 6912}. For each of these, a gate and its inverse look identical to
every measure used there.

**1. Time reversal is a mirror.** Complex conjugation of qutrit states acts on the symplectic data as τ: (x, z) ↦ (x, −z)
on every qutrit. This is the anti-symplectic outer automorphism of PSp(6, 3), and it fixes the standard split. Computed
on the 20 orbitals:
* τ **fixes 16** and **swaps exactly the two tied reversal pairs**. For those pairs Oᵀ = τ(O): the reversed relation is the
  mirror image of the relation. The two members of each pair are enantiomers.
* The third reversal pair (2304) is τ-fixed. That is why the signature already separates it.
* Consequence: any invariant blind to the sign of ω cannot separate the tied pairs. This rules out ranks, image and
  kernel points, incidences, block determinants and loop-holonomy classes. Reversing a loop acts on unipotent
  holonomies the same way τ does.

**2. The chiral invariant.** The Kashiwara–Maslov index τ(L₁, L₂, L₃) of three Lagrangians is the Witt class in
W(F₃) = ℤ/4 of Q = ω(x₁,x₂) + ω(x₂,x₃) + ω(x₃,x₁) on L₁ ⊕ L₂ ⊕ L₃, computed as rank − 2[disc = 2] mod 4. It is cyclic,
alternating, and odd under τ. For a pair (F₀, F):
* take the 64 product Lagrangians of each split;
* give each F₀-Lagrangian its *profile*: its sorted intersection dimensions with the 64 F-Lagrangians;
* count the ordered triples (L₁, L₂, M) with profile(L₁) > profile(L₂). This orientation is gauge-invariant, so
  alternation cannot cancel it.

The resulting spectrum is constant on orbitals (checked on 4 members of each). Its net chirality
χ = #{class 1} − #{class 3} is:

| orbitals | χ |
|---|---|
| all 16 τ-fixed | 0 (spectrum exactly symmetric under negation, as τ-invariance forces) |
| 256 pair | −54 / +54 |
| 6912 pair | +18 / −18 |

With Pass 11184's signature and fine invariant, **(signature, fine invariant, χ) classifies all 20 relations.**

**3. It is a Berry phase.** For the zero-shift stabilizer states |L⟩, with projector (1/27) Σ_{v∈L} D(v) and D the
symmetric Weyl operators,

    tr(P_{L₁} P_{L₂} P_{L₃}) = |tr(…)| · (−i)^{τ(L₁,L₂,L₃)}      (200/200 random mixed triples agree).

This is the finite-field case of the classical identity between the Maslov/Kashiwara index and the phase of the Weil
representation (Lion–Vergne; T. Thomas, *Weil representation and transfer factor*; de Gosson, *Expo. Math.* 2015). So
χ is minus the net imaginary part of the Pancharatnam phases of the stabilizer-state triangles that the two splits span.

**Reading.** For three qutrits, the relation between two subsystem descriptions has a direction in time exactly when it
is chiral. That direction is recorded in geometric phase and nowhere else: the entanglement data of the relation are
time-symmetric, and the orientation lives in the sign of i.

**Scope.**
* Constancy of the spectrum is checked on 4 members per orbital. It is not proved for every member. It does follow from
  H-invariance, since the spectrum is built from H-equivariant data.
* The Bargmann identity is checked numerically on random triples here; it is not proved here.
* What this pass adds is the τ-swap of the tied pairs and the chirality values; the identity itself is classical.
