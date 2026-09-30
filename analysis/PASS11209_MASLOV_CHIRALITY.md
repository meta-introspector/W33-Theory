# Pass 11209 — the direction of time between two subsystem splits is a chirality, and a Maslov/Berry phase sees it

Producer: `analysis/w33_pass11209_maslov_chirality.py`
Frozen: `data/w33_pass11209_maslov_chirality.json`
Regression: `tests/test_w33_pass11209_maslov_chirality.py`

**Numbering and relation to master.** This was first committed on the branch as Pass 11189 (2026-09-30, 19:10 UTC). It
is renumbered because master's own Pass 11189 (`analysis/PASS11189_ORIENTED_HOLONOMY.md`, from the parallel session)
settles the same question by a different invariant.

**What master's Pass 11189 has.**
* The tied relations of Pass 11184, the pairs {256, 256} and {6912, 6912}, are separated by oriented data: the
  SL(2,3)-class of oriented walk holonomies (N versus −N) for the 6912 pair, and a symplectic twist Q = ±1 for the 256
  pair.
* Complex conjugation τ = (x, z) ↦ (x, −z) exchanges exactly these two pairs and fixes every other relation.

This branch reached the same τ statement independently: τ fixes 16 of the 20 orbitals and swaps exactly the two tied
pairs, so Oᵀ = τ(O) and the members of each pair are enantiomers.

**Point of precision.** An invariant blind to the sign of ω cannot separate the tied pairs. Examples are ranks, image and
kernel points, incidences, block determinants, and the *direction-blind* multiset of loop holonomies taken over both
orientations of every loop: reversing a loop sends a holonomy to its adjugate, and on SL(2,3) conjugation by τ does the
same up to conjugacy. Master's holonomy invariant works because its walks carry an orientation, not because holonomies
separate as an unoriented multiset.

**What this pass adds.**
* A second, uniform chiral invariant, the Kashiwara–Maslov spectrum. It needs no case split between the two pairs.
* Its identification with the Bargmann (Pancharatnam–Berry) phase of stabilizer states.

Master's reservation of 11199–11206 lists "the twist as a discrete Bargmann phase". The Bargmann/Maslov identity below
is the finite-field ingredient that item needs. Whether master's twist Q equals a specific Maslov index is not decided
here.

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
