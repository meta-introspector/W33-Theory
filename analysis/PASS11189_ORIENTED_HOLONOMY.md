# Pass 11189 — the arrow entanglement cannot see: oriented holonomy on split relations

Producer: `analysis/w33_pass11189_oriented_holonomy.py`
Certificate: `data/w33_pass11189_oriented_holonomy.json`
Regression: `tests/test_w33_pass11189_oriented_holonomy.py`

**What was open.** Pass 11184 classified how two three-qutrit splits relate (the 20 orbitals of Sp(6,3) on pairs of
factorisations) up to two ties. Each tie is a relation and its time reverse (orbital sizes 256 and 6912).

**Why nothing unoriented can separate them.** The reverse relation has T⁻¹ with blocks (T⁻¹)ⱼᵢ = adj(Tᵢⱼ). Every loop
trace is therefore unchanged by reversal. Ranks, determinants, incidences and traces are all time-symmetric.

**Oriented data.** SL(2,3) fuses neither a unipotent u with u⁻¹ nor a nilpotent N with −N. These pairs are
GL(2,3)-conjugate but not SL(2,3)-conjugate. Reversing a loop replaces its holonomy
P = T_{r₁c₁} adj(T_{r₂c₁}) T_{r₂c₂} ⋯ with adj(P), which is P⁻¹ or −N. So the SL(2,3)-class of each rooted,
oriented walk holonomy is a local invariant that carries a direction.

**Result.**
* **The 6912 pair** is separated by holonomy on the shortest walks, the 4-cycles: its nilpotent holonomies are N in
  one relation and −N in its reverse.
* **The 256 pair has no holonomy at all.** Every walk of length ≤ 12 has holonomy 0, because all six rank-one blocks
  are the same nilpotent-type block. It is separated by a symplectic *twist* instead. Let u_i span row i's common
  image, R_i be the invertible block of row i, and E = u_{i'} ⊗ φ be a rank-one block of column σ(i). Then
  ψ_i = ω(u_i, R_i ·) = ρ(i,i′)φ. Define Q = ρ(0,1)ρ(1,2)ρ(2,0) ∈ {±1}.
  * Q is invariant under rescalings and local SL(2,3), and equals the product for the opposite cyclic order.
  * It changes sign under time reversal (one ω per factor, three factors).
  * Q = −1 on all 256 members of one orbital and +1 on all 256 of its reverse.
* **Holonomy multiset + twist classify all 20 relations.** This is a complete invariant of how two three-qutrit splits
  relate, *including the direction of time*.
* **Time reversal is anti-unitary.** τ = (x,z) ↦ (x,−z) is complex conjugation, an anti-unitary that fixes the
  standard split. It exchanges exactly the two pairs that no unoriented invariant separates (256 ↔ 256,
  6912 ↔ 6912). It fixes both members of the 2304 pair, which the incidences already separate, and every
  self-paired relation.

**Reading.** Entanglement data (ranks, mutual informations, incidences) classify split relations only up to time
reversal. The missing arrow is an orientation: the SL-class of a nilpotent holonomy, or the sign of a symplectic
twist. Anti-unitary time reversal flips exactly that.

**Prior art.**
* Paired orbitals (HgH vs Hg⁻¹H) are classical permutation-group theory.
* For two-qubit gates, local invariants distinguish U from U†: the two halves of the Weyl chamber (Zhang, Ye & Guo,
  PRA 71, 062331 (2005)). Operator entanglement is invariant under U → U† (Balakrishnan & Sankaranarayanan,
  arXiv:1005.2467).
* Anti-unitary Cliffords are the anti-symplectic maps (Appleby, arXiv:quant-ph/0412001).
* The repo already had transpose pairs of frame orbitals (Pass 1082) and an outer element swapping paired orbitals
  (Pass 4795).
* The holonomy/twist classification of three-qutrit split relations was not found.
