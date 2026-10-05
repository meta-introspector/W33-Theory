# Pass 11492 — the lowest time-odd invariant of a qudit state, by dimension: 9, 6, 6, 5, 5 for d = 2, 3, 5, 7, 11; and the two-qutrit register's unique one, which restricts to h₆/√108 (1/108 = 1/120 + 1/1080, derived)

Producer: `analysis/w33_pass11492_time_odd_by_dimension.py`
Certificate: `data/w33_pass11492_time_odd_by_dimension.json`
Regression: `tests/test_w33_pass11486_11492.py`

**Question.**
* Pass 11357 measured, for unitaries, the lowest Clifford witness degree of time reversal: 10, 6, 4 for d = 2, 3, 5.
* Pass 11491 gave the qutrit **state** answer: degree 6, with the full series.
* This pass asks the state question in other dimensions. What is the lowest degree k of a Clifford-invariant polynomial
  of bidegree (k, k) in (ψ, ψ̄) that changes sign under time reversal?

**Method (exact, as in Pass 11491).** This is the method of Pass 11370 (twisted Frobenius–Schur counts, Kawanaka–Matsuyama
1990), which counted τ-odd invariants of **unitaries** in the matrix coefficients of PU(d). Here it is applied to **states**:
the representation is Sym^k ⊗ Sym^k-bar, and the twist is conj(g)·g.
* Molien-type sums D_k (all invariants) and Tw_k (trace of time reversal) over the projective one-qudit Clifford group.
  Its order is 24, 216, 3000, 16 464 and 159 720 for d = 2, 3, 5, 7, 11.
* Eigenvalues are identified as N-th roots of unity, with N = 12, 12, 60, 168, 660.
* Every cyclotomic sum is evaluated at 60 digits and asserted to be an integer divisible by |G|.
* **Two independent group constructions.**
  * For d = 2 and 3: Pass 11357's BFS.
  * For d ≥ 5: W(a)·V_S, with V_S from the Weyl twirl Σ_q W(Sq) A W(q)†.
  * Both constructions give identical series at d = 5, the control.

## Results

| d | \|G\| | design strength | lowest odd degree | odd invariants there | O_k from that degree on |
|---|---|---|---|---|---|
| 2 | 24 | 3 | **9** | 1 | 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 7, … |
| 3 | 216 | 2 | **6** | 1 | 1, 2, 3, 5, 8, 11, 16, 22, … |
| 5 | 3000 | 2 | **6** | 4 | 4, 12, 31, 71, 146, 280, … |
| 7 | 16 464 | 2 | **5** | 2 | 2, 18, 72, 243, 711, 1865, … |
| 11 | 159 720 | 2 | **5** | 16 | 16, 171, 1102, 5803, 26 319 |

**The two-qutrit register (the substrate's own 3 ⊗ 3).**
* The projective two-qutrit Clifford group has 81 × 51 840 = 4 199 040 elements. It is streamed one Weyl translation at a
  time.
* There are 167 eigenvalue-ratio classes and 11 twisted classes, with N = 360.

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| all invariants D_k | 1 | 1 | 1 | 2 | 3 | 6 | 14 | 32 | 87 |
| time-odd O_k | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 6 | 25 |

* **The register's first time-odd invariant has degree 6 and is unique**, exactly as for one qutrit. A 9-dimensional
  system thus behaves like d = 3, not like a generic d ≥ 5.
* **What it is (numerically exact, mechanism open).**
  * Because the odd space is one-dimensional, the frame-potential difference
    Δ⁽²⁾(Ψ) = avg_C |⟨Ψ|CΨ⟩|¹² − avg_C |⟨Ψ|CΨ̄⟩|¹² (over all 4 199 040 Cliffords) is a multiple of |h⁽²⁾(Ψ)|².
  * On Ψ = ψ ⊗ s, with s a stabilizer state:

  > **Δ⁽²⁾(ψ ⊗ s) = Δ₆(ψ) / 108**, to 10⁻¹⁰ relative.

  * This holds for 3 different stabilizer states, for random ψ, and with the factors in either order.
  * So the register's arrow, seen from one qutrit beside a stabilizer state, **is** the single-qutrit chirality h₆ of
    Pass 11434 (Δ₆ = 2|h₆|², Pass 11419).
  * **The factor 108 is derived.**
    * **Compress.** Write A_C = (1 ⊗ ⟨s|) C (1 ⊗ |s⟩), so that ⟨ψ⊗s|C(ψ⊗s)⟩ = ⟨ψ|A_C ψ⟩ (and likewise with ψ̄, since s
      is real).
    * **Classify.** Over all 4 199 040 Cliffords (`derivation_108` in the certificate), A_C is:
      * a one-qutrit Clifford on exactly **1/120** of them;
      * 3^(−1/2) times a one-qutrit Clifford on exactly **27/40**;
      * rank one or zero on the rest.
    * **Uniformity.** In both Clifford cases the induced Cliffords cover the 216 one-qutrit Cliffords **exactly
      uniformly** (equal counts).
    * **Cancellation.** The rank-one terms cancel in the odd difference (relative sum 3·10⁻¹⁶). The symmetry pairing
      them is not written out, so this step is numerical.
    * **Conclusion.** Δ⁽²⁾(ψ⊗s) = (1/120 + (27/40)·3⁻⁶)·Δ₆(ψ) = (1/120 + 1/1080)·Δ₆(ψ) = **Δ₆(ψ)/108**.
    * The direct evaluation agrees to 1.2·10⁻¹¹.
  * On generic product states ψ ⊗ φ the relation is **open**. The values are recorded in `generic_products`; the one
    combination tried, |√Δ₆(ψ) ± √Δ₆(φ)|, does not fit. Even degree-6 invariants of each factor presumably enter as
    weights, but that is not tested.

**Closed forms, proved by the quasi-polynomial bound (period 12; (2d − 1)·12 coefficients checked).**
* **Qubit:** O(t) = t⁹ / ((1−t)(1−t⁴)(1−t⁶)), D(t) = (1 − t³ + t⁶) / ((1−t)³(1+t)(1+t²)(1+t+t²)),
  Tw(t) = (1 + t³ + t⁶) / ((1−t)²(1+t)²(1+t²)(1−t+t²)).
  * This is a control, not news. The qubit Clifford group acts on the Bloch sphere as the octahedral rotation group,
    and time reversal is a reflection. The odd invariants are the octahedral pseudoscalar xyz(x²−y²)(y²−z²)(z²−x²)
    (degree 9) times the even ones.
* **Qutrit:** reproduces Pass 11491 term for term.
* **d = 5, 7, 11:** coefficients exact to the printed degree. No closed form is claimed, because the period bound needs
  more terms than were computed (9·60 for d = 5).

## Reading

1. **Design strength re-derived as a control.** D₂ = 1 for every d, and D₃ = 1 only for d = 2. This re-derives the
   classical facts that the qubit Clifford group is a 3-design and the odd-prime ones are exactly 2-designs.
2. **The lowest odd degree falls with d: 9, 6, 6, 5, 5.**
   * It is not constant over odd primes (6 at p = 3, 5 but 5 at p = 7, 11).
   * The first odd invariant is unique only for d = 2 and 3.
   * The qutrit is the last dimension with a single, canonical time-odd state invariant: h₆ of Pass 11434, a chirality
     of the four MUBs.
   * From d = 5 on there are several independent ones.
3. **No odd invariant of degree ≤ 4 exists in any dimension computed.** A mechanism is **not** claimed.
   * A tempting one is wrong. The triangle invariants I_s = Σ_{ω(p,q)=s} χ_p χ_q χ_{−p−q} satisfy I_s = I_{−s}
     identically. This holds for mixed states as well (checked at p = 3, 5, 7, 11), because swapping two vertices
     reverses the area.
   * So it is not purity, and the cubic triangle invariants are simply even.
   * Whether 5 is the floor for all large p is **open**. p = 13 needs a 369 096-element group; it is feasible with
     streaming.
4. **Comparison with unitaries (Pass 11357: 10, 6, 4).**
   * The two ladders agree only at d = 3 (both 6).
   * They measure different objects: a state invariant has bidegree (k, k) on ℂ^d, while J_{2t} is a moment of U.
   * No relation between them is claimed.

**Prior art.**
* The counting method is Pass 11370's (for unitaries) and Kawanaka–Matsuyama's.
* The qubit case is classical octahedral invariant theory.
* Clifford-group design strengths are classical (Zhu; Webb; Kueng–Gross).
* The time-odd / time-even split of Clifford ray invariants for d ≥ 3 was not found in the corpus. The d = 3 sequences
  are not in the OEIS (searched 2026-10-04).
