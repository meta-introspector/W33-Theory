# Pass 11250 — the qubit protection limit is unconditional: both identities are theorems in the literature

Producer: `analysis/w33_pass11250_qubit_identities_proved.py`
Certificate: `data/w33_pass11250_qubit_identities_proved.json`
Regression: `tests/test_w33_pass11250_11254.py`

**Starting point.**
* Pass 11244 (parallel session) gave the qubit limit lim E_n[c] = **0.66516065**. It was conditional on two identities
  about the unipotent elements of Sp(2m,2), verified on every Jordan type with m ≤ 6 but not proved.
* Pass 11245 proved the counting lemma behind identity 2 and reduced that identity to a uniformity statement.

## Identity 1 is a theorem (cited, not re-derived)

* **Pieces are Jordan strata.** For Sp_{2n} in characteristic 2, two unipotent classes lie in the same Lusztig
  "unipotent piece" iff they have the same Jordan type (Xue, *On unipotent and nilpotent pieces for classical groups*,
  arXiv:0912.3820, Lemma 2.7, citing Lusztig).
* **Point counts are characteristic-free.** The number of F_q-points of a piece is a polynomial in q with integer
  coefficients, independent of the characteristic. This is Lusztig, *Unipotent elements in small characteristic*,
  Transform. Groups 10 (2005), property P5, proved for Sp in §3; Xue §6 restates it for all F_{p^s}.
* **Odd q.** In odd characteristic the piece of type λ is the single geometric class λ. Its F_q-points are the rational
  classes, which is the odd-q formula summed over the signs of the even parts.
* Evaluating the same polynomial at q = 2 gives identity 1.

## Identity 2 follows from Lusztig's own parametrisation

* In the proof of P5–P8 for Sp with p = 2 (op. cit. §3), the fibre E is identified with triples (w, α, j).
* Each j_n ranges **independently over all nondegenerate symmetric bilinear forms** on the multiplicity space of the
  size-n blocks.
* The classes inside a piece are distinguished exactly by which j_n are symplectic (alternating).
* Pass 11244's χ₂ = 0 says that ω(v, (u−1)v) vanishes on ker(u−1)². That is the statement that j₂ is alternating:
  * on ker N², only the size-2 blocks contribute, because a block of size k ≥ 3 involves N^{2k−3}w = 0;
  * a symmetric form in characteristic 2 is alternating iff its diagonal vanishes.
* So the χ₂ = 0 fraction is #{alternating}/#{nondegenerate symmetric} at rank m₂. Pass 11245's lemma gives this as
  2^{−m₂} for m₂ even and 0 for m₂ odd.
* More generally, for even q the fraction is q^{−m₂} for m₂ even, by the same count with Lusztig's recursion for
  symmetric forms.

**Consequence.** The qubit limit **0.66516065** (rigorous tail bracket [0.6651604993, 0.6651614201], Pass 11244) is no
longer conditional on unverified identities. It rests on the cited theorems.

## Independent check in a new field (q = 4)

Lusztig states P5 for split structures with q − 1 "sufficiently divisible". q = 2 was already checked by Pass 11244 for
m ≤ 6. As an independent test of the characteristic-free statement, all of Sp(4,4) is enumerated here by brute force.

| | value |
|---|---|
| elements enumerated | **979 200** = \|Sp(4,4)\| |
| unipotents | **65 536** = q⁸ (Steinberg) |
| identity 1, Jordan types (1⁴), (2,1²), (2²), (4) | **1, 255, 4080, 61200**, each equal to the q = 4 formula |
| identity 2, type (2,1²) | χ₂ = 0 fraction **0** (predicted 0) |
| identity 2, type (2²) | χ₂ = 0 fraction **1/16 = q⁻²** (predicted 3/48) |

**Rediscovery note.** This pass found that the two "open" identities are published. The lesson is the repo's usual
one: search the literature for the *result* (the piece / Jordan-stratum point count) before treating it as new.
