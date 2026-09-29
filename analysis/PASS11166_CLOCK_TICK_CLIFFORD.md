# Pass 11166 — the spectral clock tick is the Clifford propagator of the light-cone quadric; no translation-invariant clock is perfect

Producer: `analysis/w33_pass11166_clock_tick_clifford.py`
Regression: `tests/test_w33_pass11166_clock_tick_clifford.py`

**The clock.** This is the clock of Theorem 4.3 (`w33_20260924_null_history_spectral_clock.py`).
* The 27 events are F₃³, with edges along the 8 null vectors of q(v) = v0v2 − v1².
* The tick is U = exp(−2πiL/9), with U³ = I.
* We read the 27 events as three qutrits, one per coordinate.

**U is a Clifford gate.** It sends each of the 729 Paulis to a single Pauli. Up to phases:

    X_a → X_a  (all translations conserved)
    Z0 → Z0 X2,   Z1 → Z1 X1,   Z2 → Z2 X0

So S = [[I, 0], [M, I]] in (x | z) form, where M is the Gram matrix of q's polar form,
B(u,v) = u0v2 + u2v0 + u1v1. The tick is the **free propagator e^{−iq(p)} of the Lorentzian form**:
* it couples the two light-cone coordinates;
* it acts locally on the transverse coordinate.

**Party blocks.**

| | block ranks | block dets |
|---|---|---|
| diagonal | 2, 2, 2 | 1, 1, 1 |
| (0,2), (2,0) | 1 | 0 |
| others | 0 | 0 |

The column law Σ_i det S_ij = 1 holds. The tick is **not perfect**: it exchanges one trit between the light-cone
qutrits and none with the transverse one.

**Theorem (elementary).** A Clifford gate that commutes with all translations X_a has S = [[I, 0], [M, I]] with M
symmetric.
* Every off-diagonal party block is then [[0, 0], [m, 0]]: rank ≤ 1, det 0.
* So **no translation-invariant (spectral, Cayley-graph) clock is ever perfect**. Conserved momentum caps scrambling at
  one trit per pair per tick.
* This was checked for all 3⁶ = 729 symmetric M: maximum off-diagonal rank 1, and 0 perfect.
* It is the Clifford form of the known fact that conservation laws slow scrambling.

**Information flow.** These are brute-force Choi entropies, equal to the stabilizer formula of Pass 11165.
* The transverse qutrit keeps its whole past: I(in1 : out1) = 2.
* Each light-cone qutrit keeps one trit in its own future (I(in0 : out0) = 1), and the full two trits are held by the
  light-cone pair.
* So forgetting a partner after a clock tick erases **at most one trit**. The clock does not secret-share, in contrast
  with Pass 11165.

**Cross-reference.** The quadric q = v0v2 − v1² is the null cone of Pass 10946: its four projective null directions,
the images of P¹(F₃) under (x,y) ↦ (x², xy, y²), are the four Hesse/tetracode clocks. The spectral clock's propagator
uses the same quadric's polar form, and M = M⁻¹, so q and its dual have the same Gram matrix.

**Reading.** The paper's clock is reversible bookkeeping on the light cone, not a scrambler. The two-trit arrow needs a
tick that breaks translation invariance, i.e. an interaction.
