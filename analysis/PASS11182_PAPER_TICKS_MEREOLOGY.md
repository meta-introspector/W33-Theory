# Pass 11182 — the paper's clock is a product of three independent qutrit clocks in 108 splits, none of them the light-cone coordinates; the perfect tick is local in none

Producer: `analysis/w33_pass11182_paper_ticks_mereology.py`
Regression: `tests/test_w33_pass11182_paper_ticks_mereology.py`

Profiles over the 110 565 three-qutrit tensor factorisations (Pass 11180 method):

| tick | order | local | perfect | partial |
|---|---|---|---|---|
| clock K (Theorem 4.3), K², its Fourier-dual kick V | 3 | **108** | 8748 | 101 709 |
| perfect interacting tick V(N) K V(N), N = (x₀+x₁+x₂)² | 9 | **0** | 4860 | 105 705 |
| F₉ gate K = P² + i𝟙 | 3 | 27 | 3564 | 106 974 |

Random ticks: 2/3 of a sample of 300 are local in no split (exact value 7922/12285, Pass 11181).

**The clock's own subsystems.** In all 108 splits where the clock is local, it fixes each plane: it is literally three
independent single-qutrit clocks there.
* **None** of the 108 contains a light-cone coordinate plane, and only 3 contain the transverse one.
* Exactly 4 keep positions and momenta separate, and these are the **4 orthogonal frames of the light-cone metric q**,
  the only q-orthogonal frames of F₃³.
* The other 104 mix positions with momenta.

**Readings.**
* The paper writes its clock in light-cone coordinates, where it entangles (Pass 11166: one trit between the light-cone
  qutrits). Split along the principal axes of the Lorentzian form, or along 104 mixed splits, the same clock entangles
  nothing.
* The perfect interacting tick is local in no split at all. The interaction that makes it scramble maximally also
  removes every subsystem description in which it would be local. Its order is 9, and its class is fixed-point-free
  (Pass 11181).
