# Pass 11176 — temporal entanglement through a tick is all-or-nothing; perfect ticks cause amnesia with an echo; Bell violation survives without temporal entanglement

Producer: `analysis/w33_pass11176_temporal_echo.py`
Regression: `tests/test_w33_pass11176_temporal_echo.py`

**Setup.** One qutrit of the 27-event register is the system; the other two are its environment, maximally mixed. The
temporal negativity between times separated by n ticks is that of the pseudo-density operator (Pass 11148) of the
n-tick channel.

**Theorem (all-or-nothing).** For a Clifford tick S, a qutrit keeps temporal negativity 1 exactly when its column of
blocks is local (S_iq = 0 for i ≠ q); otherwise it has 0.
* Proof: the surviving Paulis form a subgroup of F₃² of order 1, 3 or 9. With cube-root phases, the pseudo-density
  eigenvalues are ∝ 1 + 2cos(2πm/3) ≥ 0 unless all 9 survive.
* Checked on 120 (qutrit, tick) pairs.

**Results.**

| tick | temporal negativity per qutrit | three-time sum 2N(1) + N(2) |
|---|---|---|
| identity | 1, 1, 1 | 3 |
| clock K (order 3) | light-cone 0, transverse 1; full revival at tick 3 | 0, 3, 0 |
| perfect tick V K V, N = (x₀+x₁+x₂)² (order 9) | 0 for eight ticks, 1 at tick 9 | 0 |
| random perfect ticks | 0 after one and after two ticks | 0 |

**Echo census** over all 108 perfect kick–tick–kick ticks:
* the orders run from 7 to 36;
* the transverse qutrit's entanglement with its past often returns first, e.g. at tick 3 of 12, 4 of 24, 9 of 36.

**Bell without temporal entanglement.** The two-time CGLMP value (d = 3) through one clock tick:
* transverse qutrit: 3.1628, the identity value;
* light-cone qutrits: **3.0606 > 2, although their temporal negativity is 0**;
* through a perfect tick: exactly 0, because the channel is completely depolarising.

**Reading.**
* In time, a Bell violation is not a witness of entanglement. It witnesses the setting-dependent memory left by the
  first measurement (compare Pass 11173).
* Scrambling destroys a qutrit's entanglement with its own past outright, not gradually.
* What survives, and returns as an echo, is decided by when a power of the tick becomes local again for that qutrit.
