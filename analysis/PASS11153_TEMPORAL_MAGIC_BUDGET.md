# Pass 11153 — in time the magic must sit in the measurements; in space it can sit in the state

Producer: `analysis/w33_pass11153_temporal_magic_budget.py`
Regression: `tests/test_w33_pass11153_temporal_magic_budget.py`

| where the magic sits | space (two qutrits) | time (one qutrit, two Lüders measurements) |
|---|---|---|
| **in the state** (Pauli measurements, any state) | **2.6017**, violates (11149, 11152) | **2 exactly**, no violation (exhaustive) |
| **in the measurements** (stabilizer state/start, any bases) | 2.8729 (\|Ω⟩, 11144) | 3.1628 (any pure start) |

* **Time, Pauli measurements.** For fixed bases the temporal value is linear in the start state, so its maximum is the top
  eigenvalue of a 3×3 operator. Enumerating every Pauli choice gives exactly 2: no start, however magic, violates. The
  reason is that Pauli bases are mutually unbiased, so the transition probabilities |⟨b|a⟩|² are 0, ⅓ or 1 whatever the
  state. The two-time statistics are then a classical Markov chain.
* **Time, arbitrary bases.** Every pure start is unitarily equivalent to every other (rotate all four bases), so a
  stabilizer start already reaches the full optimum 3.1628. The start's magic is irrelevant.

**Reading.** Space can store spookiness in a state; time can store it only in the observations. Pair this with Pass 11143,
where a Bell pair is one qutrit across time up to a partial transpose. The partial transpose moves magic between "state"
and "measurement": the magic of a spatial Bell state becomes, in the temporal picture, magic of the dynamics and
measurements.
