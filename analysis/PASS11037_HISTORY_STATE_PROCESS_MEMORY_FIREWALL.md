# Pass 11037 — relational cubic history is not the same resource as quantum temporal memory

Producer: `analysis/w33_pass11037_history_state_process_memory_firewall.py`
Certificate: `data/w33_pass11037_history_state_process_memory_firewall.json`
Regression: `tests/test_w33_pass11037_history_state_process_memory_firewall.py`

For the CCZ history state built from |+>³, the cubic phase has eigenspace counts 19,4,4. The conditional-history overlap is exactly 5/9, so the reduced qutrit clock spectrum is

(19/27, 4/27, 4/27),

with entropy 1.1730125853 bits.

Choosing equal amplitude from the three CCZ eigensectors makes the three orbit histories mutually orthogonal. The clock is then maximally mixed with entropy log₂3.

Those facts certify relational history encoding. They do not certify quantum memory across interventions.

We therefore implement an operational causal-break witness. A trace-and-prepare break with no hidden memory destroys reference entanglement and has negativity zero. A hidden qutrit store/retrieve channel preserves the maximally entangled qutrit pair and has negativity one. A computational-basis measure-and-prepare memory also has negativity zero.

For a depolarizing qutrit memory channel,

N(η) = max(0,(4η−1)/3),

so quantum memory survives exactly for η>1/4. Under pure coherence damping the Bell-pair negativity equals the remaining coherence λ.

This separates three resources cleanly: clock-system entanglement = relational history; the cubic echo = temporal holonomy; reference entanglement surviving a causal break = quantum temporal memory.
