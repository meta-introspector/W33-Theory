# 2026-10-01 — the G26 magic vacuum needs a frame to compute

Producer: `analysis/w33_20261001_tmagic_orbit_factory_contract.py`
Certificate: `data/w33_20261001_tmagic_orbit_factory_contract.json`
Regression: `tests/test_w33_20261001_tmagic_orbit_factory_contract.py`

## What is added to the existing results

The balanced mirror calculation owns the exact T-magic local minimum and its
72-ray Clifford orbit. Pass 411 and `w33_qutrit_t_teleportation_port.py` already
own T injection and its nine Clifford corrections. This pass connects those
objects and specifies the missing classical resource.

The 72 rays partition into **eight Pauli orbits of nine**, with exactly two
Pauli orbits over each of the four MUB axes. A known axis does not determine
a usable magic state. Even a known nine-ray orbit, if its displacement label
is uniformly forgotten, gives exactly the maximally mixed state:

    (1/9) sum_(a,b) X^a Z^b rho Z^-b X^-a = Tr(rho) I/3.

The identity is proved symbolically for an arbitrary 3-by-3 input operator.
Pauli twirling is standard; its consequence for this new geometric vacuum
family is the resource audit established here.

## Retained frame: the existing universal-computation interface works

For each ray retain a projective Clifford representative C with ray = C|T3>.
Apply C dagger, then SUM to the state and a fresh |0>:

    SUM (|T3> tensor |0>) = (I tensor T)|Phi3>.

This is the exact Choi resource consumed by the existing T port. All **648**
ray/Bell branches (72 times 9) satisfy the full operator identity

    correction_(a,b) * Kraus_(a,b) = T/3.

Each outcome has probability 1/9 independent of the input; because these are
operator identities they also apply to an input entangled with another
register. T plus the existing Clifford/SUM operations supplies the known
approximately universal qutrit gate set. This is a connection to prior
injection, not a new universality theorem.

## Forgotten frame: the apparent factory loses its quantum resource

SUM applied to I/3 tensor |0><0| instead creates

    (1/3) sum_j |jj><jj|,

which is a separable resource. Feeding it to the same Bell/correction port
produces **rho -> diag(rho)**, an entanglement-breaking dephasing channel.
The certificate replays this on all nine input matrix units, including the
six off-diagonal ones. Thus the same geometric candidate can implement a T
gate or erase coherence depending on whether its relative frame is available.

The claim concerns uniform forgotten-frame ensembles. A single fixed but
unknown frame shared across a batch may be learnable; correlated resources
and reference alignment are separate questions. Scalar G26 invariants and a
MUB-axis label alone cannot replace that alignment information.

## Physical boundary

The barrier proves an exact local stationary orbit. It does not provide a
state-preparation Hamiltonian, cooling mechanism, measured noise rate or
fault-tolerant factory. The group/orbit census is numerical; the universal
Pauli-twirl identity is symbolic. A physical implementation must both prepare
the resource and retain or establish its frame.

External check: [Glaudell et al., Qutrit Metaplectic Gates Are a Subset of
Clifford+T](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.TQC.2022.12)
records the standard qutrit T gate's injectable, approximately universal role.
