# Pass 11264 — restored Z6-I plane R charges reverse the mu-protection verdict

Producer: `analysis/w33_20261001_z6i_rcharge_lattice_quotient.py`
Raw extractor/reference: `analysis/w33_pass11264_z6i_rcharge_mu_audit.py`
Certificate: `data/w33_20261001_z6i_rcharge_lattice_quotient.json`

## Result

The missing Z6-I plane R charges were restored from the Orbifolder 1.2.1 model
archives using the same extraction formula validated against the frozen Z6-II
snapshot: `R_i = q_sh,i + N_i - Nbar_i`, with plane orders `(6,6,3)`.

The raw archive contains all **87** Z6-I models and **30,980** fields.
D-flatness depends only on the continuous FI cone, so the exact Pass 11241
classification remains 23 D-flat models and 64 non-D-flat models. Coupling
selection was then replayed with the restored R coordinates.

The exhaustive result is **6,695,116 / 6,695,116** D-flat vacua protecting
`mu = H_u H_d` while allowing at least one top Yukawa. All 23 D-flat models
contain protected vacua; in fact every enumerated vacuum is protected.

This **supersedes the Z6-I mu-protection reading of Pass 11241**, whose exact
integer lattice did not contain the absent plane R charges. It does not
invalidate that pass's D-flat enumeration.

## Why the replay is fast and still exact
The 6.7 million vacua are mostly field-choice multiplicities over only 22--108
FI-ray supports per model. For a fixed support, two field choices that generate
the same exact integer-echelon lattice make exactly the same all-orders
coupling decisions.

The new engine therefore quotients choices by the frozen `Lattice.ech`
signature, evaluates Higgs/top/exotic membership once per distinct lattice,
and restores the multiplicity afterwards. Across the entire class the
6,695,116 vacua collapse to only **1,464** support/lattice states.

This quotient was checked against the original per-vacuum engine on 20 explicit
choices in every D-flat model. Protected, clean, and clean-plus-F-flat verdicts
agree choice by choice.
## The vacuum problem is not solved

The clean count is still **0 / 6,695,116**.

So the restored R symmetry closes the mu loophole but simultaneously preserves
additional unwanted light structure. The companion exhaustive colored-mass
audit strengthens this: none of the 6,695,116 protected-Hu cases admits a
complete colored matching.

## Boundary

“Protected” here means forbidden to all orders by the exact frozen
U(1)-remnant + space-group + plane-R charge lattice while a top Yukawa remains
allowed. It does not compute the eventual mu scale after R breaking and does
not establish a phenomenologically viable vacuum.
