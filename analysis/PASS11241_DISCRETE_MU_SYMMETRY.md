# Pass 11241: μ under the full vacuum symmetry — impossible in Z6-I (exhaustive); never clean in Z6-II

Producer: `analysis/w33_pass11241_discrete_mu_symmetry.py`. It is an exact integer-lattice engine, also used by
Passes 11240 and 11246.
Certificate: `data/w33_pass11241_discrete_mu_symmetry.json`

**Criterion.** A superpotential term is allowed at some order iff its field vectors minus the W R-vector lie in the
integer span of:
* the condensing fields' vectors (U(1) charges, discrete charges/order, R charges/order);
* the integrality vectors.

Integer spans capture the discrete remnants of broken U(1)s. "Forbidden" verdicts are exact; "allowed" verdicts are
necessary conditions.

**Data.** Pass 10967's frozen charges:
* Z6-II: space group Z₆×Z₃×Z₂×Z₂ plus geometric R Z₆^R×Z₃^R×Z₂^R;
* Z6-I: point group Z₆ only (Pass 11247).

**Results (128 D-flat models).**
* **Z6-I, the W33 vacuum's class.**
  * Coverage: all 23 D-flat models, exhaustively over every FI ray and every choice of condensing field, **6 695 116
    vacua** in all.
  * **0 vacua** have any symmetry, continuous, remnant or point group, that forbids μ while allowing the top Yukawa.
  * In the W33 vacuum's class μ can only be forbidden by an R-symmetry, and its charges are not in the data (Pass
    11247).
* **Z6-II** (sampled: up to 150 ray supports × 4 field choices per model):
  * μ-protected vacua exist in **104 models** (25 215 protected cases), more than the 54 found with U(1)s alone,
    because the discrete rules forbid more;
  * **0** vacua are clean: in every one, some second Higgs pair or vector-like exotic stays light.

**Reading.** The scorecard's "μ must come from a vacuum symmetry" is sharpened.
* In Z6-I no available symmetry does it at all.
* In Z6-II the protecting symmetry always leaves light states. For continuous U(1)s this is Pass 11242's theorem; for
  discrete symmetries it is tested in Pass 11246.
* μ remains **OPEN**.
