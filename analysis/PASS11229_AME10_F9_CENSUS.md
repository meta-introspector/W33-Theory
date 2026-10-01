# Pass 11229: an independent, larger AME(10,3) census — every state found is Glynn's

Producer: `analysis/w33_pass11229_ame10_f9_census.py`. It reruns master's tabu search (Pass 11186) with fresh seeds
11229000–11229002 and applies Pass 11224's tests.
Certificate: `data/w33_pass11229_ame10_f9_census.json`
Regression: `tests/test_w33_pass11229_ame10_f9_census.py`

**Result.**
* 487 restarts were run, about 82 minutes on 3 workers.
* **138** AME(10,3) graph states were found. All 138 are new relative to master's census of 71.
* **Every one is F₉-linear**: 4 F₉-structures each, and a Schur square of dimension 10, which identifies Glynn's
  non-classical arc as in Pass 11224.
* Together with master's 71 states, that makes 209 independent hits and no non-F₉-linear state.

**Status.** This is evidence, not proof, as with any negative random search. Pass 11249 adds an exhaustive algebraic
result for states with a ten-cycle symmetry. Whether a non-F₉-linear AME(10,3) stabilizer state exists remains open.
