# Pass 11125 — the B-field map of the six: charged tachyons need both B and small radius

Producer: `analysis/w33_pass11125_bfield_charged_onset.py`
Frozen map: `data/w33_pass11125_bfield_map_six.json` (Pass 11107 enumerator + Pass 11119 torus-resolved classification)
Certificate: `data/w33_pass11125_bfield_charged_onset.json`
Regression: `tests/test_w33_pass11125_bfield_charged_onset.py`

Setup: both Wilson-line tori at T = x + iy, with x (the B-field) from 0 to 0.5 and y from 0.7 to 1.4; family torus at ρ.
The map is the same in all six (N = neutral only, X = charged and neutral, `.` = none):

```
        x = 0.0 0.1 0.2 0.3 0.4 0.5
y = 1.4     N   N   N   N   N/.  N/.     ('.' in the four golden models)
y = 1.2     N   N   N   N   N   N
y = 1.1     N   N   N   N   N   N
y = 1.0     N   N   N   X   X   X
y = 0.9     N   N   X   X   X   X
y = 0.8     N   X   X   X   X   X
y = 0.7     N   X   X   X   X   X
```

Findings:
* **Charged tachyons need both a B-field and a small radius.** The highest charged cell is y = 0.8, 0.9, 1.0, 1.0, 1.0 for
  x = 0.1 … 0.5. At B = 0 the whole range down to y = 0.7 is neutral-only.
* **The neutral instability comes first on every constant-B line.** The charged onset stays at y ≤ 1.0 (to ±0.1), at
  least 0.4 below the neutral onset at y_c = 1.51 or 1.73.
* **B can lift the neutral tachyon.** At y = 1.4 and x ≥ 0.4 there is no tachyon at all in 46043, 24165, 40521 and 5904.
  These are exactly the four models with the golden critical radius φ²/√3 (Pass 11123); the two √3 models keep theirs.

Scope:
* Equal T on both Wilson-line tori; B in [0, ½], with −B related by reflection.
* Grid spacing 0.1, so boundaries are located to ±0.1.
* The one-loop potential was not evaluated off the B = 0 diagonal here. Whether B is itself driven toward 0 is Pass 11118
  territory (model 57 only) and remains open for the six.
