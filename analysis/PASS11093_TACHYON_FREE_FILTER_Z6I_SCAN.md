# Pass 11093 — in the Z6-I W(3,3) scan, "tachyon-free" and "Standard Model" never occur together

Producer: `analysis/w33_pass11093_tachyon_free_filter_z6i_scan.py`
C++ driver: `analysis/orbifolder_n0_drivers/nsscan.cpp` (on the patched orbifolder of Pass 11092)
Inputs: `data/w33_pass11093_z6i_base_shifts.json` (the 58 base shifts), `data/w33_pass11093_scan_evidence.json`
(frozen scan logs)
Certificate: `data/w33_pass11093_tachyon_free_filter_z6i_scan.json`
Regression: `tests/test_w33_pass11093_tachyon_free_filter_z6i_scan.py`

Pass 11092 found all 87 non-supersymmetric twins of the Z6-I W(3,3) Standard Models tachyonic. It gave a
criterion: a θ fixed point of the twin is tachyonic iff the supersymmetric parent has a massless
oscillator-excited left mover there. This pass applies the criterion to the **whole** Z6-I W(3,3) scan rather
than only its Standard Models: the 58 base shifts of the Holotrade 5b3f3ad scan, with random Wilson lines
W₅ = W₆.

## A. 47 of the 58 base shifts are tachyonic whatever the Wilson lines (exact)

One of the three θ fixed points carries no Wilson line; its local shift is V itself. For **47 of the 58**
base shifts that fixed point already has ½p_sh² ∈ {7/12, 5/12, 1/4, 1/12}, so it has a tachyon. **Every model
built on those shifts is tachyonic, for any Wilson lines.** Eleven base shifts survive: 14, 23, 29, 31, 34, 35,
46, 49, 54, 55 and 56. At each of these the Wilson-line-free fixed point has minimum ½p_sh² = ¾ (for Z6I_35 no
state reaches ½p_sh² ≤ ¾).

## B. On the 11 survivors: 9,576 tachyon-free models, 705 Standard Models, and no overlap

The scan made 6000 random Wilson-line draws on each survivor (66,000 models). The C++ scanner applies the
Pass 11092 criterion to the supersymmetric parent. It builds every tachyon-free twin and runs the Standard
Model test on it. The criterion was checked on two controls: the scanner reproduces the Python verdict on the
87 known Standard-Model parents (87/87 tachyonic), and it finds the Font–Hernández model tachyon-free.

| | count |
|---|---|
| models | 66,000 |
| tachyon-free twins | **9,576** (226 inequivalent spectra) |
| supersymmetric Standard-Model parents | **705** (on 4 of the 11 shifts: 34, 46, 55, 56) |
| tachyon-free twins whose parent is a Standard Model | **0** |
| tachyon-free twins that pass the Standard-Model test themselves | **0** |

If the two properties were independent, the expected overlap would be **78.7** models, taking each base
shift's own rates. The Poisson probability of seeing 0 is **7 × 10⁻³⁵**.

92 of the 226 inequivalent tachyon-free spectra contain SU(3) and SU(2) factors, yet none passes the
Standard-Model test.

## C. Why: the Standard Models' exotics become the tachyons

For every one of the 705 Standard-Model parents, the massless oscillator-excited θ-sector states (those that
turn tachyonic in the twin) were listed. They are **always** the exotics orbifolder labels v, w, x and their
conjugates, and **never** quark, lepton or Higgs fields. No Standard-Model parent is free of them (0/705).
Pass 11094 shows that these states carry electric charge ±½ or colour. In this scan, then, the Wilson lines
that produce the Standard Model always also leave half-charged exotic states just above the massless ground
state. In the twin, those states fall below it.

## Reading

With the W(3,3) Z6-I shifts, a supersymmetry-free twin can be stable (tachyon-free) or a Standard Model, but in
66,000 draws it was never both. The anti-correlation is certified by a sample and is not proved. It is
explained by C: the Standard-Model Wilson lines in this scan always leave half-charged exotics in the
oscillator-excited θ states.

Scope. The claim covers the Z6-I twins (same shifts and Wilson lines, (−1)^F attached to θ), 6000 draws per
base shift. Other Wilson-line distributions, other twists and the SO(16)×SO(16) completion (Pass 11095) are
separate questions.
