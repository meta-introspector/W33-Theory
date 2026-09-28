# Pass 11109 — at the stabilised family modulus the Yukawas give m_c = m_u = m_t/2

Producer: `analysis/w33_pass11109_instanton_yukawas_at_rho.py`
Certificate: `data/w33_pass11109_instanton_yukawas_at_rho.json`
Regression: `tests/test_w33_pass11109_instanton_yukawas_at_rho.py`

## Setup

In the 12 survivors the three generations sit at the three θ-fixed points of the Wilson-line-free family torus, and the
light Higgs sits at one of them. Tree-level twisted Yukawa couplings factorise over the tori. On the family torus there
are two cases:
* all three fields at the same point, which gives the top: coupling Y_same;
* three distinct points, which gives the off-diagonal charm–up block: coupling Y_dist.

The other tori contribute a common factor, and the Kähler normalisations cancel in the ratio. Hence, at tree level,

    m_c = m_u,     m_{c,u}/m_t = |Y_dist(T*)/Y_same(T*)|.

## The couplings

The worldsheet instanton maps the sphere onto **two** copies of the triangle spanned by the three fixed points
(Schwarz–Christoffel). Its action is therefore 2 × area/(2πα′). In the Narain normalisation of the one-loop engine (SU(3)
point at T = ρ, torus area = 4π²α′ Im T), the couplings are the A2 theta functions of the Kähler modulus:

    Y_c(T) = Σ_{v ∈ A2 + c} e^{iπT|v|²}     (c = 0 same point; c = f, 2f distinct points)

This is the known structure: Lauer, Mas and Nilles 1989–91; Kobayashi et al., arXiv:1804.06644.

**Check.** The triple (Y_0, Y_f, Y_2f) closes under T → −1/T with weight (−iT) and the finite Fourier matrix
(1/√3)[ω^{jk}]. The single-triangle normalisation does not close, so a first estimate that used it is discarded. The
minimal distinct-point triangle is 1/6 of the cell (checked).

## Results

| family modulus T\* | \|Y_dist / Y_same\| = m_{c,u}/m_t (tree) |
|---|---|
| **ρ, the one-loop minimum in this direction (Pass 11106)** | **1/2 exactly** |
| i | (√3 − 1)/2 = 0.366 |
| 1.5i | 0.130 |
| 2i | 0.045 |
| 3i | 0.0056 |
| 4i | 0.00069 |

The observed m_c/m_t ≈ 0.0036 (at M_Z) needs **Im T\* ≈ 3.2**; 0.0027 at a high scale needs 3.35.

## Reading

The two results the survivors need are incompatible:
* the one-loop potential stabilises the family modulus at its self-dual point ρ (Pass 11106);
* the single heavy top of Pass 11101 needs a large family torus, Im T\* ≈ 3.2.

At the stabilised point the top is only twice as heavy as charm and up. Pass 11110 checks that the potential at
Im T\* = 3.2 lies above its value at ρ. V_cb is moot while the tree-level 2–3 ratio is ½.

Scope: tree-level magnitudes. Higher-order terms from scalar VEVs are O(ε_VEV²) corrections and do not change a ratio of
order ½. The exact value ½ at ρ and (√3−1)/2 at i are properties of the A2 theta functions at the fixed points of the
duality group.
