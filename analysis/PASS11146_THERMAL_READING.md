# Pass 11146 — the neutral exit read thermally: no holonomy, no Hagedorn; the thermofield double is the Choi state of Euclidean time

Producer: `analysis/w33_pass11146_thermal_reading.py`
Scan: `analysis/w33_pass11146_scan_no_holonomy.py`
Frozen: `data/w33_pass11146_no_holonomy_21.json`
Regression: `tests/test_w33_pass11146_thermal_reading.py`

**1. No holonomy, no Hagedorn.**
* With the Wilson lines switched off, no Witten-sector tachyon appears at any radius sampled (Im T from 1.5 down to 0.2),
  in all 21 neutral-exit survivors. With them on, every model is tachyonic at small radius.
* Lattice reason: a pure winding needs l² = 1 in E8×E8 + V0, and there is no such vector with the β phase. Winding–momentum
  states have p_R² ≥ 1. The Wilson line supplies l = π + V0 + N·W with l² = 1 (Pass 11134).
* A circle with a (−1)^F-type twist is a Euclidean time circle, and a Wilson line on it is an imaginary chemical
  potential. So the neutral exit is a **Hagedorn transition that exists only because of the chemical potential**.
* Its critical curve is the Γ0(3) disk (11136). Beyond it, the T-dual description is the U(16) string (11140).

**2. The thermofield double is the Choi state of Euclidean time.** For a qutrit Gibbs state, the partial transpose of
|TFD⟩ = Σᵢ √pᵢ |ii⟩ **equals** the Jamiołkowski matrix of the Euclidean half-evolution A ↦ ρ^½ A ρ^½ (difference 0 at
every β tested). Its negativity is ((Σᵢ √pᵢ)² − 1)/2.
* At β = 0 it is 1: the temporal Bell line of Pass 11143, i.e. the identity channel.
* As β → ∞ it falls to 0: a product line.

So the spatial entanglement of the two TFD copies **is** the temporal correlation of one qutrit across imaginary time
β/2. This is the space–time swap of 11143, now for Euclidean time. With the paper's Section 3 modular clock (K = −log ρ,
Connes–Rovelli thermal time), the flow that defines time is the one whose Euclidean continuation prepares this
entanglement.

**Interpretation (stated, not computed).** Near the Hagedorn point, a winding condensate on the thermal circle is the
Horowitz–Polchinski string star, the string-scale endpoint of a black hole. Its Lorentzian continuation is an eternal black
hole whose Einstein–Rosen bridge joins the two TFD copies (ER = EPR). Read this way, the neutral condensate would be a
string-scale bridge built from entanglement across Euclidean time. The internal Wilson-line torus is not spacetime time,
so this is a formal dictionary (double Wick rotation), not a claim about our universe's clock.
