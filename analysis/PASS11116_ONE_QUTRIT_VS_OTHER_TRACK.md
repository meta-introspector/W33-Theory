# Pass 11116 — the one-qutrit theorem against the other track's W(3,3)/H27/243 constructions

Producer: `analysis/w33_pass11116_one_qutrit_vs_other_track.py` (`<dump>` recomputes both tables)
Data: `data/w33_pass11116_twist_is_gauge_104.json`, `data/w33_pass11116_fi_and_hidden_gluing_104.json`
Certificate: `data/w33_pass11116_one_qutrit_vs_other_track.json`
Regression: `tests/test_w33_pass11116_one_qutrit_vs_other_track.py`

The other track built three finite objects inside E8, all stated with "no heterotic vacuum inferred":
* **(a)** the external-A2 qutrit H27 on the family index of E8 ⊃ E6 × SU(3), with centre = the FI Z3
  (`2026-09-21_physical_external_a2_h27.md`, `..._physical_fi_is_h27_center.md`);
* **(b)** its Clifford normaliser H27:SL(2,3) of order 648, which is the W33 point stabiliser (1296 with similitudes)
  (`..._physical_a2_clifford648_w33_bridge.md`, `2026-09-23_extended_clifford1296_point_stabilizer.md`);
* **(c)** the two-qutrit group 3^(1+4) = H27_int ∘ H27_ext, with the internal H27 in the E6 trinification and centre
  Z(E6). Its commutation form is W(3,3), and on (27, 3) it acts as 9 copies of its 9-dimensional irrep
  (`..._e8_trinification_two_qutrit_pauli243.md`, `..._e8_pauli243_projective_w33_bridge.md`,
  `..._e8_matter81_pauli243_restriction.md`).

## Tests on the 104 SO(16)×SO(16) A8 vacua

| test | result |
|---|---|
| is the orbifold twist (= the centre of the family H27) a **gauge element** on the massless spectrum? (lattice test: the charge lattice meets the l-axis only in integers) | **104/104** with U(1)s and non-abelian centres; 97/104 with U(1)s alone |
| is it exactly the anomalous-U(1) Z3 (the other track's "FI = H27 centre")? | 3/104 (models 4, 25, 29); no survivor |
| is there an SU(3) factor, colour or hidden, whose triality equals ±l on every field it acts on (the centre gluing the 243 group needs)? | **0/104** |

The first line means the Heisenberg centre acts on the massless spectrum as a gauge transformation, so the point-group
Z3 is redundant with gauge invariance there. The Δ(54) of Z3 orbifolds has a known gauge origin at the enhanced
point (Beye, Kobayashi, Kuwakino, arXiv:1401.2310).

## Classification

* **(a) Physical.** It is realised as the family Δ(54) of every model (Pass 11105). Its centre is the orbifold twist,
  which is always a gauge element and, in 3 models, the FI Z3 itself.
* **(b) Only ±I realised,** index 12 (Pass 11105).
* **(c) Not realised on any matter:**
  * on the generations alone, excluded by the faithful degree 3^k (Pass 11111);
  * on generation ⊗ colour, excluded because the twist is not the colour centre (Pass 11111);
  * on hidden SU(3) matter, 0/104 (here).
  Its E8 construction needs the E6 trinification centre, which these vacua lack, and no substitute exists.
