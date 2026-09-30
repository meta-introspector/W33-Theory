# Pass 11205 — the one-way clock: time-oriented relations between splits have Eisenstein-complex spectra

Producer: `analysis/w33_pass11205_one_way_clock.py`
Certificate: `data/w33_pass11205_one_way_clock.json`
Regression: `tests/test_w33_pass11205_one_way_clock.py`

**Idea.** The 110 565 three-qutrit splits carry a rank-20 coherent configuration.
* Six of its relations are directed, in reversal pairs of sizes 256, 2304 and 6912 (Pass 11189).
* A process that repeatedly re-splits the register along a directed relation R never steps back along R.
* Its spectrum is that of the adjacency matrix A_R. The matrix L_R of multiplication by A_R, whose entries are
  intersection numbers computed from the orbit labels, has the same eigenvalues.

**Result.**
* **Every directed relation has non-real eigenvalues, and its time reverse has exactly the complex-conjugate
  spectrum.** The symmetric control relation (size 36) has a real spectrum.
* The rotation frequencies lie in the **Eisenstein field Q(√−3)**, the field of the qutrit's cube roots of unity:
  * 256 relation: −14 ± 18√3·i and −8 ± 12√3·i;
  * 6912 relation: −54 ± 162√3·i and ±36√3·i;
  * 2304 relation: −30 ± 6√591·i, with multiplicity 2 (a non-commutative block).
* L₃ commutes with its reverse L₄ but not with L₁₂, so the coherent configuration is **non-commutative**.

**Reading.** A random walk that re-chooses the subsystem split one directed step at a time is a clock that runs one
way. Its non-real spectrum makes probability circulate, like a Markov chain without detailed balance, and the reverse
relation circulates the other way. The imaginary parts of the 256 and 6912 relations are integer multiples of √3, the
frequency scale of the qutrit phase e^{2πi/3}.

**Scope.** This is spectral data of a finite relation algebra. No physical time scale or Hamiltonian is claimed.
