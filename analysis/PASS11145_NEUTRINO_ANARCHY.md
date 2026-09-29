# Pass 11145 — neutrino anarchy from the condensate: tested and refuted

Producer: `analysis/w33_pass11145_neutrino_anarchy.py`
Scan: `analysis/w33_pass11145_scan_neutrino_dirac.py`
Frozen: `data/w33_pass11145_neutrino_textures.json`
Regression: `tests/test_w33_pass11145_neutrino_anarchy.py`

**Hypothesis (a new direction).** The conjugate sector is anarchic (11133, 11138), and anarchy is a known explanation of
large lepton mixing (Hall–Murayama–Weiner). If the condensate also made the neutrino Dirac matrix anarchic, the class
would predict large PMNS angles alongside hierarchical quarks.

**Test.** L·N^c·h for the three tree-level top doublets h and all 39–51 singlet fermions N, with and without the
condensate, in 36621, 24165 and 40521. The charged-lepton texture comes from conj(h) (Pass 11142).

**Result: refuted.**
* **The neutrino Dirac texture is [0, 1, 1] at tree level** for every light Higgs, the up-quark pattern: one top-like
  Dirac mass and a degenerate pair.
  * It is identical with and without the condensate: the unlocked entries of 11124/11131 lie below existing tree-level
    ones.
  * This is an SO(10)-like relation, m_D(ν) ~ m_u in texture.
* **The charged leptons are anarchic** ([2,2,2] or [3,3,3]), so lepton mixing would be large, as observed.
* **But the down quarks come from the same conjugate sector,** so quark mixing would be large too. That conflicts with
  the small CKM angles.
* A top-like Dirac neutrino mass needs a seesaw scale M_R ~ m_t²/0.05 eV ≈ 6·10¹⁴ GeV. The Majorana sector is not
  computed here.
