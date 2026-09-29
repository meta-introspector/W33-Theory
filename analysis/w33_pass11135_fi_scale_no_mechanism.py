#!/usr/bin/env python3
"""Pass 11135: could the anomalous-U(1) FI term fix the condensate size?  The would-be value lands in the needed range
-- but there is no FI mechanism for a complex scalar without supersymmetry, so this is a coincidence of scales, not a
derivation.

Data: data/w33_pass11135_fi_data_21.json (scan analysis/w33_pass11135_scan_fi.py): for each of the 21 neutral-exit
survivors, the single anomalous U(1) of the U1STD basis, Tr Q_A over left-handed fermions, |u_A|^2, and the anomalous
charge q_{T,A} of the neutral tachyon.
Results:
  * every model has exactly one anomalous U(1); Tr Q_A = 12 |u_A|^2 in 21/21 -- a TAUTOLOGY of the orbifolder's
    definition t_A = (1/12) sum_f p_f, recorded as a consistency check, not a finding;
  * the tachyon carries anomalous charge |q_{T,A}| = 3 or 4 (Pass 11117's "component along the anomalous U(1)");
  * the SUPERSYMMETRIC formula <phi>^2 = g^2 Tr Q_A M_P^2 / (192 pi^2 |q|) would give <T>/M_P = 0.146-0.203 at g^2 = 1/2
    -- inside the eps ~ 0.12-0.15 that Pass 11130 needs for m_b/m_t;
  * but that formula needs a D-term: the FI mass is g^2 q xi |phi|^2 for the scalar of a CHIRAL multiplet.  For a complex
    scalar without supersymmetry nothing distinguishes phi from phi* (the dumps list both: Hd_j = conj(Hu_k) exactly,
    Pass 11130), CPT gives particle and antiparticle the same mass, and a term linear in q is not covariant under
    phi <-> phi*.  So no FI potential fixes <T> here.
Reading: the condensate size is not fixed by D-flatness in the non-supersymmetric vacuum; the scale coincidence
(0.15-0.20 vs 0.12-0.15) is recorded, and explicitly NOT claimed as a prediction.  Convention note: the SUSY value
changes by sqrt2 between the two standard normalisations of the FI term.
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "w33_pass11135_fi_data_21.json"
OUT = ROOT / "data" / "w33_pass11135_fi_scale_no_mechanism.json"


def summarize():
    d = json.loads(DATA.read_text())
    vals = [v['T_over_MP_g2half'] for v in d.values()]
    res = dict(pass_id=11135, n_models=len(d),
               one_anomalous_u1=all(v['anom_index'] == 0 for v in d.values()),
               trq_equals_12_norm2=all(F(v['TrQA']) == 12 * F(v['uA_norm2']) for v in d.values()),
               qTA_abs=sorted({abs(int(v['qT_A'])) for v in d.values()}),
               susy_formula_T_over_MP=[min(vals), max(vals)],
               needed_eps=[0.116, 0.151], fi_mechanism_applies=False)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
