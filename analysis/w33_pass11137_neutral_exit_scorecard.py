#!/usr/bin/env python3
"""Pass 11137: synthesis -- the scorecard of the neutral-exit SO(16)xSO(16) W(3,3) vacua after Passes 11117-11136,
read from the certificates (no new computation; every entry points to its pass).

What WORKS (exact or class-wide):
  * the first instability met by the Wilson-line moduli is SM-neutral on the whole B circle, exactly (11136: disks at the
    Gamma_0(3) points; neutral above charged by >= 0.197 in Im T) and the potential drives the moduli into it steeply
    (11122, dLambda/dy ~ 600; B near-flat, 11129);
  * a heavy top from one localised light Higgs (11101); with that single Higgs the bottom and tau come from its conjugate
    and, WITH the condensate, sit at the same order eps^2 in 14/21 (11132): m_b/m_tau = O(1) as observed; the needed
    eps_S eps_T ~ m_b/m_t ~ 0.0135 (SM at 2e16 GeV);
  * the condensate never touches the up sector or mu (11131: 21/21).
What FAILS or stays OPEN:
  * light generations: m_c = m_u (point reflection, 11102/11127) and the conjugate sector is anarchic at every ratio of
    VEV scales (11133) -- no m_d << m_s << m_b;
  * the size of <T>: not fixed by D-flatness (no FI mechanism without SUSY, 11135; the SUSY-analogue scale 0.15-0.20
    M_P is a coincidence, not a prediction), and not by a tree-level quartic (11123);
  * the endpoint: the tachyon deepens along an exact winding law toward the Delta = -1/2 state of a T-dual string with no
    tachyon-free island (11128, 11134) -- the condensation leaves perturbative control;
  * the vacuum energy is positive at one loop (11106).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
OUT = D / "w33_pass11137_neutral_exit_scorecard.json"


def load(name):
    return json.loads((D / name).read_text())


def summarize():
    dom = load("w33_pass11136_tachyon_domains.json")
    bt = load("w33_pass11132_bottom_tau_across_21.json")
    un = load("w33_pass11131_unlocked_across_21.json")
    fi = load("w33_pass11135_fi_scale_no_mechanism.json")
    pc = load("w33_pass11128_pinched_cycle.json")
    sp = load("w33_pass11134_dual_tachyon_species.json")
    tx = load("w33_pass11127_lepton_texture_under_condensation.json")
    works = dict(neutral_first_exactly=dom['min_gap_neutral_over_charged'] > 0,
                 bottom_tau_same_order_models=len(bt['lowered_models']),
                 up_mu_never_unlocked=un['models_with_unlocked']['up'] == 0 and un['models_with_unlocked']['mu_HuHd'] == 0)
    fails = dict(m_c_equals_m_u_and_m_mu_equals_m_e=tx['light_lepton_reflection_symmetric'],
                 conjugate_sector_anarchic=bt['down_anarchic_every_r'] and bt['three_distinct_exponent_cases'] == 0,
                 condensate_size_unfixed=not fi['fi_mechanism_applies'],
                 endpoint_nonperturbative=pc['neutral_only_down_to_0p2'] and pc['deepest'] < -0.44,
                 dual_tachyon_one_species=sp['one_complex_species'])
    res = dict(pass_id=11137, works=works, fails_or_open=fails, eps_product_needed=0.0135,
               all_consistent=all(works.values()) and all(fails.values()))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
