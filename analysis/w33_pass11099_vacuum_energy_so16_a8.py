#!/usr/bin/env python3
"""Pass 11099: the one-loop vacuum energy of the SO(16)xSO(16) A8 models (Pass 11095) is positive and runs away.

A. The 10D SO(16)xSO(16) modular integral (w33_pass11099_lambda10_so16xso16.py): the integrand is checked modular
   invariant; its level-matched constant term is -2112 = n_B - n_F of the 10D massless spectrum (1984 - 4096); the
   integral over the fundamental domain is I = -726, so Lambda_10 = -(1/2) M^10 I = +363 M^10 > 0 (M = M_s/2pi):
   POSITIVE, as Alvarez-Gaume, Ginsparg, Moore and Vafa found.
B. Large volume.  The orbifold group is G = Z2W x Z3 (|G| = 6); Z = (1/6) sum_{g,h in G} Z(g,h).  Only the sectors with
   g, h in {1, beta} (beta = Witten element, trivial on T6) carry the T6 zero modes, and
       (1/6) sum_{g,h in {1,beta}} Z(g,h) = (1/3) Z[SO(16)xSO(16) on T6].
   Every other sector rotates all three planes (no zero modes: moduli-independent).  Hence at large T6 volume V (string
   units) Lambda_4 = (V/3) Lambda_10 + O(V^0) > 0 for every one of the 104 models: the potential grows with the
   volume, and the dilaton tadpole (proportional to Lambda) drives the dilaton to weak coupling.  No stabilised vacuum
   at large volume; the value at the orbifold point needs the full twisted-sector computation (not done here).
C. Massless Bose-Fermi counts of the 104 models (frozen from our spectra): n_B - n_F in [-500, -64], never 0.  The
   exponentially suppressed Lambda of interpolating models (Itoyama-Taylor; Abel-Dienes-Mavroudi) requires massless
   Bose-Fermi degeneracy; none of the 104 has it.  And by Groot Nibbelink et al. (arXiv:1710.09237) no symmetric
   non-supersymmetric toroidal orbifold has a vanishing one-loop Lambda by sector-wise Killing spinors.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11099_lambda10_so16xso16 as L  # noqa: E402

BF = ROOT / "data" / "w33_pass11099_massless_bose_fermi.json"
OUT = ROOT / "data" / "w33_pass11099_vacuum_energy_so16_a8.json"


def volume_sector_weight():
    """G = Z2 (beta, trivial on space) x Z3 (theta, rotates all three planes).  Fraction of (1/|G|) sum Z(g,h) carried by
    sectors with zero modes, expressed as a multiple of Z[SO16xSO16 on T6] = (1/2) sum_{g,h in Z2} Z(g,h)."""
    G = [(a, k) for a in range(2) for k in range(3)]
    vol = [(g, h) for g in G for h in G if g[1] == 0 and h[1] == 0]      # both geometrically trivial
    from fractions import Fraction
    weight = Fraction(len(vol), len(G)) / Fraction(4, 2)                  # (n_vol/|G|) / (4 sectors / |Z2|)
    rot = all(g[1] != 0 or h[1] != 0 for g in G for h in G if (g, h) not in vol)
    return dict(group_order=len(G), sectors=len(G) ** 2, zero_mode_sectors=len(vol),
                weight_of_so16_on_T6=str(weight), other_sectors_rotate_all_planes=rot)


def main():
    import numpy as np
    c = L.exact_level_matched()
    tau = complex(0.21, 1.13)
    inv = abs(L.Z(-1 / tau) - L.Z(tau)) / abs(L.Z(tau))
    lam = json.loads(L.OUT.read_text()) if L.OUT.exists() else None
    if lam is None:
        L.main()
        lam = json.loads(L.OUT.read_text())
    bf = json.loads(BF.read_text())
    vals = [v["nB_minus_nF"] for v in bf.values()]
    res = dict(pass_id=11099,
               lambda10=dict(I=lam["I"], Lambda10_over_M10=lam["Lambda10_over_M10"], sign=lam["sign"],
                             massless_constant_term=c.get(0), expected_nB_minus_nF=8 * (8 + 240) - (8 * 256 + 8 * 256),
                             modular_S_relative_error=float(inv)),
               large_volume=volume_sector_weight(),
               massless_bose_fermi=dict(models=len(vals), nB_minus_nF_min=min(vals), nB_minus_nF_max=max(vals),
                                        degenerate=sum(1 for x in vals if x == 0)))
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
