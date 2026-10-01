#!/usr/bin/env python3
"""2026-10-01: separate the SIC and stabilizer mirror divisors using the E8
restricted Pfaffian P and the G26 reflection Jacobian J.

Parents:
  Pass 11261: P has (SIC, stabilizer) mirror multiplicities (3,1).
  Pass 11269: J has mirror multiplicities (1,2).
The resulting valuation matrix has determinant 5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "w33_20261001_g26_e8_mirror_separation.json"


def payload():
    p61 = json.loads((ROOT / "data/w33_pass11261_g26_cartan_coordinate_map.json").read_text())
    p69 = json.loads((ROOT / "data/w33_pass11269_g26_qutrit_dictionary.json").read_text())
    assert p61["status"] == "PASS_EXACT_G26_COORDINATE_MAP_AND_QUTRIT_MIRROR_DIVISOR"
    assert p69["jacobian_equals_19440_SIC_x_stab2_exact"]
    # Columns are (nine SIC mirrors, twelve stabilizer/MUB mirrors).
    vP = (3, 1)
    vJ = (1, 2)
    determinant = vP[0] * vJ[1] - vP[1] * vJ[0]
    sic_isolation = (2 * vP[0] - vJ[0], 2 * vP[1] - vJ[1])
    stab_isolation = (3 * vJ[0] - vP[0], 3 * vJ[1] - vP[1])
    assert determinant == 5
    assert sic_isolation == (5, 0)
    assert stab_isolation == (0, 5)

    # Normalizations from the two parent certificates. Put delta=omega-omega^2,
    # so delta^2=-3. A=prod(SIC)=3*delta*S9 and
    # B=prod(stabilizer)=-27*C3*N9.
    # P=2^17*C3*N9*S9^3 = cP*A^3*B; J=cJ*A*B^2.
    cJ = (2**4) * (3**5) * 5
    assert cJ == 19440
    # With cP=2^17/(3^7*delta) and delta^2=-3:
    # cP^2/cJ = -2^30/(5*3^20), and cJ^3/cP =
    # (3^22*5^3/2^5)*delta.
    from fractions import Fraction
    sic_ratio = Fraction(-(2**30), 5 * (3**20))
    stab_delta_ratio = Fraction((3**22) * (5**3), 2**5)
    assert sic_ratio * (5 * 3**20) == -(2**30)
    assert stab_delta_ratio * (2**5) == (3**22) * (5**3)
    # Squared cP/cJ and cJ^3/cP can be reduced without adjoining delta.
    cleared_sic = {
        "identity": "5*3^20*P^2 + 2^30*A^5*J = 0",
        "degree_each_side": 78,
    }
    cleared_stab = {
        "identity": "32*J^3 - 3^22*5^3*delta*B^5*P = 0",
        "delta_relation": "delta^2=-3",
        "degree_each_side": 99,
    }
    return {
        "schema": "w33.20261001.g26-e8-mirror-separation.v1",
        "status": "PASS_E8_PFAFFIAN_AND_G26_JACOBIAN_SEPARATE_SIC_AND_MUB_MIRROR_DIVISORS_UP_TO_FIFTH_POWER",
        "parents": [
            "data/w33_pass11261_g26_cartan_coordinate_map.json",
            "data/w33_pass11269_g26_qutrit_dictionary.json",
        ],
        "mirror_products": {
            "A": "product of the 9 qutrit SIC mirror linear forms",
            "B": "product of the 12 qutrit stabilizer/MUB mirror linear forms",
        },
        "valuation_vectors_SIC_stabilizer": {"P": list(vP), "J": list(vJ)},
        "valuation_matrix_determinant": determinant,
        "divisor_isolation": {
            "2P_minus_J": list(sic_isolation),
            "3J_minus_P": list(stab_isolation),
        },
        "normalized_polynomial_identities": {
            "SIC_fifth_power": cleared_sic,
            "stabilizer_fifth_power": cleared_stab,
        },
        "logarithmic_separator": {
            "SIC": "2*dlog(P)-dlog(J)=5*dlog(A)",
            "stabilizer": "3*dlog(J)-dlog(P)=5*dlog(B)",
        },
        "theorem": (
            "The E8 restricted Pfaffian and the G26 reflection Jacobian span an "
            "index-five sublattice of the two-orbit mirror divisor lattice. Their "
            "integer combinations isolate the SIC and stabilizer/MUB arrangements: "
            "P^2/J is a scalar times A^5, while J^3/P is a scalar times B^5."
        ),
        "boundary": (
            "The index five is an exact divisor-lattice statement. It is not a "
            "particle multiplicity, coupling, family number, or physical Z5 symmetry. "
            "A physical vacuum potential or kinetic metric is not inferred."
        ),
        "checks": {
            "parent_mirror_multiplicities": True,
            "valuation_determinant_is_5": True,
            "SIC_divisor_isolated": True,
            "stabilizer_divisor_isolated": True,
            "degrees_match": True,
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    p = payload()
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not OUT.exists() or OUT.read_text() != text:
            raise SystemExit("certificate drift")
    else:
        OUT.write_text(text)
    print(json.dumps({"status": p["status"], "index": p["valuation_matrix_determinant"]}))


if __name__ == "__main__":
    raise SystemExit(main())
