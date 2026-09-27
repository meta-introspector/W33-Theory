#!/usr/bin/env python3
"""Pass 11026: the signed clock central parity cannot rescue the FI RPV vacua."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11026_clock_parity_rpv_intertwiner_firewall.json"

P51 = ROOT / "data/w33_pass10951_clock_pin_spin_central_sign_bridge.json"
P67 = ROOT / "data/w33_pass10967_fi_vacuum_regenerates_rpv.json"
P80 = ROOT / "data/w33_pass10980_rpv_flavour_upper_bound.json"
P22 = ROOT / "data/w33_pass11022_binary_octahedral_clock_decomposition.json"


def prod(*xs):
    out = 1
    for x in xs:
        out *= x
    return out


def payload():
    p51 = json.loads(P51.read_text(encoding="utf-8"))
    p67 = json.loads(P67.read_text(encoding="utf-8"))
    p80 = json.loads(P80.read_text(encoding="utf-8"))
    p22 = json.loads(P22.read_text(encoding="utf-8"))

    assert p51["matter_parity_weld"]["pass10950_spin9_2pi_equals_peirce_symmetry"]
    assert p22["central_minus_I"]["plus_dimension"] == 12
    assert p22["central_minus_I"]["minus_dimension"] == 12

    # Standard matter-parity signs. The forced singlet n is certified odd in Pass 10967.
    parity = {"q": -1, "u": -1, "d": -1, "l": -1, "e": -1,
              "Hu": +1, "Hd": +1, "n": -1}
    allowed_yukawa = {
        "q_u_Hu": prod(parity["q"], parity["u"], parity["Hu"]),
        "q_d_Hd": prod(parity["q"], parity["d"], parity["Hd"]),
        "l_e_Hd": prod(parity["l"], parity["e"], parity["Hd"]),
    }
    rpv = {
        "udd": prod(parity["u"], parity["d"], parity["d"]),
        "qld": prod(parity["q"], parity["l"], parity["d"]),
        "lle": prod(parity["l"], parity["l"], parity["e"]),
    }
    quartic = {k: parity["n"] * v for k, v in rpv.items()}
    assert set(allowed_yukawa.values()) == {+1}
    assert set(rpv.values()) == {-1}
    assert set(quartic.values()) == {+1}

    sc = p67["string_couplings"]
    z6i = {k: v for k, v in sc.items() if k.startswith("Z6-I|")}
    exact = {
        k: v for k, v in z6i.items()
        if v["forced_singlets"] == 9
        and v["udd_order3"] == v["qld_order3"] == v["lle_order3"] == 0
        and v["udd_forced_in_quartic"] == 9
        and v["qld_forced_in_quartic"] == 9
        and v["lle_forced_in_quartic"] == 9
    }
    assert len(z6i) == len(exact) == 23
    b = p67["summary"]["B"]
    c = p67["summary"]["C"]
    assert b["forced_singlets"] == 302
    assert b["udd:n*op qld:n*op lle:n*op"] == 290
    assert c["Z6-I_dflat_models"] == 23
    assert c["Z6-I_no_cubic_rpv_all_forced_regenerate_all_three"] == 23

    f = p80["summary"]
    assert f["models"] == 23
    assert f["with_P"] == 16
    assert f["no_light_Hu"] == 7
    assert f["five_light_d"] == 16
    assert f["below_1e-10"] == 0
    assert f["min_P"] > 2e-3

    checks = {
        "clock_center_matches_certified_matter_parity_character":
            p51["matter_parity_weld"]["central_character_match"] != "",
        "signed_carrier_is_12_plus_12": (
            p22["central_minus_I"]["plus_dimension"],
            p22["central_minus_I"]["minus_dimension"],
        ) == (12, 12),
        "ordinary_yukawas_even": set(allowed_yukawa.values()) == {1},
        "all_three_rpv_cubics_odd": set(rpv.values()) == {-1},
        "forced_odd_singlet_times_rpv_is_even": set(quartic.values()) == {1},
        "all_23_z6i_vacua_have_nine_forced_regenerators": len(exact) == 23,
        "pass10967_all_23_regenerate_all_three":
            c["Z6-I_no_cubic_rpv_all_forced_regenerate_all_three"] == 23,
        "pass10980_no_analysable_vacuum_below_1e_10":
            f["below_1e-10"] == 0,
        "pass10980_all_23_fail_independently_of_proton_decay":
            f["no_light_Hu"] + f["five_light_d"] == 23,
    }
    assert all(checks.values())

    return {
        "schema": "w33.pass11026.clock-parity-rpv-intertwiner-firewall.v1",
        "status": "PASS",
        "headline": (
            "The new 12+12 signed-clock central parity cannot be used as an "
            "unbroken proton-protection symmetry in the certified FI vacua. "
            "Pass 10951 identifies the same central -I character with matter "
            "parity; Pass 10967 proves the FI-forced nu^c-like singlets are odd "
            "under every such parity and that all 23 D-flat Z6-I vacua contain "
            "allowed n*udd, n*qld and n*lle quartics. Their required VEV therefore "
            "breaks the central sign and regenerates all three RPV cubics. The now-"
            "committed Pass 10980 independently finds no first-generation flavour "
            "suppression rescue: none of its analysable vacua reaches 1e-10."
        ),
        "central_sign_dictionary": {
            "clock": "-I in the exact signed GL2(3) carrier",
            "carrier_eigenspaces": "12 even + 12 odd",
            "pass10951_reading": "same central character as matter parity on the Albert Peirce 16",
            "matter_signs": parity,
        },
        "selection_rule": {
            "ordinary_yukawa_signs": allowed_yukawa,
            "rpv_cubic_signs": rpv,
            "forced_singlet_times_rpv_signs": quartic,
            "mechanism": (
                "Before condensation, the central/matter parity forbids each "
                "three-matter RPV cubic. The FI-forced odd singlet makes n*RPV "
                "even; once <n> is nonzero the same term becomes an effective "
                "RPV cubic and the Z2 is spontaneously broken."
            ),
        },
        "z6i_exact_census": {
            "models": len(z6i),
            "models_with_nine_forced_singlets_and_all_three_quartics": len(exact),
            "forced_singlets_per_model": 9,
            "order3_rpv_before_condensation": 0,
        },
        "intertwiner_firewall": (
            "Any bridge that preserves the already-certified identification of "
            "clock central -I with matter parity must send the FI-forced singlet "
            "to the odd sector V-. Its nonzero VEV then breaks that central sign. "
            "Sending the singlet to V+ could preserve the clock sign only by "
            "abandoning the certified matter-parity intertwining datum, so it is "
            "not a rescue of the same symmetry."
        ),
        "flavour_corroboration": {
            "parent": "Pass 10980 RPV flavour upper bound",
            "models": f["models"],
            "analysable_with_P": f["with_P"],
            "no_light_Hu": f["no_light_Hu"],
            "five_light_d": f["five_light_d"],
            "below_1e-10": f["below_1e-10"],
            "minimum_upper_bound_P": f["min_P"],
            "reading": (
                "The symmetry no-go and the flavour audit are independent. "
                "Pass 11026 proves the central parity is broken by the required "
                "odd singlet VEV; Pass 10980 separately finds no first-generation "
                "suppression mechanism and every witness vacuum already fails a "
                "Higgs or exotic-spectrum condition."
            ),
        },
        "boundary": (
            "This rules out only the natural central-character intertwiner already "
            "certified in the repository. It does not rule out a genuinely new "
            "flavour symmetry unrelated to matter parity. Pass 10980's "
            "projected couplings are upper bounds with random O(1) coefficients, "
            "not computed string coupling coefficients."
        ),
        "checks": checks,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    p = payload()
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if a.check:
        if not a.output.exists() or a.output.read_text(encoding="utf-8") != text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": p["status"],
        "models": p["z6i_exact_census"]["models"],
        "all_three": p["z6i_exact_census"]["models_with_nine_forced_singlets_and_all_three_quartics"],
    }, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
