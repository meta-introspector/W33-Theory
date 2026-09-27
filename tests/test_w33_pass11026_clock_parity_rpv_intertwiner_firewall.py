import importlib.util, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p11026",ROOT/"analysis"/"w33_pass11026_clock_parity_rpv_intertwiner_firewall.py")
P=importlib.util.module_from_spec(S); S.loader.exec_module(P)
C=json.loads((ROOT/"data"/"w33_pass11026_clock_parity_rpv_intertwiner_firewall.json").read_text())

def test_replay():
    assert P.payload()==C

def test_central_parity_selection_rule():
    s=C["selection_rule"]
    assert set(s["ordinary_yukawa_signs"].values())=={1}
    assert set(s["rpv_cubic_signs"].values())=={-1}
    assert set(s["forced_singlet_times_rpv_signs"].values())=={1}

def test_all_23_z6i_vacua_regenerate():
    z=C["z6i_exact_census"]
    assert z["models"]==23
    assert z["models_with_nine_forced_singlets_and_all_three_quartics"]==23
    assert z["forced_singlets_per_model"]==9

def test_signed_clock_is_12_plus_12():
    assert C["checks"]["signed_carrier_is_12_plus_12"]

def test_committed_pass10980_corroboration():
    f=C["flavour_corroboration"]
    assert f["models"]==23
    assert f["analysable_with_P"]==16
    assert f["no_light_Hu"]==7
    assert f["five_light_d"]==16
    assert f["below_1e-10"]==0
    assert f["minimum_upper_bound_P"] > 2e-3
