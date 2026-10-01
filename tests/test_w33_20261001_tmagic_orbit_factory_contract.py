import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location('orbit_contract',ROOT/'analysis/w33_20261001_tmagic_orbit_factory_contract.py')
M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

def test_full_orbit_and_all_input_operator_branch_replay():
    got=M.payload()
    stored=json.loads((ROOT/'data/w33_20261001_tmagic_orbit_factory_contract.json').read_text())
    for k in ('orbit','exact_theorem','frame_loss','prior_owners','checks'):
        assert got[k]==stored[k]
    assert got['retained_frame']['bell_branches_checked']==648
    assert got['retained_frame']['operator_error']<1e-12
    assert got['numerical_controls']['dephasing_channel_operator_error']<1e-12
    assert got['frame_loss']['entanglement_breaking']

def test_pauli_twirl_is_an_arbitrary_operator_identity():
    assert 'Tr(rho) I/3' in M.exact_pauli_twirl()
