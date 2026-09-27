import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11059_diagonal_background_seven_fourier_coordinates.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11059_diagonal_background_seven_fourier_coordinates.json").read_text())
def test_replay(): assert P.payload()==C
def test_sparse(): assert C["exact_formula"]["nonzero_coordinates"]==7
def test_orientations(): assert C["checks"]["plus_support7"] and C["checks"]["minus_support7"]
