import importlib.util,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("p",R/"analysis"/"w33_pass11058_triangle_factor_kernel_interference.py")
P=importlib.util.module_from_spec(S);S.loader.exec_module(P)
C=json.loads((R/"data"/"w33_pass11058_triangle_factor_kernel_interference.json").read_text())
def test_replay(): assert P.payload()==C
def test_kernel(): assert C["global_kernel"]["common_kernel_all_ten_dimension"]==1 and C["global_kernel"]["kernel_dimension_of_full_D_plus"]==15
def test_interference(): assert C["global_kernel"]["extra_dark_dimensions_created_by_interference"]==14
