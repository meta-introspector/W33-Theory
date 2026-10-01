#!/usr/bin/env python3
"""Reconstructed Cartan grade-pair metric and the zero-fit equal-mirror-spacing vacuum."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11260_exact_semisimple_cartan as P60
import w33_pass11261_g26_cartan_coordinate_map as P61
import w33_pass11258_g26_cubic_cp_interface as P58
OUT=ROOT/"data/w33_20261001_cartan_kinetic_wall_vacuum.json"

def build_pairing():
    toe=P60._load(ROOT/"tools/toe_e8_z3graded_bracket_jacobi.py","toe_cartan_metric")
    basis=np.load(ROOT/"artifacts/e6_27rep_basis_export/E6_basis_78.npy").astype(complex)
    proj=toe.E6Projector(basis)
    br=toe.E8Z3Bracket(e6_projector=proj,cubic_triads=toe._load_signed_cubic_triads(),
        scale_g1g1=1.0,scale_g2g2=-1.0/6.0,scale_e6=1.0,scale_sl3=1.0/6.0)
    sl3=[]
    for i in range(3):
        for j in range(3):
            if i!=j:
                q=np.zeros((3,3),complex);q[i,j]=1;sl3.append(q)
    sl3 += [np.diag([1,-1,0]).astype(complex),np.diag([0,1,-1]).astype(complex)]
    sl3=np.array(sl3);sf=sl3.reshape(8,-1);sleft=np.linalg.inv(sf@sf.conj().T)@sf.conj()
    zero=toe.E8Z3.zero()
    eb=[toe.E8Z3(basis[i],zero.sl3,zero.g1,zero.g2) for i in range(78)]
    eb += [toe.E8Z3(zero.e6,sl3[i],zero.g1,zero.g2) for i in range(8)]
    for gr in (1,2):
        for i in range(81):
            m=np.zeros((27,3),complex);m.flat[i]=1
            eb.append(toe.E8Z3(zero.e6,zero.sl3,m if gr==1 else zero.g1,
                               m if gr==2 else zero.g2))
    def coords(e):
        rhs=np.einsum("aij,ji->a",basis,e.e6);ce6=proj.gram_inv@rhs
        return np.r_[ce6,sleft@e.sl3.reshape(-1),e.g1.ravel(),e.g2.ravel()]
    def ad(v,gr):
        X=np.array(v,complex).reshape(27,3)
        e=toe.E8Z3(zero.e6,zero.sl3,X if gr==1 else zero.g1,X if gr==2 else zero.g2)
        raw=18*np.column_stack([coords(br.bracket(e,b)) for b in eb])
        res=float(np.max(np.abs(raw-np.rint(raw.real))))
        return np.rint(raw.real).astype(np.int64),res
    kernel,_,_=P61.slice_matrices();A1=[];A2=[];res=[]
    for q in kernel:
        u,r1=ad(q,1);v,r2=ad(q,2);A1.append(u);A2.append(v);res += [r1,r2]
    G=sp.Matrix([[sp.Rational(int(np.trace(A1[i]@A2[j])),18**2)
                  for j in range(3)] for i in range(3)])
    return G,max(res)
def local_wall_theorem(G):
    x,y,z=P61.x,P61.y,P61.z
    F=sp.expand(P61.C3*P61.N9*P61.S9)
    pt={x:0,y:0,z:1}; assert F.subs(pt)==1
    grad=[sp.diff(F,v).subs(pt) for v in (x,y,z)]
    HF=sp.Matrix([[sp.diff(F,a,b).subs(pt) for b in (x,y)] for a in (x,y)])
    assert grad==[0,0,21] and HF==sp.Matrix([[0,12],[12,0]])
    G2=G[:2,:2]; q0=G[2,2]
    Hr=2*HF-sp.Rational(42,1)*G2/q0
    Hi=-2*HF-sp.Rational(42,1)*G2/q0
    assert Hr.trace()<0 and Hr.det()>0 and Hi.trace()<0 and Hi.det()>0
    t=np.exp(1j*np.pi/9);std=np.array([1,t**2,t**16],complex)
    _,W,sing=P58.induced(std,True,True);J=P58.jarlskog(W)
    return {
      "wall_product":"F=C3*N9*S9 = product of all 21 mirror forms up to a nonzero scalar",
      "projective_functional":"Phi=|F|^2/(q^dag G q)^21",
      "candidate_cartan_ray":[0,0,1],
      "candidate_standard_qutrit":"[1,zeta9,zeta9^-1] up to global phase",
      "F_gradient_at_ray":[int(v) for v in grad],
      "F_tangent_hessian":[[int(v) for v in row] for row in HF.tolist()],
      "log_Phi_real_tangent_hessian":[[str(v) for v in row] for row in Hr.tolist()],
      "log_Phi_imag_tangent_hessian":[[str(v) for v in row] for row in Hi.tolist()],
      "real_hessian_trace":str(Hr.trace()),"real_hessian_det":str(Hr.det()),
      "imag_hessian_trace":str(Hi.trace()),"imag_hessian_det":str(Hi.det()),
      "strict_local_projective_maximum":True,
      "cp_interface_jarlskog":float(J),"cp_zero_tolerance":abs(J)<1e-14,
      "polar_abs2":np.round(np.abs(W)**2,15).tolist(),
      "singular_values":[float(v) for v in sing],
    }
def payload():
    G,res=build_pairing()
    assert G==sp.Matrix([[590,20,0],[20,980,0],[0,0,300]])
    minors=[G[:k,:k].det() for k in (1,2,3)];assert all(v>0 for v in minors)
    vac=local_wall_theorem(G)
    return {
      "schema":"w33.20261001.cartan-kinetic-wall-vacuum.v1",
      "status":"PASS_RECONSTRUCTED_E8_CARTAN_GRADE_PAIR_METRIC_AND_CP_SYMMETRIC_EQUAL_MIRROR_LOCAL_VACUUM",
      "cartan_grade_pair_metric":{
        "definition":"G_ij=Tr(ad(q_i in g1) ad(q_j in g2)) with the frozen E8 bracket",
        "matrix":[[int(v) for v in row] for row in G.tolist()],
        "principal_minors":[int(v) for v in minors],
        "eigenvalues_exact":["300","785-5*sqrt(1537)","785+5*sqrt(1537)"],
        "positive_definite":True,"max_integral_rounding_residual":res,
        "verification_scope":"adjoints reconstructed by rounding numerical bracket coordinates; exact E8 trace identity not independently certified",
        "physical_scope":"positive after the natural g1<->g2 grade-exchange identification"
      },
      "zero_fit_wall_vacuum":vac,
      "result":(
        "Integer reconstruction from numerical bracket coordinates supplies the displayed positive Cartan grade-pair Gram form. The coefficient-free "
        "equal-mirror-spacing functional has an exact strict local maximum at [0:0:1], but its qutrit "
        "image [1,zeta9,zeta9^-1] gives a CP-conserving postselected interface. Thus mirror repulsion "
        "alone does not select the observed orientation; an additional odd/oriented datum is required."
      ),
      "boundary":(
        "The local maximum is exact but global maximality is not proved. The positive form uses the "
        "grade-exchange real identification and is a kinetic candidate conditional on this identification, not "
        "a measured normalization. No Standard-Model mass or mixing angle is fitted."
      ),
      "parents":["data/w33_pass11260_exact_semisimple_cartan.json",
                 "data/w33_pass11261_g26_cartan_coordinate_map.json",
                 "data/w33_20261001_g26_e8_mirror_separation.json",
                 "data/w33_pass11258_g26_cubic_cp_interface.json"],
      "checks":{"metric_integer_reconstruction":True,"metric_positive":True,"ray_stationary_projectively":True,
                "real_tangent_negative":True,"imag_tangent_negative":True,"cp_interface_zero":vac["cp_zero_tolerance"]}
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true");a=ap.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not OUT.exists() or OUT.read_text()!=txt:raise SystemExit("certificate drift")
    else:OUT.write_text(txt)
    print(json.dumps({"status":p["status"],"G":p["cartan_grade_pair_metric"]["matrix"],
                      "J":p["zero_fit_wall_vacuum"]["cp_interface_jarlskog"]}))

if __name__=="__main__":raise SystemExit(main())
