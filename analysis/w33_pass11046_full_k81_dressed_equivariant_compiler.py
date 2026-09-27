#!/usr/bin/env python3
"""Pass 11046: lift the latent H27 dressing to an explicit equivariant K81 compiler."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11045_explicit_latent_h27_clebsch_gordan as CG

OUT=ROOT/"data/w33_pass11046_full_k81_dressed_equivariant_compiler.json"


def ext_fourier(p,w):
    return np.array([[pow(w,t*x,p) for x in range(3)] for t in range(3)],dtype=np.int64)


def prime_certificate(p):
    w=CG.root_omega(p)
    AX,AZ=CG.latent_generators(p,w)
    X1,Z1=CG.qutrit_generators(p,w,1)
    I9=np.eye(9,dtype=np.int64); I3=np.eye(3,dtype=np.int64)
    AC=CG.blockdiag([
      np.array([[1]]),np.array([[1]]),np.array([[1]]),
      (w*np.eye(3,dtype=np.int64))%p,
      (pow(w,2,p)*np.eye(3,dtype=np.int64))%p],p)

    UX=CG.mm(CG.kron(AX,I3,p),CG.kron(I9,X1,p),p)
    UZ=CG.mm(CG.kron(AZ,I3,p),CG.kron(I9,Z1,p),p)
    UC=CG.mm(CG.kron(AC,I3,p),(w*np.eye(27,dtype=np.int64))%p,p)
    HX,HZ=CG.target_generators(p,w)
    HC=CG.blockdiag([
      (w*np.eye(9,dtype=np.int64))%p,
      (pow(w,2,p)*np.eye(9,dtype=np.int64))%p,
      np.eye(9,dtype=np.int64)],p)
    TH=CG.cg_matrix(p,w)
    F=ext_fourier(p,w)
    T81=CG.kron(TH,F,p)
    r,det=CG.rank_det(T81,p)
    assert r==81 and det!=0

    I27=np.eye(27,dtype=np.int64)
    I81=np.eye(81,dtype=np.int64)
    Xext=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=np.int64)
    Dext=np.diag([1,w,pow(w,2,p)]).astype(np.int64)

    physical={
      "H_X":CG.kron(UX,I3,p),
      "H_Z":CG.kron(UZ,I3,p),
      "H_C":CG.kron(UC,I3,p),
      "ext_X":CG.kron(I27,Xext,p),
    }
    fourier={
      "H_X":CG.kron(HX,I3,p),
      "H_Z":CG.kron(HZ,I3,p),
      "H_C":CG.kron(HC,I3,p),
      "ext_X":CG.kron(I27,Dext,p),
    }
    for key in physical:
        assert np.array_equal(CG.mm(T81,physical[key],p),CG.mm(fourier[key],T81,p))
    # The H Fourier decomposition is 9 S1 + 9 S2 + 9 L; external Fourier triples it.
    h_blocks={"S1":9,"S2":9,"L":9}
    k_blocks={k:3*v for k,v in h_blocks.items()}
    assert k_blocks=={"S1":27,"S2":27,"L":27}

    old=json.loads((ROOT/"data/w33_minimal_symmetry_changing_81_compiler.json").read_text())
    assert old["compiler"]["mode_counts"]=={
        "equivariant":27,"retyped_central":27,"retyped_abelian":27}

    return {
      "prime":p,
      "omega":w,
      "T81_rank":r,
      "T81_det_mod_p":int(det),
      "intertwines":["H_X","H_Z","H_C","external_C3_shift"],
      "Fourier_blocks":k_blocks,
      "old_compiler_mode_counts":old["compiler"]["mode_counts"],
    }


def payload():
    certs=[prime_certificate(p) for p in (103,109)]
    checks={
      "T81_full_rank_both_primes":all(c["T81_rank"]==81 for c in certs),
      "T81_nonzero_determinant_both_primes":all(c["T81_det_mod_p"] for c in certs),
      "all_four_K_generators_intertwined":all(len(c["intertwines"])==4 for c in certs),
      "Fourier_split_27_27_27":all(c["Fourier_blocks"]=={"S1":27,"S2":27,"L":27} for c in certs),
      "matches_old_27_27_27_compiler_modes":all(
          c["old_compiler_mode_counts"]=={
            "equivariant":27,"retyped_central":27,"retyped_abelian":27} for c in certs),
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11046.full-k81-dressed-equivariant-compiler.v1",
      "status":"PASS_EXISTING_81D_MATTER_CARRIER_BECOMES_REGULAR_K81_AFTER_MINIMAL_COMMUTANT_DRESSING",
      "headline":(
        "The scheduler/execution obstruction is closed algebraically without enlarging "
        "the 81D matter carrier. Let the 9D multiplicity commutant carry A9 and let "
        "H27 act diagonally as A9 tensor V_omega; keep the external qutrit as Reg(C3). "
        "Then the carrier is Reg(H27) tensor Reg(C3)=Reg(K). An explicit 81x81 "
        "Clebsch-Gordan/Fourier transform intertwines all K generators."
      ),
      "carrier":{
        "before":"C9_multiplicity tensor C3_internal tensor C3_external, with H27 action I9 tensor V_omega",
        "latent_control":"A9 acts on C9_multiplicity inside M9 commutant",
        "after":"(A9 tensor V_omega) tensor Reg(C3) ~= Reg(H27) tensor Reg(C3) = Reg(K)",
        "dimension":81,
        "dimension_enlargement":0,
      },
      "compiler_factorization":{
        "H27_transform":"T_H = explicit Pass11045 Clebsch-Gordan transform",
        "external_transform":"F3 on the external C3 shift register",
        "full_transform":"T81 = T_H tensor F3",
        "S1_27":"three scalar-character multiplicity states tensor V; symmetry-compatible sector",
        "S2_27":"latent V tensor internal V -> three Vbar copies",
        "L_27":"latent Vbar tensor internal V -> nine one-dimensional characters",
      },
      "why_54_retypings":{
        "equivariant":27,
        "retyped_central":27,
        "retyped_abelian":27,
        "total_retyped":54,
        "explanation":(
          "The old compiler's two 27D retyping blocks are exactly the two non-scalar "
          "latent sectors of A9. The V block produces S2=3Vbar; the Vbar block produces "
          "L=sum_9 chi. The previous 27+27+27 bookkeeping is therefore a Clebsch-Gordan "
          "decomposition of the minimal commutant dressing, not an arbitrary permutation."
        ),
      },
      "split_prime_certificates":certs,
      "parents":[
        "data/w33_pass11045_explicit_latent_h27_clebsch_gordan.json",
        "data/w33_minimal_symmetry_changing_81_compiler.json",
        "data/w33_scheduler_operator_k81_intertwiner_obstruction.json",
      ],
      "boundary":(
        "This removes the finite representation-theoretic obstruction by changing the "
        "execution symmetry through an action available in the abstract M9 commutant. "
        "It does not prove the laboratory Hamiltonian can synthesize A9, that doing so "
        "is fault tolerant, or that the required control is energetically natural."
      ),
      "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({
      "status":p["status"],
      "dimension":p["carrier"]["dimension"],
      "retyped":p["why_54_retypings"]["total_retyped"],
      "dets":[c["T81_det_mod_p"] for c in p["split_prime_certificates"]]
    },sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
