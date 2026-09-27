#!/usr/bin/env python3
"""Pass 11051: tensor induction reduces the 27D compiler to nine qutrit Fourier blocks."""
from __future__ import annotations
import argparse,itertools,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11050_monomial_two_qutrit_latent_clifford as LAT

OUT=ROOT/"data/w33_pass11051_tensor_induction_fourier_compiler.json"
G=list(itertools.product(range(3),repeat=3));GIDX={g:i for i,g in enumerate(G)}

def mul(x,y):
    a,b,c=x;d,e,f=y
    return ((a+d)%3,(b+e)%3,(c+f-b*d)%3)

def left(g):
    M=np.zeros((27,27),dtype=np.int64)
    for x in G:M[GIDX[mul(g,x)],GIDX[x]]=1
    return M

def kron(A,B,p):return np.kron(A,B)%p
def star(A,p,w):
    w2=pow(w,2,p);out=np.zeros_like(A)
    mp={0:0,1:1,w:w2,w2:w}
    vals=set(int(x) for x in np.unique(A%p))
    assert vals<=set(mp)
    for x,y in mp.items():out[A%p==x]=y
    return out

def support_components(B,p):
    rows={i:set(np.where(B[i]%p!=0)[0].tolist()) for i in range(27)}
    cols={j:set(np.where(B[:,j]%p!=0)[0].tolist()) for j in range(27)}
    seen_r=set();seen_c=set();out=[]
    for r0 in range(27):
        if r0 in seen_r:continue
        rs={r0};cs=set();todo=[("r",r0)]
        seen_r.add(r0)
        while todo:
            typ,k=todo.pop()
            if typ=="r":
                for j in rows[k]:
                    if j not in seen_c:seen_c.add(j);cs.add(j);todo.append(("c",j))
            else:
                for i in cols[k]:
                    if i not in seen_r:seen_r.add(i);rs.add(i);todo.append(("r",i))
        out.append((sorted(rs),sorted(cs)))
    return out
def prime_theta_packet(p,theta):
    w=LAT.root_omega(p);I3=np.eye(3,dtype=np.int64)
    AX=LAT.perm9(lambda b,c:((b+1)%3,c))%p
    AC=LAT.perm9(lambda b,c:(b,(c+1)%3))%p
    S=LAT.perm9(lambda b,c:(b,(c+b)%3))%p
    AZ=(pow(w,theta,p)*S)%p
    X3=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=np.int64)
    Z3=np.diag([1,w,pow(w,2,p)]).astype(np.int64)
    UX=kron(AX,X3,p)
    UZ=kron(AZ,Z3,p)
    UC=kron(AC,(w*I3)%p,p)

    v0=np.zeros(27,dtype=np.int64);v0[:3]=1
    cols=[]
    for a,b,c in G:
        U=LAT.mpow(UZ,a,p)@LAT.mpow(UX,b,p)%p
        U=U@LAT.mpow(UC,c,p)%p
        cols.append((U@v0)%p)
    B=np.column_stack(cols)%p

    Bs=star(B,p,w).T%p
    assert np.array_equal((Bs@B)%p,(3*np.eye(27,dtype=np.int64))%p)
    for U,g in ((UX,(0,1,0)),(UZ,(1,0,0)),(UC,(0,0,1))):
        assert np.array_equal((U@B)%p,(B@left(g))%p)

    comps=support_components(B,p)
    assert len(comps)==9 and all(len(r)==len(c)==3 for r,c in comps)
    assert np.count_nonzero(B)==81
    Fexp=[[0,0,0],[0,1,2],[0,2,1]]
    block_rows=[]
    root_exp={1:0,w:1,pow(w,2,p):2}
    for rs,cs in comps:
        block=B[np.ix_(rs,cs)]%p
        assert np.count_nonzero(block)==9
        E=np.array([[root_exp[int(x)] for x in row] for row in block],dtype=int)
        D=(E-E[:,[0]]-E[[0],:]+E[0,0])%3
        assert D.tolist()==Fexp
        block_rows.append({"rows":rs,"columns":cs,"dephased_exponents":D.tolist()})

    return {
      "prime":p,"theta":theta,"omega":w,
      "rank":27,"nonzero_entries":81,"nonzeros_per_column":3,
      "Gram":"B^* B = 3 I27",
      "support_components":block_rows,
      "components":9,"component_shape":[3,3],
      "every_component_dephases_to_F3":True,
      "intertwines_left_regular_H27":True,
    }
def payload():
    cert=[prime_theta_packet(p,t) for p in (103,109) for t in range(3)]
    checks={
      "all_six_replays_rank27":all(c["rank"]==27 for c in cert),
      "all_have_81_nonzeros":all(c["nonzero_entries"]==81 for c in cert),
      "nine_disjoint_3x3_blocks":all(c["components"]==9 and c["component_shape"]==[3,3] for c in cert),
      "every_block_is_F3_up_to_monomial_phases":all(c["every_component_dephases_to_F3"] for c in cert),
      "exact_Gram_3I":all(c["Gram"]=="B^* B = 3 I27" for c in cert),
      "exact_intertwiner":all(c["intertwines_left_regular_H27"] for c in cert),
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11051.tensor-induction-fourier-compiler.v1",
      "status":"PASS_TENSOR_INDUCTION_COMPILER_IS_NINE_PARALLEL_QUTRIT_FOURIER_BLOCKS",
      "headline":(
        "The identity V tensor Ind_L(theta) ~= Ind_L(Res_L V tensor theta) "
        "~= Ind_L(Reg L) ~= Reg(H27) is now an explicit 27x27 matrix. A cyclic "
        "vector on one coset generates an orthogonal regular orbit. Its synthesis "
        "matrix has exactly 81 nonzero mu3 entries and splits into nine disjoint "
        "3x3 blocks, each dephasing exactly to the qutrit Fourier matrix F3."
      ),
      "tensor_identity":{
        "step1":"V tensor Ind_L(theta) ~= Ind_L(Res_L(V) tensor theta)",
        "step2":"Res_L(V)=Reg(C3)",
        "step3":"theta tensor Reg(C3) ~= Reg(C3)",
        "step4":"Ind_L(Reg L) ~= Reg(H27)",
      },
      "compiler":{
        "dimension":27,
        "balanced_blocks":9,
        "block_size":3,
        "active_F3_depth":1,
        "matrix_nonzeros":81,
        "generic_27x27_unitary_required":False,
        "normalization":"divide each balanced block by sqrt(3)",
      },
      "physical_reading":(
        "In a mode encoding this is nine parallel tritters plus permutations and "
        "mu3 phase gauges. In a tensor-qutrit encoding it is one qutrit Fourier "
        "operation applied uniformly over nine coset labels. The latent dressing "
        "therefore does not require a generic U(9) or U(27) synthesis."
      ),
      "split_prime_theta_certificates":cert,
      "parents":[
        "data/w33_pass11048_induced_h27_hesse_latent_modules.json",
        "data/w33_pass11050_monomial_two_qutrit_latent_clifford.json",
      ],
      "boundary":(
        "This closes the finite representation compiler for one induced H27 slice. "
        "It does not implement the E6 cubic 54D tangent pulse or provide measured "
        "loss/fidelity for a photonic tritter array."
      ),
      "checks":checks,
    }
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args();p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"F3_blocks":9,"nonzeros":81},sort_keys=True))

if __name__=="__main__":raise SystemExit(main())
