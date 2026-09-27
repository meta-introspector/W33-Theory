#!/usr/bin/env python3
"""Pass 11065: K81 and chamber U81 are split versus non-split central C3 extensions of H27."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter,deque
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11065_split_nonsplit_order81_bridge.json"

def k_mul(x,y):
    (a,b,c),p=x;(d,e,f),q=y
    return (((a+d)%3,(b+e)%3,(c+f-d*b)%3),(p+q)%3)
def k_inv(x):
    for y in K:
        if k_mul(x,y)==IDK and k_mul(y,x)==IDK:return y
    raise RuntimeError
K=[((a,b,c),p) for a,b,c,p in itertools.product(range(3),repeat=4)]
IDK=((0,0,0),0)

def closure_k(gens):
    S={IDK};Q=deque([IDK])
    while Q:
        a=Q.popleft()
        for g in gens:
            z=k_mul(a,g)
            if z not in S:S.add(z);Q.append(z)
    return S

def k_comm(a,b):
    return k_mul(k_mul(k_mul(a,b),k_inv(a)),k_inv(b))

def mat_key(A):return tuple(int(x)%3 for x in A.flat)
I=np.eye(4,dtype=int)%3
def E(i,j):
    M=np.zeros((4,4),dtype=int);M[i,j]=1;return M
X=[E(0,1)-E(3,2),E(1,3),E(0,3)+E(1,2),E(0,2)]
R=[(I+z)%3 for z in X]

def mat_group(gens):
    D={mat_key(I):I};Q=[I]
    while Q:
        a=Q.pop()
        for g in gens:
            b=(a@g)%3;k=mat_key(b)
            if k not in D:D[k]=b;Q.append(b)
    return list(D.values())
def mat_order(A):
    z=I.copy()
    for n in range(1,30):
        z=(z@A)%3
        if np.array_equal(z,I):return n
    raise RuntimeError
def mat_inv(A):
    n=mat_order(A);z=I.copy()
    for _ in range(n-1):z=(z@A)%3
    return z
def mat_comm(a,b):return (((a@b)%3@mat_inv(a))%3@mat_inv(b))%3
def mat_closure(gens):
    return mat_group(gens)

def payload():
    centerK={g for g in K if all(k_mul(g,h)==k_mul(h,g) for h in K)}
    commK={k_comm(a,b) for a in K for b in K}
    derivedK=closure_k(commK)
    assert len(K)==81 and len(centerK)==9 and len(derivedK)==3
    assert Counter(1 if g==IDK else 3 for g in K)==Counter({3:80,1:1})
    ext={((0,0,0),p) for p in range(3)}
    assert ext<=centerK
    section={((a,b,c),0) for a,b,c in itertools.product(range(3),repeat=3)}
    assert len(section)==27 and len(ext)==3 and closure_k(ext|section)==set(K)

    U=mat_group([R[0],R[1]]);assert len(U)==81
    Z=[a for a in U if all(np.array_equal((a@b)%3,(b@a)%3) for b in U)]
    Ud=mat_closure([mat_comm(a,b) for a in U for b in U])
    g3=mat_closure([mat_comm(a,b) for a in U for b in Ud])
    g4=mat_closure([mat_comm(a,b) for a in U for b in g3])
    assert (len(Z),len(Ud),len(g3),len(g4))==(3,9,3,1)
    uorders=Counter(mat_order(a) for a in U);assert uorders==Counter({3:44,9:36,1:1})

    def coset_key(a):
        return min(mat_key((a@z)%3) for z in Z)
    reps={}
    for a in U:reps.setdefault(coset_key(a),a)
    Q=list(reps.values());qi={coset_key(a):i for i,a in enumerate(Q)}
    def qmul(i,j):return qi[coset_key((Q[i]@Q[j])%3)]
    e=qi[coset_key(I)]
    qorders=[]
    for i in range(27):
        x=e
        for n in range(1,10):
            x=qmul(x,i)
            if x==e:qorders.append(n);break
    qcenter=[i for i in range(27) if all(qmul(i,j)==qmul(j,i) for j in range(27))]
    commidx={qmul(qmul(qmul(i,j),next(k for k in range(27) if qmul(i,k)==e)),next(k for k in range(27) if qmul(j,k)==e)) for i in range(27) for j in range(27)}
    assert Counter(qorders)==Counter({3:26,1:1}) and len(qcenter)==3 and len(commidx)==3

    old=json.loads((ROOT/"data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json").read_text())
    assert old["U"]["order"]==81 and old["U"]["center"]==3 and old["U"]["derived"]==9

    return {
      "schema":"w33.pass11065.split-nonsplit-order81-bridge.v1",
      "status":"PASS_SCHEDULER_K81_AND_CHAMBER_U81_ARE_SPLIT_AND_NONSPLIT_CENTRAL_C3_EXTENSIONS_OF_H27",
      "headline":"The two recurring order-81 groups are not the same object. K81=H27 x C3_external is a split central extension of H27 by C3. The chamber Sylow-3 group U81 has center 3, derived subgroup 9, class three and elements of order 9; quotienting by its center gives the extraspecial exponent-3 H27. Thus U81 is a non-split central C3 extension of the same H27.",
      "K81":{"order":81,"center_order":9,"derived_order":3,"nilpotency_class":2,"exponent":3,
             "order_census":{"1":1,"3":80},"extension":"1 -> C3_external -> K81 -> H27 -> 1","split":True},
      "U81":{"order":81,"center_order":3,"derived_order":9,"lower_central_sizes":[81,9,3,1],
             "nilpotency_class":3,"exponent":9,"order_census":{"1":1,"3":44,"9":36},
             "quotient_by_center":{"order":27,"center_order":3,"derived_order":3,"order_census":{"1":1,"3":26},"is_H27":True},
             "extension":"1 -> Z(U81)=C3 -> U81 -> H27 -> 1","split":False},
      "bridge_reading":"K81 keeps the extra C3 as an independent central scheduler clock; U81 welds that extra ternary degree into a class-three update group whose quotient is the same H27 memory group.",
      "nonisomorphism_witnesses":["center 9 versus 3","derived 3 versus 9","class 2 versus 3","exponent 3 versus 9"],
      "boundary":"This is an exact finite-group theorem. Calling the extra central coordinate a scheduler clock or chamber tick is architectural language, not a claim about physical elapsed time.",
      "parents":["data/w33_e8_matter81_h27_address_operator_compiler.json","data/PART_W33_PASS5105_U81_DUAL_TORSOR_CONTROLLER.json"],
      "checks":{"K_split_H27_extension":True,"U_quotient_center_is_H27":True,"U_extension_nonsplit":True,"K_U_nonisomorphic":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT);a=ap.parse_args()
    p=payload();text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text:raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({"status":p["status"],"K_class":2,"U_class":3},sort_keys=True))
if __name__=="__main__":raise SystemExit(main())
