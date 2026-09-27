#!/usr/bin/env python3
"""Pass 11049: twelve rank-9 right-regular H27 projectors lift the Hesse affine plane."""
from __future__ import annotations
import argparse,itertools,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11049_right_regular_hesse_projectors.json"
DIRS=[(1,0),(0,1),(1,1),(1,2)]
G=list(itertools.product(range(3),repeat=3))
GIDX={g:i for i,g in enumerate(G)}

def mul(x,y):
    a,b,c=x; d,e,f=y
    return ((a+d)%3,(b+e)%3,(c+f-b*d)%3)

def right(g):
    M=np.zeros((27,27),dtype=np.int64)
    for x in G:M[GIDX[mul(x,g)],GIDX[x]]=1
    return M

def left(g):
    M=np.zeros((27,27),dtype=np.int64)
    for x in G:M[GIDX[mul(g,x)],GIDX[x]]=1
    return M
def rank_mod(A,p):
    M=A.copy()%p; n,m=M.shape; r=0
    for c in range(m):
        z=next((i for i in range(r,n) if M[i,c]%p),None)
        if z is None:continue
        M[[r,z]]=M[[z,r]]
        iv=pow(int(M[r,c]),-1,p);M[r]=(M[r]*iv)%p
        for i in range(n):
            if i!=r and M[i,c]%p:
                q=int(M[i,c]);M[i]=(M[i]-q*M[r])%p
        r+=1
        if r==n:break
    return r

def root_omega(p):
    roots=[x for x in range(2,p) if (x*x+x+1)%p==0]
    assert len(roots)==2
    return min(roots)

def affine_line(d,theta):
    a,b=d
    return sorted((r,s) for r,s in itertools.product(range(3),repeat=2)
                  if (r*a+s*b-theta)%3==0)

def char_vector(r,s,p,w):
    return np.array([pow(w,-(r*a+s*b),p) for a,b,c in G],dtype=np.int64)
def prime_packet(p):
    w=root_omega(p); inv3=pow(3,-1,p); I=np.eye(27,dtype=np.int64)
    P={}
    for di,d in enumerate(DIRS):
        R=right((d[0],d[1],0))%p
        for theta in range(3):
            P[di,theta]=(inv3*(I
                +pow(w,-theta,p)*R
                +pow(w,-2*theta,p)*((R@R)%p)))%p
            assert rank_mod(P[di,theta],p)==9
        assert np.array_equal(sum(P[di,t] for t in range(3))%p,I%p)
        for a in range(3):
            for b in range(3):
                if a!=b:assert not np.any((P[di,a]@P[di,b])%p)

    LX,LZ,LC=(left(g)%p for g in ((0,1,0),(1,0,0),(0,0,1)))
    for Q in P.values():
        for L in (LX,LZ,LC):
            assert np.array_equal((L@Q)%p,(Q@L)%p)

    cross=[]
    for di in range(4):
      for dj in range(di+1,4):
       for ti in range(3):
        for tj in range(3):
            A=P[di,ti];B=P[dj,tj]
            inter=27-rank_mod(np.vstack([(A-I)%p,(B-I)%p]),p)
            assert inter==1
            C=(A@B@A)%p
            rankC=rank_mod(C,p)
            m1=27-rank_mod((C-I)%p,p)
            m13=27-rank_mod((C-inv3*I)%p,p)
            zero_in_A=9-rankC
            assert (rankC,m1,m13,zero_in_A)==(7,1,6,2)

            Li=set(affine_line(DIRS[di],ti))
            Lj=set(affine_line(DIRS[dj],tj))
            meet=sorted(Li&Lj)
            assert len(meet)==1
            r,s=meet[0];v=char_vector(r,s,p,w)
            assert np.array_equal((A@v)%p,v%p)
            assert np.array_equal((B@v)%p,v%p)
            cross.append({
              "directions":[di,dj],"thetas":[ti,tj],
              "affine_line_meet":[r,s],
              "intersection_dimension":inter,
              "PQP_spectrum_on_imP":{"1":m1,"1/3":m13,"0":zero_in_A},
            })
    assert len(cross)==54
    return {"prime":p,"omega":w,"rank_each":9,"cross_pairs":cross}
def payload():
    cert=[prime_packet(p) for p in (103,109)]
    checks={
      "twelve_rank9_projectors":all(c["rank_each"]==9 for c in cert),
      "four_three_way_orthogonal_decompositions":True,
      "right_projectors_commute_with_left_H27":True,
      "54_cross_pairs_each_intersection_dim1":all(
          r["intersection_dimension"]==1 for c in cert for r in c["cross_pairs"]),
      "cross_PQP_spectrum_1_1over3x6_0x2":all(
          r["PQP_spectrum_on_imP"]=={"1":1,"1/3":6,"0":2}
          for c in cert for r in c["cross_pairs"]),
      "intersection_ray_is_unique_dual_affine_point":True,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11049.right-regular-hesse-projectors.v1",
      "status":"PASS_12_RANK9_PROJECTORS_REALIZE_THE_HESSE_AFFINE_PLANE_INSIDE_REG_H27",
      "headline":(
        "For each of the four noncentral H27/Z directions and three C3 eigenvalues, "
        "the right-regular spectral projector P=(1/3)sum_k omega^(-theta k)R_g^k "
        "has rank 9 and commutes with the left regular H27 action. Each direction "
        "gives an orthogonal 9+9+9 decomposition of C[H27]. Distinct directions "
        "meet in exactly one dimension, reproducing the affine-line intersection law."
      ),
      "projector_formula":"P_(d,theta)=(I+omega^-theta R_g+omega^-2theta R_g^2)/3",
      "geometry":{
        "projectors":12,
        "parallel_classes":4,
        "projectors_per_class":3,
        "rank":9,
        "same_class_intersection":0,
        "different_class_intersection":1,
        "common_ray":"the unique one-dimensional H27 character at the intersection of the two dual affine lines",
      },
      "cross_angle_spectrum":{
        "operator":"P Q P restricted to Im(P)",
        "squared_cosines":{"1":1,"1/3":6,"0":2},
        "reading":(
          "two nonparallel 9D compiler sectors share one exact character ray, "
          "have six principal overlaps |<.|.>|^2=1/3, and two orthogonal directions. "
          "They are not isoclinic; the Hesse lift has a rigid three-angle spectrum."
        ),
      },
      "split_prime_certificates":cert,
      "parents":["data/w33_pass11048_induced_h27_hesse_latent_modules.json"],
      "boundary":(
        "The projector theorem is an exact statement in the regular H27 carrier. "
        "The principal-angle language uses the complex unitary realization; the "
        "two split-prime calculations certify the corresponding cyclotomic rank "
        "and minimal-polynomial multiplicities without selecting a physical gauge."
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
    print(json.dumps({"status":p["status"],"cross_pairs":54},sort_keys=True))

if __name__=="__main__":raise SystemExit(main())
