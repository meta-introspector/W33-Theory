#!/usr/bin/env python3
"""Balanced G26 mirror barrier, exact T-magic stationary orbit, and Cartan kinetic pullback."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass11269_g26_qutrit_dictionary as G69
import w33_pass11255_11259_g26_common as GC
OUT=ROOT/"data/w33_20261001_balanced_g26_tmagic_vacuum.json"

# Exact cyclotomic field Q(zeta_36), enough for zeta_9 and i.
r=sp.Symbol("r")
PHI=sp.Poly(sp.cyclotomic_poly(36,r),r,domain=sp.QQ)
def rem(e):
    return sp.rem(sp.Poly(sp.expand(e),r,domain=sp.QQ),PHI).as_expr()
def inv(e):
    return sp.invert(sp.Poly(rem(e),r,domain=sp.QQ),PHI).as_expr()
def red(e):
    n,d=sp.fraction(sp.cancel(e))
    return rem(rem(n)*inv(d))
def cj(e):
    return red(sp.expand(e).subs(r,r**35))
def re(e):
    return red((e+cj(e))/2)

zeta=r**4;omega=r**12;ii=r**9
m=[1,zeta,r**32]
u1=[1,omega*zeta,omega**2*r**32]
u2=[1,omega**2*zeta,omega*r**32]
TBAS=[u1,[ii*x for x in u1],u2,[ii*x for x in u2]]
SIC=[(1,-omega**k,0) for k in range(3)]+[(1,0,-omega**k) for k in range(3)]+[(0,1,-omega**k) for k in range(3)]
STAB=[(1,0,0),(0,1,0),(0,0,1)]+[(1,omega**a,omega**b) for a in range(3) for b in range(3)]
NORMALS=SIC+STAB
def inner(n,v):
    return red(sum(cj(n[j])*v[j] for j in range(3)))
def exact_local_certificate():
    ratios=[]
    for n in NORMALS:
        a=inner(n,m)
        ratios.append([red(inner(n,u)/a) for u in TBAS])
    H=sp.zeros(4)
    for j in range(4):
        for k in range(4):
            x=sp.Rational(42) if j==k else sp.Rational(0)
            for q in ratios:
                x += 4*re(q[j])*re(q[k])-2*re(cj(q[j])*q[k])
            H[j,k]=red(x)
    want=sp.Matrix([[42,0,-12,0],[0,42,0,12],[-12,0,42,0],[0,12,0,42]])
    assert H==want

    def piece(normals):
        rr=[]
        for n in normals:
            aa=inner(n,m)
            rr.append([red(inner(n,u)/aa) for u in TBAS])
        Q=sp.zeros(4)
        for j in range(4):
            for k in range(4):
                xx=sp.Rational(2*len(normals)) if j==k else sp.Rational(0)
                for q in rr:
                    xx += 4*re(q[j])*re(q[k])-2*re(cj(q[j])*q[k])
                Q[j,k]=red(xx)
        return Q
    HS=piece(SIC);HM=piece(STAB)
    assert HS==sp.Matrix([[18,0,-36,0],[0,18,0,36],[-36,0,18,0],[0,36,0,18]])
    assert HM==sp.Matrix([[24,0,24,0],[0,24,0,-24],[24,0,24,0],[0,-24,0,24]])
    assert HS+HM==H
    Sm=[0,0,0]
    for n in NORMALS:
        a=inner(n,m);den=red(cj(a)*a)
        for j in range(3):
            Sm[j]+=red(n[j]*a/den)
    assert [red(Sm[j]-7*m[j]) for j in range(3)]==[0,0,0]
    lam=sp.Symbol("lambda")
    cp=sp.expand(H.charpoly(lam).as_expr())
    assert sp.expand(cp-(lam-30)**2*(lam-54)**2)==0
    # G26 quotient coordinates of unit T state. Raw m has norm^2=3.
    x,y,z=GC.VARS
    fs=GC.invariants()
    raw=[red(f.subs({x:m[0],y:m[1],z:m[2]})) for f in fs]
    assert raw==[0,0,-27]
    unit=["0","0","-1/729"]
    return {"stationary_equation":"sum_n n<n|m>/|<n|m>|^2 = 7 m",
            "real_projective_hessian":H.tolist(),
            "hessian_eigenvalues":[30,30,54,54],
            "hessian_characteristic_polynomial":str(cp),
            "hessian_orbit_split":{
                "SIC_9":HS.tolist(),
                "SIC_eigenvalues":[-18,-18,54,54],
                "MUB_12":HM.tolist(),
                "MUB_eigenvalues":[0,0,48,48],
                "complementarity":"the MUB sector adds +48 exactly on the two SIC-unstable modes, while the SIC sector supplies +54 on the two MUB-flat modes"
            },
            "raw_basic_invariants":list(map(str,raw)),
            "unit_basic_invariants":unit}
def kinetic_certificate():
    a,w=sp.symbols("a w")
    core=sp.Matrix([[a,a**2,1],[a,w*a**2,w**2],[a,w**2*a**2,w]])
    # conjugation sends w -> w^2. Row phases t^8 and -t have unit norm and drop out.
    gram=core.applyfunc(lambda q:q.subs(w,w**2)).T*core
    def rw(q):
        return sp.rem(sp.Poly(sp.expand(q),w),sp.Poly(w**2+w+1,w)).as_expr()
    G=gram.applyfunc(rw)
    want=sp.diag(3*a**2,3*a**4,3)
    assert G==want
    return {"standard_hermitian_pullback":"3*diag(a^2,a^4,1), a^3=2",
            "rescaled_coordinates":"q=(a*x,a^2*y,z) gives 3*I",
            "unique_up_to_overall_scale":"G26 standard reflection representation is irreducible and unitary"}

def projective_group():
    w=np.exp(2j*np.pi/3)
    S=[q/np.linalg.norm(q) for q in G69.num(G69.SIC)]
    T=[q/np.linalg.norm(q) for q in G69.num(G69.STAB)]
    gens=[np.eye(3)-2*np.outer(q.conj(),q) for q in S]
    gens += [np.eye(3)+(w-1)*np.outer(q.conj(),q) for q in T]
    def key(M):
        k=int(np.argmax(np.abs(M.ravel())>1e-9));ph=M.ravel()[k]/abs(M.ravel()[k]);N=M/ph
        return tuple(np.round(np.r_[N.real.ravel(),N.imag.ravel()],9))
    D={key(np.eye(3)):np.eye(3,dtype=complex)};front=list(D.values())
    while front:
        new=[]
        for M in front:
            for g in gens:
                N=g@M;k=key(N)
                if k not in D:D[k]=N;new.append(N)
        front=new
    assert len(D)==216
    return list(D.values()),T
def orbit_certificate():
    G,stab_states=projective_group()
    zz=np.exp(2j*np.pi/9)
    v=np.array([1,zz,zz.conjugate()],complex)/np.sqrt(3)
    def rkey(q):
        j=int(np.argmax(np.abs(q)>1e-9));ph=q[j]/abs(q[j]);q=q/ph
        return tuple(np.round(np.r_[q.real,q.imag],9))
    rays={}
    for M in G:rays.setdefault(rkey(M@v),M@v)
    assert len(rays)==72
    stab=sum(abs(abs(np.vdot(v,M@v))-1)<1e-8 for M in G)
    assert stab==3
    comps=json.loads((ROOT/"data/w33_pass11262_g26_mub_s4_yukawa_bridge.json").read_text())["mub_action"]["components"]
    axes=[0,0,0,0]
    for q in rays.values():
        pur=[]
        for c in comps:
            p=np.abs(np.array(stab_states)[c].conj()@q)**2
            pur.append(float(np.sum(p*p)))
        hit=[i for i,x in enumerate(pur) if abs(x-1/3)<1e-8]
        assert len(hit)==1
        assert all(abs(pur[i]-5/9)<1e-8 for i in range(4) if i!=hit[0])
        axes[hit[0]]+=1
    assert axes==[18,18,18,18]
    return {"projective_clifford_order":216,"T_state_orbit_size":72,
            "T_state_stabilizer_order":3,"MUB_axis_counts":axes,
            "MUB_purity_pattern":["1/3","5/9","5/9","5/9"],
            "TM1_reading":"one MUB/tetrahedral phi axis is selected; the residual C3 cycles the other three axes"}
def payload():
    ex=exact_local_certificate()
    kin=kinetic_certificate()
    orb=orbit_certificate()
    out={
      "schema":"w33.20261001.balanced-g26-tmagic-vacuum.v1",
      "status":"PASS_BALANCED_21_MIRROR_BARRIER_HAS_EXACT_T_MAGIC_LOCAL_MINIMUM_AND_SELECTS_ONE_MUB_AXIS",
      "potential":{
        "balanced_projective_barrier":"V=-sum_9SIC log|<phi|v>|^2-sum_12MUB log|<s|v>|^2, ||v||=1",
        "divisor_origin":"(A*B)^5 is proportional to P*J^2, so all 21 mirrors receive equal logarithmic weight",
        "polynomial_wall_no_go":"V24=lambdaS*||v||^6*|A|^2+lambdaM*|B|^2 has zero-energy SIC-intersection-MUB lines for lambdaS,lambdaM>=0; it cannot select a regular vacuum."
      },
      "canonical_vacuum":"|T3>=(1,zeta9,zeta9^-1)/sqrt(3)",
      "exact_local":ex,
      "kinetic_metric":kin,
      "finite_orbit":orb,
      "repo_crosscheck":{
        "Pass416":"the same 72-ray T-state Clifford orbit was already exhausted as the input/output correction family in the five-qutrit distillation search",
        "qutrit_T_port":"analysis/w33_qutrit_t_teleportation_port.py uses exactly this T_PLUS ray as its Hesse-SIC factory audit reference",
        "port_boundary":"the port injects a separate T-Choi resource; selecting T_PLUS does not by itself prepare that entangled Choi pair",
        "increment":"new here is dynamical selection of that known resource orbit by the G26 mirror chamber plus its one-MUB/tetrahedral-axis signature"
      },
      "parents":[
        "data/w33_20261001_g26_e8_mirror_separation.json",
        "data/w33_pass11262_g26_mub_s4_yukawa_bridge.json",
        "data/w33_pass416_qutrit_distillation_search.json"
      ],
      "time_reversal_weld":{
        "largest_probability":"((1+2*cos(2*pi/9))/3)^2 = 0.712386014201...",
        "meaning":"the same cyclotomic level appears as the exact two-qutrit reversal fidelity level in Pass11265",
        "claim":"shared algebraic level, not by itself a dynamical identity"
      },
      "chirality_firewall":{
        "parity_map":"swap computational basis states 1 and 2 sends |T3> to its complex conjugate and is Clifford",
        "quotient_coordinates":"u6=u12=0, u18=-1/729 are therefore identical for both conjugates",
        "consequence":"no scalar G26-invariant Cartan potential can choose the CP/time orientation inside this orbit"
      },
      "tm1_bridge":{
        "selected_object":"exactly one of four MUBs, hence one tetrahedral phi axis",
        "residual":"C3 vertex stabilizer cycles the remaining three MUBs",
        "next_selector":"one ternary choice among the remaining three axes gives an ordered MUB pair and hence an incident oriented chord; Pass11267 classifies every such incident chord as TM1"
      },
      "boundary":"The local minimum and Hessian are exact; the 72-ray orbit/MUB census is finite numerical group enumeration. Global optimality of the log barrier is not proved. A physical radial scale, chord-selector dynamics and chirality-breaking interface remain separate requirements.",
      "checks":{"exact_stationary":True,"hessian_positive":[30,30,54,54],"quotient_coordinates_exact":True,
                "kinetic_pullback_exact":True,"orbit72":True,"one_MUB_axis_per_ray":True}
    }
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true")
    a=ap.parse_args();p=payload();txt=json.dumps(p,indent=2,sort_keys=True,default=str)+"\n"
    if a.check:
        if not OUT.exists() or OUT.read_text()!=txt:raise SystemExit("certificate drift")
    else:OUT.write_text(txt)
    print(json.dumps({"status":p["status"],"hessian":p["exact_local"]["hessian_eigenvalues"],"orbit":72}))

if __name__=="__main__":
    raise SystemExit(main())
