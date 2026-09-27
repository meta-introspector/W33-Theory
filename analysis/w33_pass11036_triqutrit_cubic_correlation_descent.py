#!/usr/bin/env python3
"""Pass 11036: exact triqutrit cubic derivative ladder and correlation-sector firewall."""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11036_triqutrit_cubic_correlation_descent.json"
F=range(3)


def f(a,b,c):
    return (a*b*c)%3


def delta(vals,axis):
    out={}
    for x in itertools.product(F,repeat=3):
        y=list(x); y[axis]=(y[axis]+1)%3
        out[x]=(vals[tuple(y)]-vals[x])%3
    return out


def depends_on(vals):
    active=[]
    for ax in range(3):
        found=False
        for x in itertools.product(F,repeat=3):
            y=list(x); y[ax]=(y[ax]+1)%3
            if vals[tuple(y)]!=vals[x]:
                found=True; break
        if found: active.append(ax)
    return active
def contraction_rank():
    # U_i is 2D. Fix one nonzero covector lambda_i=(1,0) on each.
    # The direct-sum map T -> (lambda1 T, lambda2 T, lambda3 T)
    # goes C3 (dim 8) -> three pair blocks (dim 12).
    triples=list(itertools.product(range(2),repeat=3))
    pair12=[(ax,p,q) for ax in range(3)
            for p,q in itertools.product(range(2),repeat=2)]
    M=sp.zeros(12,8)
    row={x:i for i,x in enumerate(pair12)}
    for col,(i,j,k) in enumerate(triples):
        if i==0: M[row[(0,j,k)],col]=1
        if j==0: M[row[(1,i,k)],col]=1
        if k==0: M[row[(2,i,j)],col]=1
    r=M.rank()
    ker=M.nullspace()
    assert r==7 and len(ker)==1
    # Each individual contraction is surjective onto its 4D pair block.
    block_ranks=[M[4*a:4*a+4,:].rank() for a in range(3)]
    assert block_ranks==[4,4,4]
    return r,block_ranks,[list(map(int,v)) for v in ker]


def payload():
    v0={x:f(*x) for x in itertools.product(F,repeat=3)}
    v1=delta(v0,0)
    v2=delta(v1,1)
    v3=delta(v2,2)
    # abc -> bc -> c -> 1, with the chosen forward-difference convention.
    assert all(v1[x]==x[1]*x[2]%3 for x in v1)
    assert all(v2[x]==x[2]%3 for x in v2)
    assert set(v3.values())=={1}
    deps=[depends_on(v) for v in (v0,v1,v2,v3)]
    assert deps==[[0,1,2],[1,2],[2],[]]
    dims={"C0":1,"C1":6,"C2":12,"C3":8}
    assert sum(dims.values())==27
    cr,blocks,ker=contraction_rank()

    checks={
        "ccz_cubic_phase_27_histories":len(v0)==27,
        "first_difference_is_pairwise_bc":all(v1[x]==x[1]*x[2]%3 for x in v1),
        "second_difference_is_local_c":all(v2[x]==x[2]%3 for x in v2),
        "third_difference_is_central_one":set(v3.values())=={1},
        "correlation_dimensions_1_6_12_8":dims=={"C0":1,"C1":6,"C2":12,"C3":8},
        "each_slant_block_rank4":blocks==[4,4,4],
        "fixed_three_slants_joint_rank7_not12":cr==7,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11036.triqutrit-cubic-correlation-descent.v1",
      "status":"PASS",
      "headline":(
        "The qutrit cubic phase f(a,b,c)=abc has the exact finite-difference "
        "ladder abc -> bc -> c -> 1, i.e. CCZ -> controlled phase -> local "
        "phase -> central ternary phase. Separately, a 3x3x3 carrier splits "
        "by correlation order as 1+6+12+8. The crucial correction is that "
        "8->12 cannot be a single surjective linear differential: three fixed "
        "slant contractions jointly have rank 7, although each separately "
        "surjects onto its own 4D pair block."
      ),
      "cubic_echo":{
        "phase":"omega^(a*b*c)",
        "successive_exponents":["abc","bc","c","1"],
        "active_coordinate_counts":[3,2,1,0],
        "operator_reading":["CCZ","CZ-like","Z-like","central omega I"],
      },
      "correlation_decomposition":{
        "dimensions":dims,
        "identity":"27 = 1 + 6 + 12 + 8",
        "pairwise_blocks":"C2 = (U1 tensor U2) + (U1 tensor U3) + (U2 tensor U3), 4+4+4",
      },
      "slant_rank_audit":{
        "source_dimension":8,
        "target_direct_sum_dimension":12,
        "individual_pair_block_ranks":blocks,
        "joint_fixed_three_slants_rank":cr,
        "joint_kernel_dimension":1,
        "kernel_basis":ker,
      },
      "missing_12_firewall":(
        "The numbers 8,12,6,1 are correlation-sector dimensions, not the "
        "ranks of a surjective chain complex. A single cubic tensor can have "
        "nonzero projections into all three pairwise blocks, but three fixed "
        "contractions carry only seven independent combinations into their "
        "12D direct sum. Therefore the claim that a diagonal weld restores "
        "the entire missing 12 requires extra background degrees of freedom "
        "or an explicit repository intertwiner; it does not follow from CCZ "
        "or the dimension ladder alone."
      ),
      "boundary":(
        "This establishes the finite gate identity and the linear-algebra "
        "firewall. It does not yet identify the E6 diagonal phase-weld fields "
        "with a basis of all three 4D pairwise sectors."
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
        if not a.output.exists() or a.output.read_text()!=text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],
        "ladder":p["cubic_echo"]["successive_exponents"],
        "slant_rank":p["slant_rank_audit"]["joint_fixed_three_slants_rank"]},
        sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
