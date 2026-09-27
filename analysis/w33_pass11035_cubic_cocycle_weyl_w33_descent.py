#!/usr/bin/env python3
"""Pass 11035: cubic temporal cocycle -> Weyl carrier -> doubled W33."""
from __future__ import annotations
import argparse, itertools, json, sys
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_pass10970_projective_clock_code_lattice_tower as P70

OUT=ROOT/"data/w33_pass11035_cubic_cocycle_weyl_w33_descent.json"
F=range(3)
E=((1,0,0),(0,1,0),(0,0,1))


def det3(a,b,c):
    return (
        a[0]*(b[1]*c[2]-b[2]*c[1])
        - a[1]*(b[0]*c[2]-b[2]*c[0])
        + a[2]*(b[0]*c[1]-b[1]*c[0])
    )%3


def add(a,b):
    return tuple((x+y)%3 for x,y in zip(a,b))


def points3():
    return tuple(itertools.product(F, repeat=3))
def cocycle_check():
    V=points3()
    bad=0
    for g in V:
      for h in V:
       for k in V:
        for l in V:
         x=(det3(h,k,l)+det3(g,add(h,k),l)+det3(g,h,k)
            -det3(add(g,h),k,l)-det3(g,h,add(k,l)))%3
         bad += x!=0
    return bad


def complement(i):
    return [E[j] for j in range(3) if j!=i]


def slant_packet(i):
    t=E[i]
    b=complement(i)
    H=[((a*b[0][0]+c*b[1][0])%3,
        (a*b[0][1]+c*b[1][1])%3,
        (a*b[0][2]+c*b[1][2])%3) for a in F for c in F]
    form=[[det3(t,h,k) for k in H] for h in H]
    radical=[h for h,row in zip(H,form) if all(x==0 for x in row)]
    assert radical==[(0,0,0)]
    return {"contracted_axis":i,"plane_size":9,"radical_size":len(radical),
            "basis_form":[[det3(t,b[r],b[c]) for c in range(2)] for r in range(2)]}
def normproj(v):
    for x in v:
        if x%3:
            inv=pow(x,-1,3)
            return tuple((inv*y)%3 for y in v)
    raise ValueError


def omega4(x,y):
    # future plane minus past plane
    return (x[0]*y[1]-x[1]*y[0]-x[2]*y[3]+x[3]*y[2])%3


def w33_packet():
    pts=sorted({normproj(v) for v in itertools.product(F,repeat=4)
                if any(v)})
    assert len(pts)==40
    lines=set()
    for i,a in enumerate(pts):
      for b in pts[i+1:]:
        if b==a or omega4(a,b): continue
        span={normproj(tuple((r*x+s*y)%3 for x,y in zip(a,b)))
              for r,s in itertools.product(F,repeat=2) if r or s}
        if len(span)==4: lines.add(frozenset(span))
    assert len(lines)==40
    inc={p:sum(p in L for L in lines) for p in pts}
    assert set(inc.values())=={4}
    delta=frozenset(normproj((u[0],u[1],u[0],u[1]))
                    for u in itertools.product(F,repeat=2) if any(u))
    assert len(delta)==4 and delta in lines
    transverse=[L for L in lines if L.isdisjoint(delta)]
    assert len(transverse)==27
    return {"projective_points":40,"lagrangian_lines":40,
            "points_per_line":4,"lines_per_point":4,
            "bell_diagonal_points":4,
            "transverse_lagrangians":len(transverse)}
def weyl_packet():
    w=sp.exp(2*sp.pi*sp.I/3)
    X=sp.zeros(3); Z=sp.zeros(3)
    for j in range(3):
        X[(j+1)%3,j]=1
        Z[j,j]=w**j
    assert sp.simplify(Z*X-w*X*Z)==sp.zeros(3)
    mats=[]
    for a,b in itertools.product(F,repeat=2):
        mats.append(X**a*Z**b)
    vecs=[sp.Matrix(M).reshape(9,1) for M in mats]
    r=sp.Matrix.hstack(*vecs).rank()
    assert r==9
    return {"weyl_relation":"ZX = omega XZ","monomial_span_dimension":r,
            "twisted_algebra":"M3(C)"}


def tetracode_packet():
    C=P70.codewords(3)
    nz=[c for c in C if any(c)]
    assert len(nz)==8 and set(map(P70.hamming_weight,nz))=={3}
    omitted={str(i):[] for i in range(4)}
    for c in nz:
        z=[i for i,x in enumerate(c) if x==0]
        assert len(z)==1
        omitted[str(z[0])].append(list(c))
    assert all(len(v)==2 for v in omitted.values())
    assert all(tuple((-x)%3 for x in v[0])==tuple(v[1]) or
               tuple((-x)%3 for x in v[1])==tuple(v[0])
               for v in omitted.values())
    return {"nonzero_words":8,"weight":3,"omission_classes":omitted,
            "law":"4 omitted clocks x 2 opposite codewords"}
def payload():
    bad=cocycle_check()
    assert bad==0
    slants=[slant_packet(i) for i in range(3)]
    w33=w33_packet(); weyl=weyl_packet(); tetra=tetracode_packet()
    checks={
      "determinant_phase_is_3cocycle":bad==0,
      "all_three_slants_nondegenerate":all(x["radical_size"]==1 for x in slants),
      "weyl_algebra_full_M3":weyl["monomial_span_dimension"]==9,
      "doubled_geometry_W33_counts":w33["projective_points"]==w33["lagrangian_lines"]==40,
      "bell_diagonal_lagrangian":w33["bell_diagonal_points"]==4,
      "twenty_seven_transverse_lagrangians":w33["transverse_lagrangians"]==27,
      "tetracode_four_by_two_orientation_index":tetra["nonzero_words"]==8,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11035.cubic-cocycle-weyl-w33-descent.v1",
      "status":"PASS",
      "headline":(
        "The alternating cubic phase nu(g,h,k)=omega^det(g,h,k) is an exact "
        "normalized 3-cocycle on F3^3. Contracting any one axis gives a "
        "nondegenerate 2-cocycle on the complementary F3^2 plane, hence the "
        "qutrit Weyl algebra M3(C). Doubling two such planes with opposite "
        "orientation gives exactly the 40-point/40-line W(3,3) symplectic "
        "geometry; the diagonal Bell line cancels the opposite forms and has "
        "exactly 27 transverse Lagrangian complements."
      ),
      "cubic_class":{
        "group":"F3^3","exponent":"det(g,h,k) mod 3",
        "cocycle_failures_over_all_27^4_quadruples":bad,
        "slant_contractions":slants,
        "nontriviality_witness":(
          "Each contraction has nondegenerate alternating commutator; a "
          "coboundary on an abelian plane would have trivial projective "
          "commutator, so the cubic class cannot be globally trivial."
        ),
      },
      "boundary_carrier":weyl,
      "past_future_doubling":w33,
      "tetracode_indexing":tetra,
      "interpretation":(
        "This makes the PDF's anomaly-to-qutrit idea exact after choosing the "
        "alternating determinant representative: one cubic composition class "
        "has three nondegenerate slant reductions, while opposite temporal "
        "orientations glue into the repository's native symplectic W33 label "
        "space. The tetracode supplies exactly the same 4 x 2 omission/chirality "
        "index set."
      ),
      "boundary":(
        "The 4x2 tetracode match is an exact combinatorial bijection, not yet "
        "an intertwiner between cocycle chirality and the committed E8 tensor "
        "dualizations. 'Temporal anomaly' remains an interpretation of the "
        "finite cocycle, not a derivation of physical time."
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
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text)
    print(json.dumps({"status":p["status"],"w33":p["past_future_doubling"],
                      "tetracode":p["tetracode_indexing"]["nonzero_words"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
