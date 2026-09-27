#!/usr/bin/env python3
"""Pass 11040: objectwise CCZ finite-difference / E8 cubic-bracket temporal echo."""
from __future__ import annotations
import argparse, itertools, json, sys
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import w33_20260925_cubic_clock_bracket_ladder as E8

OUT=ROOT/"data/w33_pass11040_exceptional_cubic_temporal_echo.json"


def delta(vals,axis):
    out={}
    for x in itertools.product(range(3),repeat=3):
        y=list(x); y[axis]=(y[axis]+1)%3
        out[x]=(vals[tuple(y)]-vals[x])%3
    return out


def lift(z):
    z=int(z)%3
    return 0 if z==0 else (1 if z==1 else -1)


def cubic_echo(sign):
    vals={x:(sign*x[0]*x[1]*x[2])%3
          for x in itertools.product(range(3),repeat=3)}
    for axis in (0,1,2):
        vals=delta(vals,axis)
    assert len(set(vals.values()))==1
    return lift(next(iter(vals.values())))
def payload():
    cubic_path=ROOT/"extracted_v13/W33-Theory-master/artifacts/canonical_su3_gauge_and_cubic.json"
    raw=json.loads(cubic_path.read_text())
    triads={tuple(r["triple"]):int(r["sign"]) for r in raw["solution"]["d_triples"]}
    table,_=E8.pair_cross_table(triads)
    assert len(triads)==45
    echo_by_sign={s:cubic_echo(s) for s in (-1,1)}
    assert echo_by_sign=={-1:-1,1:1}

    checks_count=0
    nonzero=0
    reversal=0
    mismatch=0
    rows=[]
    for tri,d in sorted(triads.items()):
        i,j,k=tri
        gate=cubic_echo(d)
        assert gate==d
        for a,b,c in itertools.product(range(2),repeat=3):
            eps=int(E8.EPS[a,b])
            bracket=E8.nested_left(i,a,j,b,k,c,table)
            expected=np.zeros(2,dtype=int)
            expected[c]=-eps*gate
            if not np.array_equal(bracket,expected):
                mismatch+=1
            if eps:
                nonzero+=1
                rev=E8.nested_left(i,b,j,a,k,c,table)
                assert np.array_equal(rev,-bracket)
                reversal+=1
            checks_count+=1
        rows.append({"triad":list(tri),"cubic_sign":d,"gate_echo":gate})
    assert mismatch==0
    assert checks_count==45*8
    assert nonzero==45*4
    assert reversal==nonzero

    # Bare scalar finite differences commute: orientation does NOT come from
    # permuting the three delta operators.
    orientation_blind=0
    for d in (-1,1):
        base={x:(d*x[0]*x[1]*x[2])%3
              for x in itertools.product(range(3),repeat=3)}
        vals=[]
        for order in itertools.permutations(range(3)):
            z=base
            for ax in order:z=delta(z,ax)
            vals.append(lift(next(iter(z.values()))))
        assert len(set(vals))==1 and vals[0]==d
        orientation_blind+=1
    checks={
      "45_signed_E6_triads":len(triads)==45,
      "triple_finite_difference_returns_cubic_sign":echo_by_sign=={-1:-1,1:1},
      "all_360_objectwise_matches":mismatch==0 and checks_count==360,
      "180_nonzero_oriented_echoes":nonzero==180,
      "temporal_basis_swap_reverses_nested_bracket":reversal==180,
      "bare_finite_difference_order_is_orientation_blind":orientation_blind==2,
    }
    assert all(checks.values())
    return {
      "schema":"w33.pass11040.exceptional-cubic-temporal-echo.v1",
      "status":"PASS",
      "headline":(
        "The qutrit cubic echo and the committed E8 bracket ladder now match "
        "objectwise on all 45 signed E6 triads. For f_ijk(x,y,z)=d_ijk xyz, "
        "Delta_x Delta_y Delta_z f_ijk=d_ijk, while the E8 nested bracket is "
        "[[e_(i,a),e_(j,b)],e_(k,c)]=-eps(a,b)d_ijk T_c. All 360 basis cases "
        "match the single formula -eps(a,b)*(triple finite difference)*T_c."
      ),
      "gate_level":{
        "phase":"omega^(d_ijk*x*y*z)",
        "triple_forward_difference":"d_ijk mod 3",
        "signed_echo_values":{str(k):v for k,v in echo_by_sign.items()},
        "difference_orders_tested":6,
        "orientation_from_difference_order":False,
      },
      "exceptional_level":{
        "triads":len(triads),
        "basis_cases":checks_count,
        "nonzero_cases":nonzero,
        "identity":"[[e_(i,a),e_(j,b)],e_(k,c)] = -eps(a,b) Delta_x Delta_y Delta_z(d_ijk xyz) T_c",
        "temporal_reversal":"swap a and b; eps flips sign",
      },
      "important_correction":(
        "The bare scalar triple finite difference does not itself know temporal "
        "orientation: finite-difference operators commute. The orientation sign "
        "comes from the alternating temporal two-form eps (equivalently, from the "
        "choice of conjugate temporal orientation/antiunitary phase convention). "
        "Thus a cubic temporal witness should include that antisymmetric carrier, "
        "not identify time reversal with merely permuting Delta operators."
      ),
      "triad_table":rows,
      "boundary":(
        "This proves coefficient-level and basis-level equality between the "
        "qutrit cubic derivative and the repository's abstract E8 nilpotent "
        "bracket model. It still does not derive a Hamiltonian that implements "
        "these bracket generators as laboratory interventions or fix an energy/time scale."
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
    print(json.dumps({"status":p["status"],"cases":p["exceptional_level"]["basis_cases"],
      "nonzero":p["exceptional_level"]["nonzero_cases"],
      "orientation_source":"eps"},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
