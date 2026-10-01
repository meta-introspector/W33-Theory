#!/usr/bin/env python3
"""Pass 11242: why a U(1)' that protects mu keeps exotics light -- an anomaly theorem, checked on every vacuum.

THEOREM.  Let U(1)' be an anomaly-free U(1) (no SU(3)^2-U(1)' anomaly) under which H_u H_d is charged:
t(H_u) + t(H_d) != 0.  Then the colored fields cannot all become massive through U(1)'-allowed mass terms -- vector-like
masses X Xbar <S...>, up-type Yukawas q bu H_u and bq u H_d, down-type Yukawas q bd H_d and bq d H_u.

Proof.  Split colored states by electric charge.  Up sector: triplets q_up (one per q), u; antitriplets bq_up, bu.
Down sector: q_down, d; bq_down, bd.  A complete set of masses is a perfect matching in each sector.  Every allowed
matched pair (a, b) has t_a + t_b = -t(H) with H = H_u, H_d or nothing:
    up:   (q,bu)[H_u] n1,  (q,bq)[-] n2,  (u,bu)[-] n3,  (u,bq)[H_d] n4
    down: (q,bd)[H_d] m1,  (q,bq)[-] m2,  (d,bd)[-] m3,  (d,bq)[H_u] m4.
Counting q and bq in the up sector gives n1 - n4 = #q - #bq = 3; in the down sector m1 - m4 = 3.  So the numbers of H_u-
and H_d-pairs are equal: a = n1 + m4 = b = m1 + n4 = 3 + n4 + m4 >= 3.  The SU(3)^2-U(1)' anomaly is the sum of t over
all triplet and antitriplet components = sum over matched pairs of (t_a + t_b) = -a (t(H_u) + t(H_d)).  It vanishes, so
t(H_u) + t(H_d) = 0, contradicting the hypothesis.  QED.  (Generation charges need not be universal; no assumption on
Yukawa textures beyond the perfect matching.  The classic family-universal special case -- U(1)' forbidding mu needs
exotics -- is known from U(1)'-extended MSSM model building; this is the matching form used here.)

So an unbroken non-anomalous U(1)' can protect mu only if some colored state stays massless (a light exotic or a
massless quark).  In heterotic orbifolds every U(1) except the single anomalous one is anomaly-free, and the anomalous
one is broken by the FI vacuum, so the theorem applies to every vacuum U(1)' of Pass 11232.  Discrete Z_N analogues hold
mod N up to Green-Schwarz shifts and are NOT claimed (that is how Z4^R evades it in the literature).

CHECKS (this producer): SU(3)^2-U(1) anomaly of every non-anomalous U(1) generator of all 215 models is zero; over every
FI-cancelling vacuum span of Pass 11232 and every mu-protected H_u, NO choice of H_d gives a perfect colored matching
(the theorem's prediction); the escape used (which sector fails) is tallied.
"""
from __future__ import annotations

import json
import re
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass10960_matter_even_dflat_closure as P  # noqa: E402
import w33_pass11232_mu_vacuum_symmetry as MU  # noqa: E402

OUT = ROOT / "data" / "w33_pass11242_u1prime_mu_anomaly_theorem.json"
ab = lambda s: abs(int(re.match(r"(-?\d+)", s).group(1)))


def su3_anomaly_vector(model):
    """A_3[i] = sum over colored fields of (SU(2) x hidden multiplicity) q_i, with sign +1 for 3 and -1 for 3bar
    (index 1/2 dropped).  Must vanish for every non-anomalous U(1) direction."""
    left = model["left"]
    q0 = left[[f["name"] for f in left].index(next(f["name"] for f in left if f["name"].startswith("q_")))]
    c3 = [i for i, x in enumerate(q0["dim"].split(",")) if ab(x) == 3][0]
    A = [Fraction(0)] * model["nu1"]
    colored = 0
    for f in left:
        dims = f["dim"].split(",")
        d3 = int(re.match(r"(-?\d+)", dims[c3]).group(1))
        if abs(d3) != 3:
            continue
        colored += 1
        mult = 1
        for i, x in enumerate(dims):
            if i != c3:
                mult *= ab(x)
        sgn = 1 if d3 > 0 else -1
        for i, v in enumerate(f["q"]):
            A[i] += mult * Fraction(v)        # chirality, not 3 vs 3bar, enters the anomaly with the same sign
    return A, colored


def perfect(rows, cols, ok):
    rows, cols = list(rows), list(cols)
    return len(rows) == len(cols) and MU.structural_rank(rows, cols, ok) == len(rows)


def check_model(name, model):
    left = MU.fields(model)
    by = {}
    for f in left:
        by.setdefault(f["base"], []).append(f)
    tr = [sum((f["dimprod"] * f["q"][i] for f in left), Fraction(0)) for i in range(model["nu1"])]
    rec = dict(anomalous=tr[0] != 0)
    if tr[0] == 0:
        return rec
    lab = [f for f in left if f["base"] in P.SMY]
    y = P.solve_affine([f["q"] for f in lab], [P.SMY[f["base"]] for f in lab])[0]
    sing = [f for f in left if f["trivial"] and P.dot(y, f["q"]) == 0]
    types = {}
    for f in sing:
        types.setdefault(tuple(f["q"]), []).append(f["name"])
    T = list(types)
    s = 1 if tr[0] > 0 else -1
    C = [s * t[0] for t in T]
    Mx = [list(t[1:]) for t in T]
    if P.farkas(C, Mx) is not None:
        rec["dflat"] = False
        return rec
    rays = [r for r in P.cone_rays(Mx) if P.dot(C, r) < 0] if len(T) <= 44 else MU.sampled_vertices(C, Mx)
    add = lambda *vs: [sum(x) for x in zip(*vs)]
    g = lambda b: by.get(b, [])
    hu, hd = g("bl"), g("l")
    raw = {f["name"]: f["dim"].split(",") for f in model["left"]}
    q0 = raw[g("q")[0]["name"]]
    l0 = raw[g("l")[0]["name"]]
    c3 = [i for i, x in enumerate(q0) if ab(x) == 3][0]
    c2 = [i for i, x in enumerate(l0) if ab(x) == 2 and ab(q0[i]) == 2][0]

    def hkey(f):
        """hidden-sector representation up to conjugation (permissive: a necessary condition for a mass pairing)"""
        return tuple(ab(x) for i, x in enumerate(raw[f["name"]]) if i not in (c3, c2))

    def comps(b):
        out = []
        for f in g(b):
            m = 1
            for x in hkey(f):
                m *= x
            out += [(b, f)] * m
        return out
    seen, prot, violations, escapes = set(), 0, 0, {"up sector": 0, "down sector": 0, "both": 0}
    for r in rays:
        sup = [T[k] for k, v in enumerate(r) if v != 0]
        V = MU.Span([list(t) for t in sup])
        key = V.key()
        if key in seen:
            continue
        seen.add(key)
        inV = lambda *fs: add(*[f["q"] for f in fs]) in V
        for H in hu:
            if any(inV(H, L) for L in hd):
                continue
            if not any(inV(qa, H, ub) for qa in g("q") for ub in g("bu")):
                continue
            prot += 1
            up_rows = comps("q") + comps("u")
            up_cols = comps("bq") + comps("bu")
            dn_rows = comps("q") + comps("d")
            dn_cols = comps("bq") + comps("bd")
            up_ok_any = dn_ok_any = False
            for Hd in hd:
                def up(i, j, H=H, Hd=Hd):
                    (ta, fa), (tb, fb) = up_rows[i], up_cols[j]
                    if hkey(fa) != hkey(fb):
                        return False
                    if (ta, tb) == ("q", "bu"):
                        return inV(fa, fb, H)
                    if (ta, tb) == ("u", "bq"):
                        return inV(fa, fb, Hd)
                    return inV(fa, fb)

                def dn(i, j, H=H, Hd=Hd):
                    (ta, fa), (tb, fb) = dn_rows[i], dn_cols[j]
                    if hkey(fa) != hkey(fb):
                        return False
                    if (ta, tb) == ("q", "bd"):
                        return inV(fa, fb, Hd)
                    if (ta, tb) == ("d", "bq"):
                        return inV(fa, fb, H)
                    return inV(fa, fb)
                pu = perfect(range(len(up_rows)), range(len(up_cols)), up)
                pd = perfect(range(len(dn_rows)), range(len(dn_cols)), dn)
                up_ok_any |= pu
                dn_ok_any |= pd
                if pu and pd:
                    violations += 1
            escapes["both" if not up_ok_any and not dn_ok_any else ("up sector" if not up_ok_any else
                                                                      ("down sector" if not dn_ok_any else "both"))] += 1
    rec.update(dflat=True, spans=len(seen), protected=prot, theorem_violations=violations, escapes=escapes,
               component_counts=dict(up=[len(comps("q")) + len(comps("u")), len(comps("bq")) + len(comps("bu"))],
                                     down=[len(comps("q")) + len(comps("d")), len(comps("bq")) + len(comps("bd"))]))
    return rec


def run():
    ledger, sha = P.load_ledger()
    res, anom = {}, {}
    t0 = time.time()
    for n in sorted(ledger):
        m = ledger[n]
        A, nc = su3_anomaly_vector(m)
        anom[n] = dict(A3=[str(a) for a in A], nonanomalous_zero=all(a == 0 for a in A[1:]))
        res[n] = check_model(n, m)
        print(n, res[n], flush=True)
    summary = dict(models=len(res), su3_anomaly_zero_on_all_nonanomalous_U1s=all(v["nonanomalous_zero"] for v in
                                                                                  anom.values()),
                   protected_cases=sum(v.get("protected", 0) for v in res.values()),
                   theorem_violations=sum(v.get("theorem_violations", 0) for v in res.values()),
                   escapes={k: sum(v.get("escapes", {}).get(k, 0) for v in res.values())
                            for k in ("up sector", "down sector", "both")},
                   ledger_sha256=sha, seconds=round(time.time() - t0, 1))
    return dict(pass_id=11242, summary=summary, su3_anomaly=anom, models=res)


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res["summary"], indent=1))


if __name__ == "__main__":
    main()
