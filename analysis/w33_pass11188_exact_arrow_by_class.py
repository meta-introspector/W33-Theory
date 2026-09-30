"""Pass 11188: the intrinsic arrow of time, exactly, on every conjugacy class of PSp(4,3) and PSp(6,3).

Pass 11183 defined the intrinsic arrow of a Clifford tick S as

    A(S) = min over all stabilizer tensor factorisations F of  sum_q export_q(S, F),
    export_q(S, F) = rank S_F[others, q]  (= 2 - I(q_in : q_out) in trits),

and found A in {0, 2} exactly for two qutrits and A in {0, 2, 3} on a 240-element sample for three qutrits.

A is a class function: conjugating S by g carries the factorisation F to gF and S_{gF}(gSg^-1) = S_F(S), and the sign
-S changes no rank.  So A is a function on the conjugacy classes of PSp(2n,3).  GAP (analysis/gap/
w33_pass11188_class_reps.g) prints one matrix per class in an adapted symplectic basis together with the class size
and the number of factorisations the class fixes; here every representative is

  * checked to be symplectic for our form and of the GAP projective order,
  * checked to be local in exactly GAP's number of fixed factorisations (the Pass 11180 cross-check),
  * scanned over all 45 (n = 2) or 110565 (n = 3) factorisations for its exports.

This turns the three-qutrit sample of Pass 11183 into an exact distribution with class sizes, not a sample.
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11180_mereology as M  # noqa: E402
import w33_pass11182_paper_ticks_mereology as PT  # noqa: E402
import w33_pass11183_intrinsic_arrow as AR  # noqa: E402

GAPOUT = ROOT / "data" / "w33_pass11188_gap_class_reps.txt"
OUT = ROOT / "data" / "w33_pass11188_exact_arrow_by_class.json"
ORDER = {2: 25920, 3: 4585351680}


def parse_classes(path=GAPOUT):
    text = path.read_text().replace("\\\n", "")
    text = re.sub(r"\s+", " ", text)
    out = {2: [], 3: []}
    pat = r"CLASS n=(\d) order=(\d+) size=(\d+) fixed=(\d+) M=\s*(\[ \[.*?\] \])"
    for n, order, size, fixed, mat in re.findall(pat, text):
        rows = re.findall(r"\[([^\[\]]*)\]", mat)
        Mg = np.array([[int(x) for x in r.split(",")] for r in rows], np.int64)
        out[int(n)].append(dict(order=int(order), size=int(size), fixed=int(fixed), S=Mg.T % 3))
    return out


def facs(n):
    return np.array(M.factorisations(2)) if n == 2 else PT.load_facs()


def run():
    classes = parse_classes()
    res = {}
    for n in (2, 3):
        Bs = facs(n)
        J = M.form(n)
        rows = []
        for c in classes[n]:
            S = c["S"]
            assert np.array_equal((S.T @ J @ S) % 3, J % 3), "not symplectic in our basis"
            assert M.porder(S, n) == c["order"], (M.porder(S, n), c["order"])
            prof = M.profile(S, Bs, n)
            assert prof["local"] == c["fixed"], (prof["local"], c["fixed"])
            ex = AR.exports(S, Bs, n)
            tot = ex.sum(1)
            A = int(tot.min())
            rows.append(dict(order=c["order"], size=c["size"], fixed=c["fixed"], perfect=prof["perfect"], A=A,
                             argmin=int((tot == A).sum()), max_export=int(tot.max()),
                             export_hist={int(k): int(v) for k, v in sorted(Counter(tot.tolist()).items())}))
        assert sum(r["size"] for r in rows) == ORDER[n]
        dist = Counter()
        for r in rows:
            dist[r["A"]] += r["size"]
        res[n] = dict(classes=rows,
                      A_distribution={int(a): dict(elements=int(m), fraction=str(Fraction(m, ORDER[n])))
                                      for a, m in sorted(dist.items())})
    return res


def main():
    res = run()
    for n in (2, 3):
        print(f"n={n}: A distribution over PSp({2*n},3):", res[n]["A_distribution"])
    json.dump({str(k): v for k, v in res.items()}, open(OUT, "w"), indent=1)
    print("wrote", OUT.name)


if __name__ == "__main__":
    main()
