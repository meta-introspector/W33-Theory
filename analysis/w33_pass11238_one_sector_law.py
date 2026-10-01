"""Pass 11238: is the one-sector law exact? Every word up to length 5 over {cubic phase on one qutrit, SUM, SUM^dag}.

Pass 11227 found that (T(x)T) SUM breaks substrate time-reversal symmetry while every single-sector placement tried did
not.  This pass tests the statement exhaustively on finite families (words of length <= L, deduplicated as matrices up
to phase), with the exact best time-reversal fidelity over all 51840 x 81 anti-unitary two-qutrit Cliffords:
  * target sector: gates {I(x)T, SUM, SUM^dag};
  * control sector: gates {T(x)I, SUM, SUM^dag};
  * both sectors (control family): {T(x)I, I(x)T, SUM, SUM^dag}.
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11227_two_qutrit_cp_mixing as P2  # noqa: E402
import w33_pass11235_two_qutrit_t_spectrum as S2  # noqa: E402

OUT = ROOT / "data" / "w33_pass11238_one_sector_law.json"
QUANTUM = (1 + 2 * np.cos(2 * np.pi / 9)) / 3


def key(U):
    f = U.ravel()
    k = int(np.argmax(np.abs(f) > 1e-9))
    V = U * (abs(f[k]) / f[k])
    return tuple(np.round(V.real, 6).ravel() + 0.0) + tuple(np.round(V.imag, 6).ravel() + 0.0)


def words(gates, L):
    seen = {}
    for n in range(1, L + 1):
        for combo in itertools.product(range(len(gates)), repeat=n):
            U = np.eye(9, dtype=complex)
            for g in combo:
                U = gates[g][1] @ U
            if not any("SUM" not in gates[g][0] for g in combo):       # at least one cubic gate
                continue
            k = key(U)
            if k not in seen:
                seen[k] = (combo, U)
    return list(seen.values())


def run(L=5):
    reps = np.load(P2.CACHE)
    T1 = np.kron(P2.T, P2.I3)
    T2 = np.kron(P2.I3, P2.T)
    SUM, SUMd = P2.SUM, P2.SUM.conj().T
    fams = {
        "target_sector": [("I(x)T", T2), ("SUM", SUM), ("SUM^dag", SUMd)],
        "control_sector": [("T(x)I", T1), ("SUM", SUM), ("SUM^dag", SUMd)],
        "both_sectors": [("T(x)I", T1), ("I(x)T", T2), ("SUM", SUM), ("SUM^dag", SUMd)],
    }
    res = dict(pass_id=11238, max_length=L)
    for name, gates in fams.items():
        Ws = words(gates, L if name != "both_sectors" else 4)
        verdicts, fid, examples = Counter(), Counter(), []
        for combo, U in Ws:
            f = S2.fidelity_fast(U, reps)
            v = "reversible" if f > 1 - 1e-9 else "violating"
            verdicts[v] += 1
            if v == "violating":
                fid[round(f, 9)] += 1
                if len(examples) < 5:
                    examples.append(" ".join(gates[g][0] for g in combo))
        res[name] = dict(distinct_words=len(Ws), verdicts=dict(verdicts),
                         violator_fidelities={str(k): v for k, v in fid.items()}, violating_examples=examples)
        print(name, res[name], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
