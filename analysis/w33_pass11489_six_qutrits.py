"""Pass 11489: six qutrits, one magic gate -- the bad-class fraction at n = 6.

Same estimator as Passes 11437/11460 (exact cell shares x sampled cell fractions), with Pass 11421's sparse decider at
n = 6: 729-dimensional states, 531,441 frames per class.  Sequence so far: P_1 = P_2 = 1/8, P_3 = 437/3276 (exact),
P_4 = 0.128 +- 0.005, P_5 = 0.138 +- 0.004.

Streams the class jobs (imap_unordered) and checkpoints every 100 classes, so a long run is never lost; rerunning resumes
from the checkpoint.  A first unstreamed attempt (6 workers, 5000 classes) was stopped after ~52,000 CPU-seconds without
output because it starved the other computations of memory; nothing from it is used.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11437_four_qutrit_cells as C  # noqa: E402

OUT = ROOT / "data" / "w33_pass11489_six_qutrits.json"
CKPT = ROOT / "data" / "w33_pass11489_six_qutrits_checkpoint.json"
SEED0 = 114890000


def run(per_kind=1200, nproc=4):
    done = json.load(open(CKPT)) if CKPT.exists() else {}
    jobs = [(SEED0 + i, kind) for i in range(per_kind) for kind in ("collinear", "non-collinear")]   # interleaved
    todo = [j for j in jobs if f"{j[0]}|{j[1]}" not in done]
    print(f"{len(done)} done, {len(todo)} to do", flush=True)
    if not todo:
        nproc = 1
    with Pool(nproc, initializer=C._init, initargs=(6,)) as pool:
        for n, (job, r) in enumerate(zip(todo, pool.imap(C._job, todo, chunksize=2)), 1):
            done[f"{job[0]}|{job[1]}"] = dict(cell=r["cell"], status="undecided" if r.get("undecided") else
                                              ("bad" if r["bad"] else "good"))
            if n % 100 == 0:
                json.dump(done, open(CKPT, "w"))
                print(n, "classes", flush=True)
    json.dump(done, open(CKPT, "w"))
    cells = defaultdict(Counter)
    for key, r in done.items():                          # every completed class counts (interleaving keeps kinds balanced)
        cells[r["cell"]][r["status"]] += 1
    cells = {k: dict(v) for k, v in cells.items()}
    kinds = Counter(k.split("|")[1] for k in done)
    res = dict(pass_id=11489, completed_by_kind=dict(kinds), n6=dict(sampled_cells=cells, estimate=C.estimate(6, cells)))
    p, s = res["n6"]["estimate"]["P_estimate"], res["n6"]["estimate"]["stderr"]
    res["sequence"] = {"1": 0.125, "2": 0.125, "3": 437 / 3276, "4": [0.1284, 0.0045], "5": [0.13824, 0.00405],
                       "6": [p, s]}
    res["n6_z_vs_one_eighth"] = (p - 0.125) / s
    res["n6_z_vs_n5"] = (p - 0.13824) / (s ** 2 + 0.00405 ** 2) ** 0.5
    print(json.dumps(res, indent=1, default=str), flush=True)
    return res


def main():
    per_kind = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    res = run(per_kind, nproc)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
