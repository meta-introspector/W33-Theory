"""Pass 11501: six qutrits, one magic gate, extended to the n = 5 precision.

Continues Pass 11489's streamed, checkpointed census (same seeds, same interleaving: the first 1000 classes are 11489's)
to 2000 classes per kind, and writes its own certificate.  Pass 11489's committed certificate is left unchanged.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11437_four_qutrit_cells as C  # noqa: E402
import w33_pass11489_six_qutrits as S  # noqa: E402

OUT = ROOT / "data" / "w33_pass11501_six_qutrits_extended.json"


def summarise(per_kind):
    done = json.load(open(S.CKPT))
    keep = {k: v for k, v in done.items() if int(k.split("|")[0]) - S.SEED0 < per_kind}
    cells = defaultdict(Counter)
    for r in keep.values():
        cells[r["cell"]][r["status"]] += 1
    cells = {k: dict(v) for k, v in cells.items()}
    est = C.estimate(6, cells)
    p, s = est["P_estimate"], est["stderr"]
    res = dict(pass_id=11501, classes=len(keep), completed_by_kind=dict(Counter(k.split("|")[1] for k in keep)),
               sampled_cells=cells, estimate=est,
               z_vs_one_eighth=(p - 0.125) / s, z_vs_n5=(p - 0.13824) / (s ** 2 + 0.00405 ** 2) ** 0.5,
               sequence={"1": 0.125, "2": 0.125, "3": 437 / 3276, "4": [0.1284, 0.0045], "5": [0.13824, 0.00405],
                         "6": [p, s]})
    return res


def main():
    per_kind = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    if "--summarise-only" not in sys.argv:
        S.run(per_kind, nproc)
    res = summarise(per_kind)
    print(json.dumps(res, indent=1, default=str), flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
