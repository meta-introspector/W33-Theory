"""Pass 11514: six qutrits, one magic gate, to 4000 classes (the n = 5 precision).

Continues the streamed, checkpointed census of Passes 11489/11501 (same seeds, same interleaving; 2500 classes already done)
to 2000 classes per kind, and writes its own certificate and a frozen copy of the per-class checkpoint.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11489_six_qutrits as S  # noqa: E402
import w33_pass11501_six_qutrits_extended as X  # noqa: E402

OUT = ROOT / "data" / "w33_pass11514_six_qutrits_4000.json"
FROZEN = ROOT / "data" / "w33_pass11514_checkpoint.json"


def main():
    per_kind = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    if "--summarise-only" not in sys.argv:
        S.run(per_kind, nproc)
    res = X.summarise(per_kind)
    res["pass_id"] = 11514
    shutil.copyfile(S.CKPT, FROZEN)
    print(json.dumps(res, indent=1, default=str), flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
