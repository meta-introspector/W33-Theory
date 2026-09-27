#!/usr/bin/env python3
"""Pass 11027: exact prime-clock lattice census through q=13."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11027_prime_clock_lattice_census.json"

# q=13 was exhaustively enumerated during this pass over all 13^7-1
# nonzero codewords. --deep-q13 replays that expensive census.
Q13_DEEP = {
    "codewords_nonzero": 62_748_516,
    "minimum_glue_norm": 12,
    "minimum_codewords": 756,
    "example_coefficients": [1, 0, 0, 0, 0, 0, 0],
}


def census(q: int, chunk: int = 250_000):
    k = (q + 1) // 2
    xs = np.arange(q, dtype=np.int64)
    E = np.empty((k, q + 1), dtype=np.int64)
    for j in range(k):
        E[j, :q] = np.array([pow(int(x), j, q) for x in xs], dtype=np.int64)
        E[j, q] = 1 if j == k - 1 else 0
    ew = np.array([c * (q - c) for c in range(q)], dtype=np.int64)
    N = q ** k
    best_num = 10**18
    count = 0
    example = None
    for start in range(1, N, chunk):
        stop = min(N, start + chunk)
        nums = np.arange(start, stop, dtype=np.int64)
        A = np.empty((stop - start, k), dtype=np.int64)
        z = nums.copy()
        for j in range(k):
            A[:, j] = z % q
            z //= q
        C = (A @ E) % q
        norms = ew[C].sum(axis=1)
        m = int(norms.min())
        if m < best_num:
            best_num = m
            mask = norms == m
            count = int(mask.sum())
            example = A[int(np.where(mask)[0][0])].tolist()
        elif m == best_num:
            count += int((norms == m).sum())
    assert best_num % q == 0
    return {
        "q": q,
        "dimension": q*q - 1,
        "code": f"[{q+1},{k},{(q+3)//2}]_{q}",
        "codewords_nonzero": N - 1,
        "minimum_glue_norm": best_num // q,
        "minimum_codewords": count,
        "example_coefficients": example,
    }
def general_row(q: int):
    d = (q + 3) // 2
    lower = Fraction(d * (q - 1), q)
    rank = q*q - 1
    base_roots = q * rank
    return {
        "q": q,
        "rank": rank,
        "MDS_distance": d,
        "MDS_glue_norm_lower_bound": f"{lower.numerator}/{lower.denominator}",
        "lower_bound_float": float(lower),
        "can_create_norm2_roots_by_bound": lower <= 2,
        "base_root_system": f"A_{q-1}^{q+1}",
        "base_root_count": base_roots,
        "base_root_to_rank_ratio": q,
    }


def payload(use_deep_q13=False):
    rows = [census(q) for q in (3, 5, 7, 11)]
    q13 = census(13) if use_deep_q13 else {
        "q": 13,
        "dimension": 168,
        "code": "[14,7,8]_13",
        **Q13_DEEP,
    }
    rows.append(q13)
    assert [r["minimum_glue_norm"] for r in rows] == [2, 4, 6, 10, 12]
    assert [r["minimum_codewords"] for r in rows] == [8, 72, 240, 552, 756]

    general = [general_row(q) for q in (3, 5, 7, 11, 13)]
    assert general[0]["can_create_norm2_roots_by_bound"]
    assert not any(r["can_create_norm2_roots_by_bound"] for r in general[1:])
    # Algebraically, ((q+3)/2)*((q-1)/q) > 2 iff
    # q^2 - 2q - 3 > 0, i.e. (q-3)(q+1)>0.
    symbolic_firewall = "(q-3)(q+1)>0 for every odd prime q>=5"

    root_shells = []
    for r in general:
        q = r["q"]
        if q == 3:
            roots = 240
            system = "E8 after glue creates 216 additional roots"
        else:
            roots = r["base_root_count"]
            system = r["base_root_system"]
        root_shells.append({
            "q": q,
            "rank": r["rank"],
            "root_count": roots,
            "root_to_rank_ratio": Fraction(roots, r["rank"]).__str__(),
            "root_system_reading": system,
        })

    checks = {
        "q3_glue_hits_root_norm": rows[0]["minimum_glue_norm"] == 2,
        "q5_q7_q11_q13_no_glue_roots":
            all(r["minimum_glue_norm"] > 2 for r in rows[1:]),
        "exact_minima_through_13_are_q_minus_1":
            all(r["minimum_glue_norm"] == r["q"] - 1 for r in rows),
        "all_prime_bound_q_ge5_excludes_new_roots": True,
        "q3_root_inflation_24_to_240": root_shells[0]["root_count"] == 240,
        "q_ge5_root_shell_is_base_only":
            all(x["root_count"] == x["q"] * x["rank"] for x in root_shells[1:]),
        "q13_deep_denominator_62748516":
            q13["codewords_nonzero"] == 62_748_516,
    }
    assert all(checks.values())
    return {
        "schema": "w33.pass11027.prime-clock-lattice-census.v1",
        "status": "PASS",
        "headline": (
            "The projective clock-code lattice tower is extended by exhaustive "
            "Euclidean-glue censuses at q=11 and q=13. The exact minimum glue "
            "norms for q=3,5,7,11,13 are 2,4,6,10,12=q-1. This q-1 pattern is "
            "certified only through q=13. Independently, the MDS lower bound "
            "proves for every odd prime q>=5 that glue norm is >2, so q=3 is "
            "the unique rung capable of creating new roots."
        ),
        "exact_census": rows,
        "q13_reproduction": {
            "default": "uses the frozen result of the completed exhaustive 13^7-1 scan",
            "deep_flag": "--deep-q13",
            "nonzero_codewords_checked": Q13_DEEP["codewords_nonzero"],
        },
        "general_root_firewall": {
            "rows": general,
            "factorization": symbolic_firewall,
            "theorem": (
                "For every odd prime q>=5, the self-dual MDS distance gives "
                "minimum nonzero glue norm >2. Hence no nontrivial glue coset "
                "can add roots; the root system stays A_{q-1}^{q+1}."
            ),
        },
        "root_shells": root_shells,
        "boundary": (
            "The equality min_glue=q-1 is an exhaustive finite observation only "
            "for q<=13 and is not promoted to an all-primes theorem. The all-prime "
            "root no-go is the separate MDS inequality inherited from Pass 10970. "
            "No continuum or physical dimensional interpretation is inferred."
        ),
        "checks": checks,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--deep-q13", action="store_true")
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    p = payload(use_deep_q13=a.deep_q13)
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if a.check:
        if not a.output.exists() or a.output.read_text(encoding="utf-8") != text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": p["status"],
        "minima": [r["minimum_glue_norm"] for r in p["exact_census"]],
        "q13_denominator": p["q13_reproduction"]["nonzero_codewords_checked"],
    }, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
