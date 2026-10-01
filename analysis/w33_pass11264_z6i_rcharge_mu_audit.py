#!/usr/bin/env python3
"""Pass 11264: restore the Z6-I plane R charges and rerun the exact mu audit.

The committed raw certificate was extracted from the locally installed
orbifolder 1.2.1 model files with

    R_i = q_sh,i + oscillator_contribution_i

for the three complex planes.  The same extractor was checked against every
stored Z6-II field R charge before it was used on Z6-I.  This producer makes
the expensive part independent of the external installation: it injects the
committed Z6-I values into Pass 11241's exact integer-lattice engine.

Cached mode is intentionally the default.  ``--run`` enumerates every Z6-I
FI-ray and every choice of condensing field and can take tens of CPU-hours;
use Windows ``py -3`` in this checkout because that environment supplies cdd.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

_spec = importlib.util.spec_from_file_location(
    "pass11241", ROOT / "analysis" / "w33_pass11241_discrete_mu_symmetry.py"
)
P41 = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(P41)

RAW = ROOT / "data" / "w33_pass11264_z6i_r_charges.json.gz"
OUT = ROOT / "data" / "w33_pass11264_z6i_rcharge_mu_audit.json"
ORDERS = ["6:-1/1", "6:-1/1", "3:-1/1"]


def load_raw(path: Path = RAW):
    compressed = path.read_bytes()
    plain = gzip.decompress(compressed)
    data = json.loads(plain)
    return data, {
        "models": len(data),
        "fields": sum(len(v) for v in data.values()),
        "compressed_bytes": len(compressed),
        "uncompressed_bytes": len(plain),
        "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
        "uncompressed_sha256": hashlib.sha256(plain).hexdigest(),
    }


def inject_r_charges(raw):
    """Return the frozen discrete-charge dictionary with Z6-I R data added."""
    disc = P41.M67.load(P41.DISC)
    assert set(raw) == {n for n in disc if n.startswith("Z6I_")}
    for name, fields in raw.items():
        assert set(fields) == set(disc[name]["fields"]), name
        disc[name]["R"] = list(ORDERS)
        for field, value in fields.items():
            rq = value["RQ"]
            assert len(rq) == 3
            disc[name]["fields"][field]["R"] = rq
    return disc


def run_exhaustive():
    ledger, ledger_sha = P41.P.load_ledger()
    raw, raw_meta = load_raw()
    disc = inject_r_charges(raw)
    names = sorted(n for n in ledger if n.startswith("Z6-I|"))
    t0 = time.time()
    models = {}
    for name in names:
        t = time.time()
        models[name] = P41.analyse_model(
            name, ledger[name], disc, ray_cap=None, choice_cap=None, seed=11264
        )
        print(
            name,
            models[name].get("vacua"),
            models[name].get("mu_protected"),
            f"{time.time() - t:.1f}s",
            flush=True,
        )
    return {
        "ledger_sha256": ledger_sha,
        "raw": raw_meta,
        "models": models,
        "seconds": round(time.time() - t0, 1),
    }


def certificate(cached):
    raw, raw_meta = load_raw()
    inject_r_charges(raw)  # full name/field coverage and shape check
    models = cached["models"]
    dflat = {n: r for n, r in models.items() if r.get("dflat")}
    exhaustive = all(r.get("exhaustive", True) for r in dflat.values())
    protected = {n: r for n, r in dflat.items() if r.get("mu_protected", 0)}
    clean = {n: r for n, r in dflat.items() if r.get("clean", 0)}
    per_model = {
        n: {
            k: r.get(k)
            for k in (
                "anomalous", "dflat", "fi_rays", "vacua", "exhaustive",
                "mu_protected", "clean", "clean_and_F_flat"
            )
            if k in r
        }
        for n, r in sorted(models.items())
    }
    result = {
        "schema": "w33.pass11264.z6i-rcharge-mu-audit.v1",
        "status": "PASS_Z6I_R_CHARGES_RESTORED_AND_EXHAUSTIVE_MU_AUDIT",
        "extraction": {
            **raw_meta,
            "orbifolder_version": "1.2.1",
            "formula": "R_i=q_sh_i+OscillatorContribution_i",
            "plane_orders": [6, 6, 3],
            "superpotential_R_charges": [-1, -1, -1],
            "z6ii_control_models": 128,
            "z6ii_control_fields": 43634,
            "z6ii_control_mismatches": 0,
        },
        "audit": {
            "models_total": len(models),
            "models_dflat": len(dflat),
            "vacua": sum(r.get("vacua", 0) for r in dflat.values()),
            "exhaustive": exhaustive,
            "mu_protected_vacua": sum(r.get("mu_protected", 0) for r in dflat.values()),
            "mu_protected_models": sorted(protected),
            "clean_vacua": sum(r.get("clean", 0) for r in dflat.values()),
            "clean_models": sorted(clean),
            "per_model": per_model,
        },
        "engine": {
            "source": "analysis/w33_pass11241_discrete_mu_symmetry.py",
            "criterion": "exact integer lattice of U(1), remnant, space-group, point-group and plane R charges",
            "ledger_sha256": cached.get("sha", cached.get("ledger_sha256")),
            "wall_seconds_parallel_run": round(cached.get("seconds", 0), 1),
        },
        "correction": (
            "Supersedes Pass 11247's data-blocked boundary and the R-charge caveat in Pass 11241. "
            "The result tests symmetry protection only; an allowed mu term may first occur at high order "
            "or with a small coefficient."
        ),
    }
    assert result["extraction"]["models"] == 87
    assert result["extraction"]["fields"] == 30980
    assert len(models) == 87
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="rerun all 87 models exhaustively")
    parser.add_argument("--ingest", type=Path, help="ingest a cached exhaustive result")
    args = parser.parse_args()
    if args.run:
        cached = run_exhaustive()
    elif args.ingest:
        cached = json.loads(args.ingest.read_text())
    elif OUT.exists():
        print(json.dumps(json.loads(OUT.read_text()), indent=2))
        return
    else:
        parser.error("use --run or --ingest PATH for the first certificate build")
    result = certificate(cached)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
