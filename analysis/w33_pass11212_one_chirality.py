#!/usr/bin/env python3
"""Pass 11212: one chirality -- the Gauss-sum twist of Pass 11201 and the Maslov chirality of Pass 11209 agree, and the
label-free form of the Choi-state invariant separates both time-reversal pairs.

Two tau-odd ("chiral") invariants of three-qutrit split relations were found in parallel:
  * Pass 11201 (master): I(sigma) = tr[rho^{(x)3} prod_p P_{sigma_p}], a third-moment invariant of the gate's Choi
    state and an exact quadratic Gauss sum.  For one labelled pattern CHIRAL it is Q i sqrt3/243 on the 256 pair, with Q
    the twist of Pass 11189.  Its stated caveat: the sign needs registers labelled compatibly with the relation.
  * Pass 11209 (branch): the profile-oriented Kashiwara-Maslov count chi = #{class 1} - #{class 3}, -+54 on the 256
    pair and +-18 on the 6912 pair.  It sums over all product Lagrangians of both splits, so it needs no labelling.

This pass compares them on the same relations and closes the labelling caveat.
 1. ONE CHIRALITY ON THE 256 PAIR.  On all 512 members of the two 256 orbitals: chi = 54 Q and
    Im I(CHIRAL) = Q sqrt3/243, so chi = 54 * 243/sqrt3 * Im I(CHIRAL) exactly.
 2. THE 6912 PATTERN HAS THE RIGHT SIGN BUT IS NOT AN INVARIANT.  Pass 11201's pattern for the 6912 pair was evaluated
    on one representative per orbital.  Across members of the same orbital (which carry different labellings of F's
    factors) its value is 0 on most members of the chiral pair and -sign(chi) on the rest, while on tau-FIXED
    (achiral) 6912 orbitals it takes both signs.  The sign is the relation's; the value is the labelling's.
 3. THE LABEL-FREE CHOI INVARIANT.  Tag each relabelling (pi_in, pi_out) of the parties by the block-rank code matrix
    of the relation read in the pattern's slots (codes: 0 zero, 1 rank one, 2/3 invertible with det 1/2; local-Clifford
    invariant and tau-blind).  The multiset over all 36 relabellings of (tag, I) does not depend on how the relation's
    factors are labelled; it is constant on orbitals and conjugated by tau.  It is chiral on BOTH reversal pairs and
    achiral on every tau-fixed orbital, for CHIRAL and for the 6912 pattern.  The Choi-state phase therefore measures
    the arrow on unlabelled registers too, once the relabellings are sorted by block type -- the same device (a
    gauge-invariant orientation) that makes the Maslov count label-free.
"""
from __future__ import annotations

import itertools
import json
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11182_paper_ticks_mereology as P  # noqa: E402
import w33_pass11184_merged_orbitals as O  # noqa: E402
import w33_pass11189_oriented_holonomy as H  # noqa: E402
import w33_pass11201_twist_chiral_invariant as G  # noqa: E402
import w33_pass11209_maslov_chirality as M  # noqa: E402

OUT = ROOT / "data" / "w33_pass11212_one_chirality.json"
CACHE = Path(tempfile.gettempdir()) / "w33_pass11212_orbit_labels.npy"
TAU = M.TAU
SQ3 = np.sqrt(3)
PAT_6912 = [tuple(int(c) for c in s) for s in json.loads(
    (ROOT / "data" / "w33_pass11201_twist_chiral_invariant.json").read_text())["pair_6912"]["pattern"]]
UNITS = {"CHIRAL": 243 / SQ3, "PAT_6912": 729 / SQ3}
PATTERNS = {"CHIRAL": G.CHIRAL, "PAT_6912": PAT_6912}


def orbit_data(Bs):
    """orbital label of every split, plus reverse and tau-image orbitals"""
    orbs, index = O.orbits(Bs, return_index=True)
    orb_of = np.zeros(len(Bs), np.int64)
    for k, o in enumerate(orbs):
        orb_of[o] = k
    rows = []
    for k, o in enumerate(orbs):
        rev = int(O.reverse_orbit(Bs[o[0]], index, orb_of))
        tau = int(orb_of[index[O.fac_key((TAU @ Bs[o[0]]) % 3)]])
        rows.append(dict(orbital=k, size=len(o), reverse=rev, tau_image=tau))
    return orb_of, rows


def load_orbits(Bs):
    if CACHE.exists():
        d = np.load(CACHE, allow_pickle=True).item()
        return d["orb_of"], d["rows"]
    orb_of, rows = orbit_data(Bs)
    np.save(CACHE, dict(orb_of=orb_of, rows=rows), allow_pickle=True)
    return orb_of, rows


def value(B, pattern, unit):
    v = G.gauss_inv(G.lagrangian(B % 3), pattern)
    return round(v.imag * unit, 6), round(v.real * 3 ** 6, 6)


def tagged_multiset(B, pattern, unit):
    """Counter over the 36 relabellings of (block codes read in the pattern's slots, I in units).  Parties 0-2 of the
    Choi state are the inputs (columns of T = the factors of F), parties 3-5 the outputs (rows of T = the factors of F0),
    so the block between input slot a and output slot b is code[row of b, column of a]."""
    code = H.codes(H.blocks(B))
    L = G.lagrangian(B % 3)
    out = Counter()
    for pi in itertools.permutations(range(3)):
        for po in itertools.permutations(range(3)):
            s = [pattern[pi[q]] for q in range(3)] + [pattern[3 + po[q]] for q in range(3)]
            v = G.gauss_inv(L, s)
            pinv, poinv = np.argsort(pi), np.argsort(po)
            tag = tuple(int(code[poinv[b], pinv[a]]) for a in range(3) for b in range(3))
            out[(tag, round(v.real * 3 ** 6, 6), round(v.imag * unit, 6))] += 1
    return out


def conj(ms):
    return Counter({(t, re, -im if im else 0.0): c for (t, re, im), c in ms.items()})


def representatives(Bs, orb_of, rows):
    """one adapted basis per tau-swapped orbital and per 6912 orbital (lets the regression recompute without the
    110565-split enumeration)"""
    return {str(r["orbital"]): dict(size=r["size"], tau_image=r["tau_image"],
                                    B=Bs[int(np.flatnonzero(orb_of == r["orbital"])[0])].tolist())
            for r in rows if r["tau_image"] != r["orbital"] or r["size"] == 6912}


def run(members_chi=None, members_label=64, members_tagged=4, seed=11212):
    t0 = time.time()
    Bs = P.load_facs()
    orb_of, rows = load_orbits(Bs)
    by = {r["orbital"]: r for r in rows}
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11212, orbitals=len(rows), orbit_sizes=sorted(r["size"] for r in rows))
    chiral_pairs = sorted({tuple(sorted((r["orbital"], r["tau_image"]))) for r in rows
                           if r["tau_image"] != r["orbital"]})
    res["tau_swapped_pairs"] = [[by[a]["size"], a, b] for a, b in chiral_pairs]
    res["tau_swapped_equals_reversal"] = all(by[a]["reverse"] == b for a, b in chiral_pairs)

    # 1. the 256 pair, every member
    pair256 = next(p for p in chiral_pairs if by[p[0]]["size"] == 256)
    table = Counter()
    for o in pair256:
        mem = np.flatnonzero(orb_of == o)
        if members_chi is not None:
            mem = mem[:members_chi]
        for i in mem:
            B = Bs[i]
            chi = M.chirality(M.spectrum(B))
            q = H.twist(B)
            q = None if q is None else (1 if q == 1 else -1)
            im, re = value(B, G.CHIRAL, UNITS["CHIRAL"])
            table[(o, chi, q, im, re)] += 1
    res["pair_256"] = [dict(orbital=o, chi=chi, Q=q, imag_units_i_sqrt3_over_243=im, real_times_729=re, members=c)
                       for (o, chi, q, im, re), c in sorted(table.items())]
    res["chi_equals_54Q"] = all(r["chi"] == 54 * r["Q"] for r in res["pair_256"])
    res["imag_equals_Q"] = all(r["imag_units_i_sqrt3_over_243"] == r["Q"] for r in res["pair_256"])
    res["pair_256_members_checked"] = sum(table.values())

    # 2. Pass 11201's 6912 pattern, member by member, on every 6912 orbital
    lab = {}
    for r in rows:
        if r["size"] != 6912:
            continue
        mem = np.flatnonzero(orb_of == r["orbital"])
        pick = rng.choice(mem, min(members_label, len(mem)), replace=False)
        c = Counter(value(Bs[i], PAT_6912, UNITS["PAT_6912"])[0] for i in pick)
        chis = Counter(M.chirality(M.spectrum(Bs[i])) for i in pick[:4])
        lab[str(r["orbital"])] = dict(tau_fixed=r["tau_image"] == r["orbital"], chi=sorted(chis),
                                      imag_values_units_i_sqrt3_over_729={str(k): v for k, v in sorted(c.items())},
                                      members=len(pick))
    res["pattern_6912_by_member"] = lab
    res["pattern_6912_constant_on_orbitals"] = all(len(v["imag_values_units_i_sqrt3_over_729"]) == 1
                                                   for v in lab.values())
    res["pattern_6912_nonzero_on_achiral_orbital"] = any(
        v["tau_fixed"] and any(k != "0.0" for k in v["imag_values_units_i_sqrt3_over_729"]) for v in lab.values())

    # 3. the label-free (tagged) Choi multiset on every orbital
    tagged = []
    for r in rows:
        mem = np.flatnonzero(orb_of == r["orbital"])
        pick = [mem[0]] + list(rng.choice(mem, min(members_tagged, len(mem)) - 1, replace=False)) \
            if len(mem) > 1 else [mem[0]]
        row = dict(orbital=r["orbital"], size=r["size"], tau_fixed=r["tau_image"] == r["orbital"])
        for name, pat in PATTERNS.items():
            ms = [tagged_multiset(Bs[i], pat, UNITS[name]) for i in pick]
            row[name] = dict(constant=all(m == ms[0] for m in ms), chiral=ms[0] != conj(ms[0]),
                             imaginary_entries=sum(c for (t, re, im), c in ms[0].items() if im),
                             key=sorted([list(t), re, im, c] for (t, re, im), c in ms[0].items()))
        row["members"] = len(pick)
        tagged.append(row)
    keys = {r["orbital"]: {n: json.dumps(r[n]["key"]) for n in PATTERNS} for r in tagged}
    res["tagged"] = [{k: v for k, v in r.items() if k not in PATTERNS} |
                     {n: {kk: vv for kk, vv in r[n].items() if kk != "key"} for n in PATTERNS} for r in tagged]
    res["tagged_constant_on_orbitals"] = all(r[n]["constant"] for r in tagged for n in PATTERNS)
    res["tagged_chiral_iff_tau_swapped"] = {n: all(r[n]["chiral"] == (not r["tau_fixed"]) for r in tagged)
                                            for n in PATTERNS}
    res["tagged_chiral_on"] = {n: sorted(r["size"] for r in tagged if r[n]["chiral"]) for n in PATTERNS}
    ta = {n: {} for n in PATTERNS}
    for a, b in chiral_pairs:
        for n in PATTERNS:
            ma = next(r for r in tagged if r["orbital"] == a)[n]["key"]
            mb = next(r for r in tagged if r["orbital"] == b)[n]["key"]
            ta[n][f"{by[a]['size']}:{a}-{b}"] = dict(separated=ma != mb)
    res["tagged_separates_pairs"] = ta
    res["tagged_classifies"] = {n: len({keys[o][n] for o in keys}) for n in PATTERNS}
    res["tagged_both_classify"] = len({(keys[o]["CHIRAL"], keys[o]["PAT_6912"]) for o in keys})
    res["representatives"] = representatives(Bs, orb_of, rows)
    res["seconds"] = round(time.time() - t0, 1)
    return res


def main():
    res = run()
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: v for k, v in res.items() if k not in ("tagged", "pair_256", "representatives")}, indent=1))
    for r in res["pair_256"]:
        print(r)


if __name__ == "__main__":
    main()
