"""Pass 11373: the exact three-qutrit one-magic-gate fraction, from orbits of the magic-axis stabiliser.

Setting (Passes 11309/11330/11352).  U = W(a) V_M (T (x) I (x) I), M in Sp(6,3), a in F3^6 (729 frames).  Pass 11352
sampled the violating fraction 0.1333 +- 0.0019 and excluded the two-qutrit value 1/8 at 4.4 sigma.

Reduction (exact).  Let H = Stab(z1) in Sp(2n,3) (z1 = the magic axis Z_1).  Every N in H lifts to a Clifford D with
D Z_1 D^dag = Z_1 exactly (a Pauli fixes the frame), hence D T_1 D^dag = T_1, and
    D (W(a) V_M T_1) D^dag = W(a') V_{N M N^-1} T_1        (a -> a' a bijection of frames).
Reversibility is conjugation-invariant, so the number of violating frames is constant on H-conjugation orbits of Sp.
GAP enumerates the orbits as C_G(g) \\ G / H double cosets over conjugacy-class representatives g (row convention:
GAP's matrices are the transposes of ours, and the column stabiliser corresponds to the row stabiliser OnRight).
Each orbit representative is decided by the F3-linear decider of Pass 11350; classes whose (S) solution space
exceeds its cap (M = +-I and a few others) by the exact Weyl criterion (Pass 11252) on all frames.

Slow classes (the (S) solution space has affine dimension > 11, all of high symmetry): for N in C_H(M), the
centraliser of M in H, V_N (W(a) V_M T_1) V_N^dag = W(N a) V_M T_1 up to phase (the Weil representation is projective
and V_N T_1 V_N^dag = T_1), so the violating frames are a union of orbits of the LINEAR action a -> N a of C_H(M).
GAP supplies generators of C_H(M) (checked: |H : C_H(M)| = orbit size); one frame per orbit is decided by the exact
Weyl criterion.  Control: for every slow class two random non-representative frames are decided directly and must
agree with their orbit's verdict.

Validation: the same pipeline at n = 2 (200 orbits) reproduces the exhaustive 223/2430 and 6480 bad classes of
Passes 11309/11330, and the n = 3 orbit sizes sum to |Sp(6,3)| = 9,170,703,360.
"""

from __future__ import annotations

import ast
import json
import subprocess
import sys
from collections import defaultdict
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11252_exact_reversibility as R  # noqa: E402
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402

OUT = ROOT / "data" / "w33_pass11373_three_qutrit_exact_fraction.json"
ORB = {n: ROOT / "data" / f"w33_pass11373_orbits_n{n}.txt" for n in (2, 3)}
CEN = {n: ROOT / "data" / f"w33_pass11373_centralisers_n{n}.txt" for n in (2, 3)}
SP_ORDER = {2: 51840, 3: 9170703360}

GAP_SCRIPT = r"""
Orbs := function(gens, n, fname)
  local G, z1, H, cc, g, C, dc, d, y, out, tot;
  G := Group(gens);
  z1 := List([1..2*n], i -> 0*Z(3)); z1[2] := Z(3)^0;
  H := Stabilizer(G, z1, OnRight);
  cc := ConjugacyClasses(G);
  out := OutputTextFile(fname, false);
  SetPrintFormattingStatus(out, false);
  tot := 0;
  for g in List(cc, Representative) do
    C := Centralizer(G, g);
    dc := DoubleCosetRepsAndSizes(G, C, H);
    for d in dc do
      y := d[1]^-1 * g * d[1];
      AppendTo(out, List(y, r -> List(r, IntFFE)), ";", d[2] / Size(C), "\n");
      tot := tot + d[2] / Size(C);
    od;
  od;
  CloseStream(out);
  Print("total ", tot, " ", tot = Size(G), "\n");
end;
"""


GAP_CENTRALISERS = r"""
G := Group(gens);
z1 := List([1..2*n], i -> 0*Z(3)); z1[2] := Z(3)^0;
H := Stabilizer(G, z1, OnRight);
out := OutputTextFile(fname, false);
SetPrintFormattingStatus(out, false);
for j in [1..Length(reps)] do
  C := Centralizer(H, reps[j]);
  if Size(H) / Size(C) <> sizes[j] then Error("size mismatch"); fi;
  AppendTo(out, idx[j], ";", Size(C), ";", List(GeneratorsOfGroup(C), g -> List(g, r -> List(r, IntFFE))), "\n");
od;
CloseStream(out);
"""


def gap_generators(n):
    gens = GEO.gen_mats(n)
    return "gens := [\n" + ",\n".join(
        "[" + ",".join("[" + ",".join(str(int(x)) for x in row) + "]" for row in g.T) + "]*Z(3)^0" for g in gens) + "];\n"


def wsl_path(p):
    p = str(p)
    return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")


def build_orbits(n):
    """run GAP in WSL (hours at n = 3); the committed orbit list is reused when present"""
    if ORB[n].exists():
        return
    work = ROOT / "data" / f"_gap_11373_n{n}.g"
    target = wsl_path(ORB[n])
    work.write_text(GAP_SCRIPT + gap_generators(n) + f'Orbs(gens, {n}, "{target}");\nQUIT;\n')
    wpath = wsl_path(work)
    subprocess.run(["wsl", "-e", "bash", "-lc", f'gap -q -o 8g "{wpath}"'], check=True)
    work.unlink()


def build_centralisers(n, slow_idx):
    """GAP: generators of C_H(M) for the slow orbit representatives (reused when the file exists)"""
    if CEN[n].exists():
        return
    lines = ORB[n].read_text().splitlines()
    work = ROOT / "data" / f"_gap_11373_cen_n{n}.g"
    reps = ",\n".join(lines[i].rsplit(";", 1)[0] + "*Z(3)^0" for i in slow_idx)
    body = (gap_generators(n) + f"n := {n};\n" + f"reps := [{reps}];\n"
            + "sizes := [" + ",".join(lines[i].rsplit(";", 1)[1] for i in slow_idx) + "];\n"
            + "idx := [" + ",".join(str(i) for i in slow_idx) + "];\n"
            + f'fname := "{wsl_path(CEN[n])}";\n' + GAP_CENTRALISERS + "QUIT;\n")
    work.write_text(body)
    subprocess.run(["wsl", "-e", "bash", "-lc", f'gap -q -o 8g "{wsl_path(work)}" < /dev/null'], check=True, stdin=subprocess.DEVNULL)
    work.unlink()


def read_orbits(n):
    rows = []
    for line in ORB[n].read_text().splitlines():
        m, s = line.strip().rsplit(";", 1)
        rows.append((np.array(ast.literal_eval(m), dtype=np.int64).T % 3, int(s)))
    return rows


_S = {}


def read_centralisers(n):
    """index -> (|C_H(M)|, generators as COLUMN-convention matrices N = K^T)"""
    out = {}
    for line in CEN[n].read_text().splitlines():
        i, size, gens = line.split(";", 2)
        out[int(i)] = (int(size), [np.array(K, dtype=np.int64).T % 3 for K in ast.literal_eval(gens)])
    return out


def slow(D, M):
    for k in range(3):
        sol = D.symplectic_solutions(M, k)
        if sol is not None and 3 ** len(sol[1]) > L.CAP:
            return True
    return False


def frame_orbits(D, gens):
    labels = D.wl.labels.astype(np.int64)
    pow3 = 3 ** np.arange(labels.shape[1])[::-1]
    idx = {int(l @ pow3): j for j, l in enumerate(labels)}
    imgs = [np.array([idx[int(((N @ labels.T) % 3)[:, j] @ pow3)] for j in range(len(labels))]) for N in gens]
    seen = -np.ones(len(labels), dtype=int)
    orbits = []
    for j in range(len(labels)):
        if seen[j] >= 0:
            continue
        orb, stack = [j], [j]
        seen[j] = len(orbits)
        while stack:
            x = stack.pop()
            for im in imgs:
                y = int(im[x])
                if seen[y] < 0:
                    seen[y] = len(orbits)
                    orb.append(y)
                    stack.append(y)
        orbits.append(orb)
    return orbits, imgs


def _init(n):
    L.CAP = 3 ** 11
    _S["D"] = L.Decider(n)


def _fast(args):
    i, M = args
    D = _S["D"]
    g = D.good_frames(M)
    return i, int((~g).sum()), "linear", GEO.cell(M, D.z1, D.wl.Om)


def _slow_orbit(args):
    i, M, a = args
    D = _S["D"]
    U = D.wl.W[a] @ D.weil(M) @ D.T1
    v = R.decide(U, D.n, np.random.default_rng(a))[0]
    assert v is not None, "Weyl enumeration limit hit"
    return i, a, v


def census(n, nproc=8):
    rows = read_orbits(n)
    total = sum(s for _, s in rows)
    assert total == SP_ORDER[n], (n, total)
    Dm = L.Decider(n)
    L.CAP = 3 ** 11
    slow_idx = [i for i, (M, _) in enumerate(rows) if slow(Dm, M)]
    build_centralisers(n, slow_idx)
    cen = read_centralisers(n)
    assert set(slow_idx) <= set(cen), "missing centraliser generators"
    H_order = {2: 648, 3: 12597120}[n]
    assert all(H_order // cen[i][0] == rows[i][1] for i in slow_idx)
    nfr = 3 ** (2 * n)
    with Pool(nproc, initializer=_init, initargs=(n,)) as pool:
        dec = {i: (b, how, cell) for i, b, how, cell in
               pool.imap_unordered(_fast, [(i, M) for i, (M, _) in enumerate(rows) if i not in set(slow_idx)],
                                   chunksize=8)}
        jobs, orbit_of, extra = [], {}, []
        rng = np.random.default_rng(11373)
        for i in slow_idx:
            orbits, _ = frame_orbits(Dm, cen[i][1])
            orbit_of[i] = orbits
            jobs += [(i, rows[i][0], orb[0]) for orb in orbits]
            # control: two random frames that are NOT orbit representatives, decided directly
            others = [(a, o) for o, orb in enumerate(orbits) for a in orb[1:]]
            for t in rng.choice(len(others), size=min(2, len(others)), replace=False) if others else []:
                extra.append((i, others[t][1], others[t][0]))
                jobs.append((i, rows[i][0], others[t][0]))
        verdict = {}
        for i, a, v in pool.imap_unordered(_slow_orbit, jobs, chunksize=1):
            verdict[(i, a)] = v
    slow_stats = []
    for i in slow_idx:
        b = sum(len(orb) for orb in orbit_of[i] if verdict[(i, orb[0])] is False)
        dec[i] = (b, "weyl-orbits", GEO.cell(rows[i][0], Dm.z1, Dm.wl.Om))
        slow_stats.append(dict(i=i, size=rows[i][1], frame_orbits=len(orbit_of[i]), bad_frames=b))
    # control: directly decided non-representative frames agree with their orbit's verdict
    ctrl = 0
    for i, o, a in extra:
        assert verdict[(i, a)] == verdict[(i, orbit_of[i][o][0])], f"orbit invariance fails at class {i}, frame {a}"
        ctrl += 1
    bad_frames = sum(rows[i][1] * dec[i][0] for i in dec)
    bad_classes = sum(rows[i][1] for i in dec if dec[i][0])
    cells = defaultdict(lambda: defaultdict(int))
    frame_hist = defaultdict(int)
    for i, (b, _, cell) in dec.items():
        s = rows[i][1]
        cells[cell]["classes"] += s
        cells[cell]["bad_classes"] += s * (b > 0)
        cells[cell]["bad_frames"] += s * b
        frame_hist[b] += s
    P = Fraction(bad_frames, total * nfr)
    return dict(n=n, orbits=len(rows), sum_of_orbit_sizes=total, slow_classes_by_frame_orbits=len(slow_idx),
                slow_detail=slow_stats, invariance_controls=ctrl,
                violating_fraction=str(P), violating_fraction_float=float(P),
                bad_class_fraction=str(Fraction(bad_classes, total)),
                bad_frames_per_class_histogram={str(k): v for k, v in sorted(frame_hist.items())},
                cells={c: dict(v, bad_class_fraction=str(Fraction(v["bad_classes"], v["classes"])),
                               violating_fraction=str(Fraction(v["bad_frames"], v["classes"] * nfr)))
                       for c, v in sorted(cells.items())})


def reversal_route(n=2, tried=40):
    """PROOF ROUTE for the universal rule 'Mz1 = -z1 => reversible in every frame'.  At k = 0, (S) reads
    Q^-1 M Q = J M^-1 J with Q z1 = z1; any anti-symplectic involution sigma with sigma M sigma = M^-1 (Wonenburger: such
    sigma exist for every symplectic M) and sigma z1 = -z1 gives the solution Q = sigma J.  Counted at n = 2 over the
    whole reversed cell: existence of such sigma, and how many frames a single sigma already makes reversible."""
    import w33_pass11330_orbit_census as O
    D = L.Decider(n)
    Ms, _ = O.all_symplectic(D.wl)
    Ms = np.array(Ms) % 3
    J, z = D.J, D.z1
    Sig = (J[None] @ Ms) % 3
    ok = ((np.einsum('aij,ajk->aik', Sig, Sig) % 3) == np.eye(2 * n, dtype=np.int64)).all(axis=(1, 2))
    Sig = Sig[ok & (((Sig @ z) % 3) == ((-z) % 3)).all(axis=1)]
    labels = D.wl.labels.astype(np.int64)
    N2 = 2 * n
    e1 = np.eye(N2, dtype=np.int64)[0]

    def frames_from_Q(M, Q):
        Minv = R._inv_mod3(M)
        rQ = D.r_frame(Q)
        Lmat = (Minv + Q @ J) % 3
        rhs0 = (-(D.gframe(e1 * 0) % 3) - rQ) % 3
        Y = ((np.eye(N2, dtype=np.int64) - Minv) % 3)[:, 1:]
        _, Nb = L.solve_affine(Y.T, np.zeros(N2 - 1, dtype=np.int64))
        if not Nb:
            return np.ones(len(labels), bool)
        vals = (np.array(Nb) @ ((rhs0[None, :] - labels @ Lmat.T) % 3).T) % 3
        return (vals == 0).all(axis=0)

    stats = defaultdict(int)
    for M in Ms:
        if not ((M @ z + z) % 3 == 0).all():
            continue
        Minv = R._inv_mod3(M)
        rev = [sg for sg in Sig if ((sg @ M @ sg - Minv) % 3 == 0).all()]
        best = 0
        for sg in rev[:tried]:
            Q = (sg @ J) % 3
            assert ((M @ Q @ (J @ M @ J) - Q) % 3 == 0).all()
            best = max(best, int(frames_from_Q(M, Q).sum()))
        stats[f"sigma exists: {bool(rev)}, best single-sigma frames: {best}"] += 1
    return dict(n=n, involutions_with_sigma_z1_eq_minus_z1=int(len(Sig)), classes=dict(stats))


def run():
    res = dict(pass_id=11373)
    for n in (2, 3):
        build_orbits(n)
        res[f"n={n}"] = census(n)
        print(json.dumps(res[f"n={n}"], indent=1), flush=True)
    res["n2_matches_223_2430"] = res["n=2"]["violating_fraction"] == "223/2430"
    res["n2_bad_classes_is_one_eighth"] = res["n=2"]["bad_class_fraction"] == "1/8"
    b3 = float(Fraction(res["n=3"]["bad_class_fraction"]))
    res["n3_bad_class_fraction_float"] = b3
    res["n3_vs_pass11352_sample_z"] = (b3 - 0.1333) / 0.0019         # Pass 11352 sampled the bad-CLASS fraction
    return res


def main():
    if "--route" in sys.argv:
        res = json.load(open(OUT))
        res["reversal_route_n2"] = reversal_route()
        print(res["reversal_route_n2"])
        json.dump(res, open(OUT, "w"), indent=1, default=str)
        return
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
