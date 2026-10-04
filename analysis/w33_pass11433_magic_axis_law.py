"""Pass 11433: the magic-axis law F' of Pass 11420 -- the grading structure behind it (proved), and an exhaustive
computer verification of its fixed-axis case at n = 4.

F' (Pass 11420): if M^2 z1 = z1, then v = z1 + M z1 is fixed by M, and every frame a with omega(a, v) != 0 violates.

GRADING LEMMA (proved).  Let M^2 z1 = z1 and v = z1 + M z1 != 0.  Then Pi = W(v) commutes with V_M (Mv = v) and with T1
(omega(z1, M z1) = 0 forces (M z1)_x1 = 0, so W(v) carries no X on qutrit 1).  Hence for U = W(a) V_M T1:
        Pi U Pi^-1 = omega^{c} U,   c = omega(v, a) (up to the sign convention),
i.e. U maps the Pi-eigenspace E_j to E_{j+c}.  Consequences, for c != 0 (exactly the frames F' calls violating):
  (G1) spec(U) is invariant under multiplication by omega (U is a cyclic block operator);
  (G2) tr(Pi^m U^r) = 0 unless 3 | r;
  (G3) any reversal Theta (antiunitary Clifford with Theta U Theta^-1 ~ U^-1) maps Pi to a Pauli that grades U^-1 with
       the same degree: Theta Pi Theta^-1 = Pi^T-type element times a Pauli commuting with U.
The grading alone does NOT force violation (the c != 0 condition is necessary in F', not yet sufficient by proof).
F' remains a conjecture; its proof needs the cubic phase of T1 across the Pi-blocks (L2 of Pass 11420 shows the phases
cancel around each cycle, so the obstruction is not a scalar holonomy).

EXHAUSTIVE AT n = 4 (fixed axis).  The fixed-axis cell IS H = Stab(z1) in Sp(8,3) (|H| = 20,056,328,248,320), and its
H-conjugation orbits are the conjugacy classes of H (GAP).  Every class representative is decided by the sparse decider
of Pass 11421; F' asks bad frames to contain {a : a_x1 != 0}.  Classes beyond the decider's cap are reported.
"""

from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11350_linear_decider as L  # noqa: E402
import w33_pass11352_three_qutrit_geometry as GEO  # noqa: E402
import w33_pass11421_four_qutrits as F4  # noqa: E402

OUT = ROOT / "data" / "w33_pass11433_magic_axis_law.json"
CLASSES = ROOT / "data" / "w33_pass11433_stab_z1_classes_n4.txt"
_S = {}


def grading_checks(n=2):
    """the grading lemma and (G1) on every n = 2 class of both 'bad' cells, all frames"""
    import w33_pass11330_orbit_census as O
    D = L.Decider(n)
    Ms, _ = O.all_symplectic(D.wl)
    z, Om, W, lab = D.z1, D.wl.Om, D.wl.W, D.wl.labels.astype(np.int64)
    w = np.exp(2j * np.pi / 3)
    st = Counter()
    for M in Ms:
        cell = GEO.cell(M, z, Om)
        if cell not in ("Mz1 = z1", "same line, M^2 z1 = z1"):
            continue
        v = (z + M @ z) % 3
        Pi = W[D.wl.index(v)]
        G = D.weil(M) @ D.T1
        st["Pi commutes with V_M T1"] += bool(np.allclose(Pi @ G, G @ Pi, atol=1e-9))
        for a in range(len(lab)):
            c = int(lab[a] @ Om @ v) % 3
            if c == 0:
                continue
            U = W[a] @ G
            X = Pi @ U @ Pi.conj().T
            graded = any(np.allclose(X, w ** k * U, atol=1e-9) for k in (1, 2))
            ev = np.linalg.eigvals(U)
            spec_inv = all(np.min(np.abs(ev - w * e)) < 1e-8 for e in ev)
            st["frames with c != 0"] += 1
            st["U graded by Pi"] += graded
            st["spec(U) omega-invariant"] += spec_inv
        st["classes"] += 1
    return dict(st)


def _init(cap=3 ** 11):
    L.CAP = cap
    _S["D"] = F4.SparseDecider(4)


def _job(args):
    i, line = args
    m, size = line.rsplit(";", 1)
    M = np.array(ast.literal_eval(m), dtype=np.int64).T % 3
    D = _S["D"]
    assert ((M @ D.z1 - D.z1) % 3 == 0).all()
    g = D.good_frames(M)
    if g is None:
        return dict(i=i, size=int(size), undecided=True)
    bad = ~g
    xne = D.wl.labels[:, 0] != 0
    return dict(i=i, size=int(size), verdict="bad == {a_x1 != 0}" if (bad == xne).all() else
                ("bad contains it" if (bad >= xne).all() else "FAILS"), nbad=int(bad.sum()))


def exhaustive_n4(nproc=6, second_cap=None):
    lines = CLASSES.read_text().splitlines()
    with Pool(nproc, initializer=_init) as pool:
        rows = pool.map(_job, list(enumerate(lines)), chunksize=4)
    und = [r["i"] for r in rows if r.get("undecided")]
    if second_cap:                       # optional second pass with a larger cap (3^14 ran > 2.5 CPU-hours per worker
        with Pool(nproc, initializer=_init, initargs=(second_cap,)) as pool:          # without finishing; off by default)
            again = {r["i"]: r for r in pool.map(_job, [(i, lines[i]) for i in und], chunksize=1)}
        rows = [again.get(r["i"], r) for r in rows]
    total = sum(r["size"] for r in rows)
    st = Counter()
    mass = Counter()
    for r in rows:
        key = "undecided" if r.get("undecided") else r["verdict"]
        st[key] += 1
        mass[key] += r["size"]
    return dict(classes=len(rows), group_order=total, first_pass_undecided=len(und), second_cap=second_cap,
                verdict_classes=dict(st),
                verdict_mass_fraction={k: v / total for k, v in mass.items()})


def run():
    res = dict(pass_id=11433)
    res["grading_n2"] = grading_checks()
    print(res["grading_n2"], flush=True)
    if CLASSES.exists():
        res["fixed_axis_n4_exhaustive"] = exhaustive_n4()
        print(res["fixed_axis_n4_exhaustive"], flush=True)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
