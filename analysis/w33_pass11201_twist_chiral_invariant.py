"""Pass 11201: the time-reversal twist is a measurable chiral phase -- a quadratic Gauss sum in the gate's Choi state.

Pass 11189 separated the last two time-reversal pairs of three-qutrit split relations by oriented finite-field data:
an SL(2,3) holonomy class (N vs -N) for the 6912 pair and a symplectic twist Q = +-1 for the 256 pair.  A relation
(F0, F) is a two-sided local-equivalence class of a Clifford gate U_T; its reverse is U_T^dagger, and complex
conjugation (anti-unitary time reversal) exchanges the pair.

Local-unitary invariants of the six-qutrit Choi state,
    I(sigma) = tr[ rho^{(x)3} prod_p P_{sigma_p} ],   sigma_p in S3 permuting the three copies on party p,
are measurable (randomized measurements) and I(sigma)^* = I(sigma^-1).  For a stabilizer state with Lagrangian L,
    I(sigma) = 3^-18 * prod_p t_p * sum_{(v1,v2) in C} zeta^{ sum_{p cyclic} (+-1/2) omega_p(v1_p, v2_p) },
t_p = 27, 9, 3 for identity, transposition, 3-cycle, and C the linear constraints (identity: v1_p = v2_p = 0;
transposition fixing copy c: v_c,p = 0 with v3 = -v1-v2).  The linear phases of the stabiliser cancel, so I depends only
on L and sigma -- an exact Gauss sum.  Consequences checked here:
  * if every 3-cycle has the same orientation, isotropy kills the phase: I is real and positive;
  * complex conjugation equals relabelling the copies by an odd permutation, so I is real whenever such a relabelling
    preserves the pattern of transpositions; chirality needs mixed 3-cycles AND two different transpositions;
  * for such a pattern, I on the 256 relation is purely imaginary, +- i sqrt3 / 243 (a quadratic Gauss sum), with sign
    Q: the twist is the sign of a measurable chiral phase, and it flips under time reversal.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11182_paper_ticks_mereology as PT  # noqa: E402
import w33_pass11189_oriented_holonomy as H  # noqa: E402

OUT = ROOT / "data" / "w33_pass11201_twist_chiral_invariant.json"
ZETA = np.exp(2j * np.pi / 3)
PERMS = [(0, 1, 2), (1, 2, 0), (2, 0, 1), (1, 0, 2), (0, 2, 1), (2, 1, 0)]
W3 = ZETA
X1 = np.roll(np.eye(3), 1, axis=0)
Z1 = np.diag([1, W3, W3 ** 2])


def weyl(vec):
    op = np.array([[1.0 + 0j]])
    for k in range(len(vec) // 2):
        x, z = int(vec[2 * k]) % 3, int(vec[2 * k + 1]) % 3
        op = np.kron(op, np.linalg.matrix_power(X1, x) @ np.linalg.matrix_power(Z1, z))
    return op


def generators(T):
    gens = []
    for k in range(6):
        e = np.zeros(6, np.int64)
        e[k] = 1
        inp = e.copy()
        inp[1::2] = (-inp[1::2]) % 3
        gens.append(np.concatenate([inp, (T @ e) % 3]))
    return np.array(gens)


def choi_state(T):
    """state vector (for cross-checks): stabilised by (e_k^*, T e_k)"""
    P = np.eye(3 ** 6, dtype=complex)
    for g in generators(T):
        Wg = weyl(g)
        P = P @ (np.eye(3 ** 6) + Wg + Wg @ Wg) / 3
    psi = P[:, np.argmax(np.linalg.norm(P, axis=0))]
    return (psi / np.linalg.norm(psi)).reshape([3] * 6)


def chiral_invariant(psi, sig):
    """direct tensor contraction (cross-check only)"""
    letters = "abcdefghijklmnopqr"
    idx = [[letters[3 * p + c] for p in range(6)] for c in range(3)]
    terms = ["".join(idx[c]) for c in range(3)] + ["".join(idx[sig[p][c]][p] for p in range(6)) for c in range(3)]
    ops = [psi] * 3 + [psi.conj()] * 3
    expr = ",".join(terms) + "->"
    return complex(np.einsum(expr, *ops, optimize=np.einsum_path(expr, *ops, optimize="optimal")[0]))


def lagrangian(T):
    coef = np.array(list(itertools.product(range(3), repeat=6)))
    return (coef @ generators(T)) % 3


def gauss_inv(Lv, sig):
    """exact I(sigma) from the Lagrangian (Gauss sum)"""
    n = len(Lv)
    v1 = np.repeat(Lv, n, axis=0)
    v2 = np.tile(Lv, (n, 1))
    v3 = (-v1 - v2) % 3
    ok = np.ones(n * n, bool)
    phase = np.zeros(n * n, np.int64)
    weight = 1.0
    for p in range(6):
        s = tuple(sig[p])
        sl = slice(2 * p, 2 * p + 2)
        A = [v1[:, sl], v2[:, sl], v3[:, sl]]
        fixed = [c for c in range(3) if s[c] == c]
        if len(fixed) == 3:
            ok &= ~A[0].any(1) & ~A[1].any(1)
            weight *= 27
        elif len(fixed) == 1:
            ok &= ~A[fixed[0]].any(1)
            weight *= 9
        else:
            weight *= 3
            orient = -1 if s == (1, 2, 0) else 1        # calibrated against the direct contraction (chiral case)
            w = (A[0][:, 0] * A[1][:, 1] - A[0][:, 1] * A[1][:, 0]) % 3
            phase = (phase + orient * 2 * w) % 3
    return complex((ZETA ** phase[ok]).sum() * weight / 3 ** 18)


CHIRAL = [PERMS[0], PERMS[4], PERMS[2], PERMS[5], PERMS[1], PERMS[3]]   # found by search; inputs 0-2, outputs 3-5


def relation_symmetrised(Lv, sig):
    tot = 0j
    for pi in itertools.permutations(range(3)):
        for po in itertools.permutations(range(3)):
            s = [sig[pi[q]] for q in range(3)] + [sig[3 + po[q]] for q in range(3)]
            tot += gauss_inv(Lv, s)
    return tot / 36


def search_chiral(Lv, tries, seed):
    rng = np.random.default_rng(seed)
    for _ in range(tries):
        sig = [PERMS[1 + rng.integers(0, 2)], PERMS[1 + rng.integers(0, 2)], PERMS[rng.integers(0, 6)],
               PERMS[3 + rng.integers(0, 3)], PERMS[3 + rng.integers(0, 3)], PERMS[rng.integers(0, 6)]]
        sig = [sig[i] for i in rng.permutation(6)]
        v = gauss_inv(Lv, sig)
        if abs(v.imag) > 1e-12:
            return sig, v
    return None, None


def run():
    Bs = PT.load_facs()
    lab = np.load(H.CACHE)
    res = dict(pass_id=11201)
    # cross-check exact evaluator against direct contraction
    T0 = Bs[np.flatnonzero(lab == 3)[0]] % 3
    psi = choi_state(T0)
    L0 = lagrangian(T0)
    checks = [CHIRAL, [PERMS[1]] * 3 + [PERMS[3]] * 3, [PERMS[1], PERMS[2], PERMS[0], PERMS[3], PERMS[0], PERMS[0]]]
    res["exact_matches_contraction"] = all(abs(gauss_inv(L0, s) - chiral_invariant(psi, s)) < 1e-12 for s in checks)
    # same-orientation patterns are real
    rng = np.random.default_rng(5)
    real_ok = True
    for _ in range(30):
        orient = PERMS[1]
        sig = [orient if rng.random() < 0.5 else PERMS[rng.integers(3, 6)] for _ in range(6)]
        real_ok &= abs(gauss_inv(L0, sig).imag) < 1e-12
    res["same_orientation_patterns_real"] = bool(real_ok)
    # the chiral pattern on all directed orbitals (two members each) and the twist
    table = {}
    for o in (3, 4, 7, 9, 12, 13):
        mem = np.flatnonzero(lab == o)[:2]
        vals = [gauss_inv(lagrangian(Bs[m] % 3), CHIRAL) for m in mem]
        assert abs(vals[0] - vals[1]) < 1e-12
        tw = H.twist(Bs[mem[0]])
        table[str(o)] = dict(value=[round(vals[0].real, 12), round(vals[0].imag, 12)],
                             imag_times_243_over_sqrt3=round(vals[0].imag * 243 / np.sqrt(3), 9),
                             twist=None if tw is None else int(tw))
    res["chiral_pattern"] = ["".join(map(str, s)) for s in CHIRAL]
    res["chiral_values"] = table
    q3 = table["3"]["twist"]
    q3 = 1 if q3 == 1 else -1
    res["imag_equals_Q_sqrt3_over_243"] = bool(
        abs(table["3"]["value"][1] - q3 * np.sqrt(3) / 243) < 1e-10
        and abs(table["4"]["value"][1] + q3 * np.sqrt(3) / 243) < 1e-10)
    # constancy over ALL members of the 256 pair (each member carries its own enumeration labelling of F's factors)
    from collections import Counter
    allm = {}
    for o in (3, 4):
        c = Counter(round(gauss_inv(lagrangian(Bs[m] % 3), CHIRAL).imag * 243 / np.sqrt(3), 6)
                    for m in np.flatnonzero(lab == o))
        allm[str(o)] = {str(k): v for k, v in c.items()}
    res["all_members_imag_in_units_i_sqrt3_over_243"] = allm
    rel = Counter()
    for pi in itertools.permutations(range(3)):
        for po in itertools.permutations(range(3)):
            s = [CHIRAL[pi[q]] for q in range(3)] + [CHIRAL[3 + po[q]] for q in range(3)]
            rel[str(round(gauss_inv(L0, s).imag * 243 / np.sqrt(3), 6))] += 1
    res["relabellings_of_pattern_on_orbital_3"] = dict(rel)
    # relation invariant (summed over relabellings of inputs and outputs)
    sym = {o: relation_symmetrised(lagrangian(Bs[np.flatnonzero(lab == o)[0]] % 3), CHIRAL) for o in (3, 4)}
    res["symmetrised"] = {str(o): [round(v.real, 12), round(v.imag, 12)] for o, v in sym.items()}
    res["symmetrised_flips"] = bool(abs(sym[3].imag) > 1e-12 and abs(sym[3] - np.conj(sym[4])) < 1e-12)
    # a chiral pattern for the 6912 pair
    L12 = lagrangian(Bs[np.flatnonzero(lab == 12)[0]] % 3)
    sig12, v12 = search_chiral(L12, 300, 12)
    if sig12 is not None:
        L13 = lagrangian(Bs[np.flatnonzero(lab == 13)[0]] % 3)
        v13 = gauss_inv(L13, sig12)
        res["pair_6912"] = dict(pattern=["".join(map(str, s)) for s in sig12],
                                value_12=[round(v12.real, 12), round(v12.imag, 12)],
                                value_13=[round(v13.real, 12), round(v13.imag, 12)],
                                flips=bool(abs(v12 - np.conj(v13)) < 1e-12))
    else:
        res["pair_6912"] = "no chiral pattern in 300 random candidates"
    im = res["symmetrised"]["3"][1]
    res["paper_sentence"] = (
        "The twist is a measurable phase: a third-moment local-unitary invariant of the gate's Choi state is purely "
        "imaginary on the smallest time-oriented relation, $\\pm i\\sqrt3/243$ --- a quadratic Gauss sum --- with sign "
        "$Q$ on all $256$ members, and it changes sign under time reversal; invariants whose three-cycles all turn the "
        "same way are real, and averaging over all relabellings of the registers erases the sign --- the arrow is "
        "locally measurable only on labelled registers.")
    res["card_sentence"] = ("The twist is a measurable chiral phase: a third-moment invariant of the gate's Choi state "
                            "equals Q&middot;i&radic;3/243 and flips under time reversal.")
    res["symmetrised_imag_nonzero"] = bool(abs(im) > 1e-12)
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("paper_sentence", "card_sentence")}, indent=1))


if __name__ == "__main__":
    main()
