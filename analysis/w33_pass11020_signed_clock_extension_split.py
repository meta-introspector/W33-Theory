#!/usr/bin/env python3
"""Pass 11020: the signed E6 cubic has a hidden GL2(3) monomial symmetry.

Pass 10975 showed that all 48 GL2(3) support operations preserve the cubic
support and each admits 64 diagonal sign repairs, while only the identity
preserves the frozen signs by a bare permutation.

This packet proves that the 64-fold repair fibres can be chosen coherently:
the signed monomial extension contains a subgroup of order 48 projecting
bijectively to the support GL2(3). The extension therefore splits.

A second exact check shows why this does not trivialize the clock reduction:
no signed lift of the central -I preserves the four six-point fibre-constant
clock subspace.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "analysis") not in sys.path:
    sys.path.insert(0, str(ROOT / "analysis"))

import w33_pass10947_five_front_execution as P47
import w33_pass10975_e6_cubic_reynolds_clock_weld as P75

OUT = ROOT / "data" / "w33_pass11020_signed_clock_extension_split.json"
P10975 = ROOT / "data" / "w33_pass10975_e6_cubic_reynolds_clock_weld.json"
def comp(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inv(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def xor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def pull(p, s):
    return tuple(s[p[i]] for i in range(len(p)))


def perm_order(p):
    ident = tuple(range(len(p)))
    x = ident
    for n in range(1, 200):
        x = comp(x, p)
        if x == ident:
            return n
    raise AssertionError("permutation order bound")


def perm_closure(gens, n, cap=5000):
    ident = tuple(range(n))
    group = {ident}
    front = [ident]
    while front:
        a = front.pop()
        for b in gens:
            c = comp(a, b)
            if c not in group:
                group.add(c)
                front.append(c)
                if len(group) > cap:
                    return group
    return group
def rref_nullspace(matrix):
    m = [row[:] for row in matrix]
    nr = len(m)
    nc = len(m[0])
    pivots = []
    r = 0
    for c in range(nc):
        i = next((i for i in range(r, nr) if m[i][c]), None)
        if i is None:
            continue
        m[r], m[i] = m[i], m[r]
        for j in range(nr):
            if j != r and m[j][c]:
                m[j] = [x ^ y for x, y in zip(m[j], m[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(nc) if c not in pivots]
    basis = []
    for f in free:
        v = [0] * nc
        v[f] = 1
        for rr, c in enumerate(pivots):
            if m[rr][f]:
                v[c] = 1
        basis.append(tuple(v))
    return r, tuple(basis)


def in_kernel(A, s):
    return all(sum(row[i] * s[i] for i in range(len(s))) % 2 == 0 for row in A)
def signed_setup():
    signed, e6_to_h, parent, affine, p72 = P75.load_inputs()
    dirs = tuple(tuple(map(int, d)) for d in parent["carrier"]["labels"])
    rows = P75.support_action_packet(signed, e6_to_h, dirs, parent)
    triads = tuple(sorted(signed))
    perms = tuple(r[3] for r in rows)
    assert len(set(perms)) == 48

    A = []
    for tri in triads:
        row = [0] * 27
        for i in tri:
            row[i] = 1
        A.append(row)
    rank, kernel_basis = rref_nullspace(A)
    assert (rank, len(kernel_basis)) == (21, 6)

    def particular(p):
        equations = []
        for tri, row in zip(triads, A):
            image = tuple(sorted(p[i] for i in tri))
            rhs = 0 if signed[image] == signed[tri] else 1
            equations.append(row + [rhs])
        r, nullity, sol = P47.solve_linear_mod(equations, 2)
        assert (r, nullity) == (21, 6)
        return tuple(sol)

    return signed, e6_to_h, dirs, rows, triads, perms, A, kernel_basis, particular
def mon_comp(x, y):
    """x after y for monomial maps e_i -> (-1)^s_i e_p(i)."""
    p, s = x
    q, t = y
    return comp(p, q), xor(t, pull(q, s))


def mon_order(x):
    p, _ = x
    ident = (tuple(range(len(p))), (0,) * len(p))
    y = ident
    for n in range(1, 5000):
        y = mon_comp(x, y)
        if y == ident:
            return n
    raise AssertionError("monomial order bound")


def mon_closure(gens, n, cap=5000):
    ident = (tuple(range(n)), (0,) * n)
    group = {ident}
    front = [ident]
    while front:
        a = front.pop()
        for b in gens:
            c = mon_comp(a, b)
            if c not in group:
                group.add(c)
                front.append(c)
                if len(group) > cap:
                    return group
    return group


def all_sign_lifts(p, particular, kernel_basis):
    s0 = particular(p)
    out = []
    for coeff in itertools.product((0, 1), repeat=len(kernel_basis)):
        s = s0
        for bit, k in zip(coeff, kernel_basis):
            if bit:
                s = xor(s, k)
        out.append((p, s))
    assert len(set(out)) == 64
    return tuple(out)
def preserves_signed_cubic(element, signed, triads):
    p, s = element
    for tri in triads:
        image = tuple(sorted(p[i] for i in tri))
        parity = sum(s[i] for i in tri) % 2
        expected = signed[tri] if parity == 0 else -signed[tri]
        if signed[image] != expected:
            return False
    return True


def find_support_generators(perms):
    for a in sorted(perms):
        for b in sorted(perms):
            if len(perm_closure((a, b), 27, cap=100)) == 48:
                return a, b
    raise AssertionError("no two-generator support pair")


def find_split_pair(g0, g1, particular, kernel_basis):
    lifts0 = all_sign_lifts(g0, particular, kernel_basis)
    lifts1 = all_sign_lifts(g1, particular, kernel_basis)
    size_hist = Counter()
    tested = 0
    for a in lifts0:
        for b in lifts1:
            tested += 1
            h = mon_closure((a, b), 27, cap=400)
            size_hist[len(h)] += 1
            if len(h) == 48 and len({p for p, _ in h}) == 48:
                return a, b, h, tested, size_hist, lifts0, lifts1
    raise AssertionError("signed extension did not split")
def fibrewise_constant_sign(s, e6_to_h, dirs):
    fibres = []
    for d in dirs:
        fibres.append([
            i for i, h in e6_to_h.items()
            if h[:2] != (0, 0) and P75.norm_dir(h[:2]) == d
        ])
    return all(len({s[i] for i in fibre}) == 1 for fibre in fibres)


def payload():
    parent = json.loads(P10975.read_text(encoding="utf-8"))
    assert parent["signed_gauge_firewall"]["repairs_per_support_operation"] == 64

    (
        signed, e6_to_h, dirs, rows, triads, perms,
        A, kernel_basis, particular,
    ) = signed_setup()

    support_hist = Counter(perm_order(p) for p in perms)
    assert support_hist == Counter({1: 1, 2: 13, 3: 8, 4: 6, 6: 8, 8: 12})

    # Kernel normality under every support pullback.
    kernel_normal = all(
        in_kernel(A, pull(p, k))
        for p in perms
        for k in kernel_basis
    )
    assert kernel_normal

    g0, g1 = find_support_generators(perms)
    split0, split1, section, tested, size_hist, lifts0, lifts1 = find_split_pair(
        g0, g1, particular, kernel_basis
    )
    assert len(section) == 48
    assert {p for p, _ in section} == set(perms)
    assert all(preserves_signed_cubic(x, signed, triads) for x in section)

    p_to_matrix = {r[3]: tuple(r[0]) for r in rows}
    gen_matrices = [list(p_to_matrix[g0]), list(p_to_matrix[g1])]
    gen_support_orders = [perm_order(g0), perm_order(g1)]
    gen_lift_orders = [mon_order(split0), mon_order(split1)]
    assert gen_support_orders == gen_lift_orders

    extension_order = len(perms) * (2 ** len(kernel_basis))
    assert extension_order == 3072

    minus_row = next(r for r in rows if tuple(r[0]) == (2, 0, 0, 2))
    pminus = minus_row[3]
    minus_lifts = all_sign_lifts(pminus, particular, kernel_basis)
    fibrewise_minus = sum(
        fibrewise_constant_sign(s, e6_to_h, dirs)
        for _, s in minus_lifts
    )
    assert fibrewise_minus == 0

    section_minus = next(x for x in section if x[0] == pminus)
    assert mon_order(section_minus) == 2
    section_fibrewise = sum(
        fibrewise_constant_sign(s, e6_to_h, dirs)
        for _, s in section
    )
    assert section_fibrewise == 1

    checks = {
        "parent_64_repairs":
            parent["signed_gauge_firewall"]["repairs_per_support_operation"] == 64,
        "support_GL23_order48": len(perms) == 48,
        "sign_kernel_dimension6": len(kernel_basis) == 6,
        "sign_kernel_order64": 2 ** len(kernel_basis) == 64,
        "kernel_normal_under_support": kernel_normal,
        "full_signed_extension_order3072": extension_order == 3072,
        "two_support_generators_generate48":
            len(perm_closure((g0, g1), 27, cap=100)) == 48,
        "split_pair_found": len(section) == 48,
        "split_projection_bijective": {p for p, _ in section} == set(perms),
        "split_section_preserves_signed_cubic":
            all(preserves_signed_cubic(x, signed, triads) for x in section),
        "generator_orders_preserved": gen_support_orders == gen_lift_orders,
        "minusI_has_64_repairs": len(minus_lifts) == 64,
        "minusI_no_fibrewise_constant_repair": fibrewise_minus == 0,
        "split_section_only_identity_preserves_fibre_constant_module":
            section_fibrewise == 1,
    }
    assert all(checks.values())

    return {
        "schema": "w33.pass11020.signed-clock-extension-split.v1",
        "status": "PASS",
        "headline": (
            "The 64-fold diagonal sign-repair fibres over the 48 exact GL2(3) "
            "clock-support operations form a 3072-element signed monomial extension "
            "that splits: an explicit 48-element subgroup projects bijectively to "
            "GL2(3) and preserves the frozen signed E6 cubic exactly."
        ),
        "extension": {
            "support_group": "GL2(3)",
            "support_order": 48,
            "support_element_order_histogram": {
                str(k): int(v) for k, v in sorted(support_hist.items())
            },
            "sign_kernel": "(C2)^6",
            "sign_kernel_dimension": 6,
            "sign_kernel_order": 64,
            "full_monomial_extension_order": extension_order,
            "exact_sequence": "1 -> (C2)^6 -> E_3072 -> GL2(3) -> 1",
            "splits": True,
            "split_subgroup_order": len(section),
        },
        "explicit_split_witness": {
            "support_generator_matrices_mod3": gen_matrices,
            "support_generator_orders": gen_support_orders,
            "lift_generator_orders": gen_lift_orders,
            "lift_generator_sign_vectors": [
                list(split0[1]), list(split1[1])
            ],
            "lift_pairs_tested_before_first_split": tested,
            "closure_size_histogram_before_first_split": {
                str(k): int(v) for k, v in sorted(size_hist.items())
            },
            "section_sign_weight_histogram": {
                str(k): int(v)
                for k, v in sorted(
                    Counter(sum(s) for _, s in section).items()
                )
            },
        },
        "clock_reduction_firewall": {
            "central_support_element": "-I in GL2(3)",
            "projective_action_on_four_clock_labels": "trivial",
            "signed_repairs_checked": len(minus_lifts),
            "fibrewise_constant_signed_repairs": fibrewise_minus,
            "split_section_minusI_order": mon_order(section_minus),
            "split_section_elements_preserving_fibre_constant_module":
                section_fibrewise,
            "conclusion": (
                "The exact signed GL2(3) symmetry acts on the full 27-coordinate "
                "cubic carrier but does not descend to the naive four-dimensional "
                "fibre-constant clock field. Every lift of central -I introduces "
                "nonconstant signs inside at least one six-point fibre."
            ),
        },
        "what_changed": (
            "Pass 10975's sign-repair statement is upgraded from elementwise "
            "compatibility to a coherent group action: the obstruction is gauge, "
            "not a nonsplit extension. The remaining obstruction is instead the "
            "reduction from the 27-coordinate signed carrier to four fibre amplitudes."
        ),
        "reynolds_consequence": (
            "The S4 Reynolds projection of Pass 10975 averages a genuine exact "
            "support-symmetry shadow of a signed GL2(3) monomial action. Because "
            "that signed action does not preserve the fibre-constant field, the "
            "Reynolds reduction remains an effective projection rather than the "
            "restriction of the 27-dimensional signed action."
        ),
        "boundary": (
            "This is a finite monomial-symmetry theorem. It does not derive a "
            "Hamiltonian that performs group averaging, an energy scale, a continuum "
            "field, or a physical vacuum. The split section is an explicit witness, "
            "not asserted to be unique or dynamically preferred."
        ),
        "checks": checks,
    }
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path, default=OUT)
    a = ap.parse_args()
    p = payload()
    text = json.dumps(p, indent=2, sort_keys=True) + "\n"
    if a.check:
        if not a.output.exists() or a.output.read_text(encoding="utf-8") != text:
            raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": p["status"],
        "extension_order": p["extension"]["full_monomial_extension_order"],
        "split": p["extension"]["splits"],
        "section_order": p["extension"]["split_subgroup_order"],
        "minusI_fibrewise_repairs":
            p["clock_reduction_firewall"]["fibrewise_constant_signed_repairs"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
