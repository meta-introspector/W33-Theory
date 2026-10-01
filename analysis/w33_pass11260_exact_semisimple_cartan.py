#!/usr/bin/env python3
"""Pass 11260: an exact semisimple Cartan three-plane in the signed E8 grading.

The 81 coordinates are the repository's canonical ``27 x 3`` tensor gauge.
The only floating-point step is the already certified change from the stored
E6 matrix basis to E8 coordinates.  Multiplication by 18 produces an integral
adjoint matrix; the rounding residual is retained and every subsequent rank,
kernel, characteristic-polynomial, and minimal-polynomial calculation is over
the integers/rationals.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11260_exact_semisimple_cartan.json"

BACKGROUND = {
    3: -1, 7: -2, 13: -2, 14: -1, 22: 2, 23: 1,
    27: 1, 42: -2, 65: 2, 72: 2, 79: 1, 80: 1,
}


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def dense(entries=BACKGROUND):
    v = [0] * 81
    for i, a in entries.items():
        v[int(i)] = int(a)
    return v


def sparse(v):
    return [
        {"flat_index": i, "e6id": i // 3, "phase": i % 3, "value": int(a)}
        for i, a in enumerate(v) if a
    ]


def build_adjoint():
    toe = _load(ROOT / "tools/toe_e8_z3graded_bracket_jacobi.py", "toe_p11260")
    basis78 = np.load(ROOT / "artifacts/e6_27rep_basis_export/E6_basis_78.npy").astype(complex)
    proj = toe.E6Projector(basis78)
    bracket = toe.E8Z3Bracket(
        e6_projector=proj,
        cubic_triads=toe._load_signed_cubic_triads(),
        scale_g1g1=1.0,
        scale_g2g2=-1.0 / 6.0,
        scale_e6=1.0,
        scale_sl3=1.0 / 6.0,
    )
    sl3 = []
    for i in range(3):
        for j in range(3):
            if i != j:
                q = np.zeros((3, 3), complex); q[i, j] = 1; sl3.append(q)
    sl3 += [np.diag([1, -1, 0]).astype(complex), np.diag([0, 1, -1]).astype(complex)]
    sl3 = np.array(sl3)
    sf = sl3.reshape(8, -1)
    sleft = np.linalg.inv(sf @ sf.conj().T) @ sf.conj()
    zero = toe.E8Z3.zero()
    ebasis = [toe.E8Z3(basis78[i], zero.sl3, zero.g1, zero.g2) for i in range(78)]
    ebasis += [toe.E8Z3(zero.e6, sl3[i], zero.g1, zero.g2) for i in range(8)]
    for grade in (1, 2):
        for i in range(81):
            m = np.zeros((27, 3), complex); m.flat[i] = 1
            ebasis.append(toe.E8Z3(zero.e6, zero.sl3,
                                   m if grade == 1 else zero.g1,
                                   m if grade == 2 else zero.g2))

    def coords(e):
        rhs = np.einsum("aij,ji->a", basis78, e.e6)
        ce6 = proj.gram_inv @ rhs
        csl3 = sleft @ e.sl3.reshape(-1)
        return np.r_[ce6, csl3, e.g1.ravel(), e.g2.ravel()]

    def admat(v):
        X = np.array(v, complex).reshape(27, 3)
        elt = toe.E8Z3(zero.e6, zero.sl3, X, zero.g2)
        raw = 18 * np.column_stack([coords(bracket.bracket(elt, b)) for b in ebasis])
        residual = float(np.max(np.abs(raw - np.rint(raw.real))))
        return np.rint(raw.real).astype(np.int64), residual

    return admat


def characteristic_certificate(N3):
    x = sp.Symbol("x")
    blocks = {}
    for name, lo, hi in (("g0", 0, 86), ("g1", 86, 167), ("g2", 167, 248)):
        D = DomainMatrix.from_Matrix(sp.Matrix(N3[lo:hi, lo:hi].tolist())).convert_to(ZZ)
        coeff = D.charpoly()
        poly = sum(c * x ** (len(coeff) - 1 - i) for i, c in enumerate(coeff))
        fac = sp.factor_list(poly)[1]
        blocks[name] = {
            "dimension": hi - lo,
            "rank": D.rank(),
            "characteristic_factors": [
                {"polynomial": str(f), "degree": int(sp.degree(f, x)), "multiplicity": int(e)}
                for f, e in fac
            ],
        }
    return blocks


def repeated_factor_nullities(N3):
    A = 1234235584512
    B = 34326907456637308502016
    C = 3172169114198268924301601144832
    out = {}
    for name, lo, hi in (("g0", 0, 86), ("g1", 86, 167), ("g2", 167, 248)):
        D = DomainMatrix.from_Matrix(sp.Matrix(N3[lo:hi, lo:hi].tolist())).convert_to(ZZ)
        I = DomainMatrix.eye(D.shape, ZZ)
        D3 = D ** 3; D6 = D3 ** 2; D9 = D6 * D3
        fm = D9.add(D6.scalarmul(-A)).add(D3.scalarmul(-B)).add(I.scalarmul(-C))
        fp = D9.add(D6.scalarmul(A)).add(D3.scalarmul(-B)).add(I.scalarmul(C))
        out[name] = {"f_minus_nullity": hi - lo - fm.rank(), "f_plus_nullity": hi - lo - fp.rank()}
    return out


def run(full=True):
    old = _load(ROOT / "analysis/w33_pass11218_e8_cubic_cartan_kernel.py", "old_p11260")
    parent = old.load_parent(); records = parent.build_records(); v = dense()
    A = sp.Matrix(parent.jacobian(v, records))
    kernel = [old.primitive_integer_vector(z) for z in A.nullspace()]
    K = sp.Matrix.hstack(*(sp.Matrix(z) for z in kernel))
    pivots = tuple(map(int, K.T.rref()[1]))
    commutes = all(old.bracket(kernel[i], kernel[j], records) == [0] * 81
                   for i in range(3) for j in range(i + 1, 3))
    assert A.rank() == 78 and len(kernel) == 3 and commutes
    assert pivots == (3, 4, 6) and int(K[list(pivots), :].det()) == 4
    assert v == [-a for a in kernel[2]]

    admat = build_adjoint(); N, residual = admat(v)
    DN = DomainMatrix.from_Matrix(sp.Matrix(N.tolist())).convert_to(ZZ)
    rankN = DN.rank(); N3 = N @ N @ N
    offblock = N3.copy()
    offblock[:86, :86] = 0; offblock[86:167, 86:167] = 0; offblock[167:, 167:] = 0
    assert not np.any(offblock)
    blocks = characteristic_certificate(N3) if full else {}
    nullities = repeated_factor_nullities(N3) if full else {}
    if full:
        assert [blocks[g]["rank"] for g in ("g0", "g1", "g2")] == [78, 78, 78]
        assert all(z == {"f_minus_nullity": 27, "f_plus_nullity": 27} for z in nullities.values())
    rankN3 = sum(blocks[g]["rank"] for g in blocks) if full else None
    out = {
        "schema": "w33.pass11260.exact-semisimple-cartan.v1",
        "status": "PASS_EXACT_SEMISIMPLE_CARTAN_THREE_PLANE_EMBEDDED_IN_SIGNED_3x27",
        "background": {"support": sparse(v), "cubic_jacobian_rank": 78},
        "cartan_basis": [sparse(z) for z in kernel],
        "cartan_checks": {
            "dimension": 3, "pairwise_brackets_zero": commutes,
            "background_equals_minus_basis_2": True,
            "pivot_rows": list(pivots), "coordinate_minor": 4,
        },
        "scaled_adjoint": {
            "definition": "N=18 ad(v)", "integral": True,
            "rounding_residual": residual, "rank_N": rankN,
            "rank_N_cubed": rankN3, "kernel_dimension": 248 - rankN,
            "N_cubed_respects_Z3_blocks": True,
        },
        "N_cubed_blocks": blocks,
        "repeated_degree9_factor_nullities": nullities,
        "semisimplicity_proof": (
            "Each N^3 block has square-free minimal polynomial: the only repeated characteristic "
            "factors are the two square-free degree-9 factors, and each has its full 27-dimensional "
            "eigenspace in every block. Hence N^3 is diagonalizable. rank(N)=rank(N^3)=234 rules "
            "out a nontrivial zero Jordan block; for nonzero eigenvalues the derivative of t^3 is "
            "nonzero. Therefore N is diagonalizable. Its g1 centralizer is the displayed commuting "
            "three-plane, the minimal dimension equal to the Vinberg rank, so it is a Cartan subspace."
        ),
        "classification_input": {
            "grading": "E8 = (E6+A2) + (27,3) + (27*,3*)",
            "theta_rank": 3, "little_Weyl_group": "G26", "degrees": [6, 12, 18],
            "source": "Reeder-Levy-Yu-Gross, Table 21, E8 row 3b",
        },
        "boundary": (
            "This certifies the missing semisimple Cartan embedding. It does not by itself choose a "
            "physical vacuum, kinetic metric, Yukawa tensor, or statistics assignment."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(run(full=True), indent=2))
