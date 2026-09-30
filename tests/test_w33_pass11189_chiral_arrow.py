"""Regression for Pass 11189: time reversal swaps exactly the tied reversal pairs; the profile-oriented Maslov spectrum is
chiral exactly there; the stabilizer Bargmann phase is (-i)^Maslov.  Frozen data plus live checks on the stored
representatives of the four chiral orbitals and a few Bargmann triples."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11189_chiral_arrow as C  # noqa: E402

D = json.loads((ROOT / "data" / "w33_pass11189_chiral_arrow.json").read_text())


def test_frozen():
    assert D['orbitals'] == 20 and D['tau_fixed'] == 16
    assert D['reversal_pairs'] == [[5, 11], [6, 7], [9, 10]]
    assert D['tau_swapped_pairs'] == [[6, 7], [9, 10]] and D['tied_pairs_are_mirror_pairs']
    assert D['spectrum_constant_on_orbitals'] and D['achiral_iff_tau_fixed'] and D['complete_classification']
    assert sorted(D['chirality_by_size'].values()) == [-54, -18, 18, 54]
    assert D['bargmann']['agree'] == D['bargmann']['total'] == 200


def test_live_chirality_and_mirror():
    rows = {r['orbital']: r for r in D['rows']}
    for k in (6, 9):
        B = np.array(rows[k]['rep'], np.int64)
        mu = C.spectrum(B)
        assert C.chirality(mu) == rows[k]['chirality'] != 0
        # the mirror image (tau B tau is an adapted basis of tau F) has the negated spectrum
        Bt = (C.TAU @ B @ C.TAU) % 3
        assert C.spectrum(Bt) == C.negated(mu)


def test_live_maslov_axioms_and_bargmann():
    rng = np.random.default_rng(7)
    P0 = C.product_lagrangians(np.eye(6, dtype=np.int64))
    B = np.array(D['rows'][9]['rep'], np.int64)
    PF = C.product_lagrangians(B)
    x, y, z = P0[rng.integers(0, 64, 200)], P0[rng.integers(0, 64, 200)], PF[rng.integers(0, 64, 200)]
    k1 = C.kashiwara(x, y, z)
    assert (k1 == C.kashiwara(y, z, x)).all() and ((k1 + C.kashiwara(y, x, z)) % 4 == 0).all()
    for t in range(6):
        c = int(C.kashiwara(x[t][None], y[t][None], z[t][None])[0])
        pr = [C.stabilizer_projector(L) for L in (x[t], y[t], z[t])]
        tr = np.trace(pr[0] @ pr[1] @ pr[2])
        assert np.isclose(tr / abs(tr), (-1j) ** c)
