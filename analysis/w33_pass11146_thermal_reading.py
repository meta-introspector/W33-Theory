#!/usr/bin/env python3
"""Pass 11146: the neutral exit read thermally -- the Hagedorn instability needs the Wilson-line holonomy, and the
thermofield double is the Choi state of Euclidean time.

(1) No holonomy, no Hagedorn.  With the Wilson lines switched off, the Witten-twisted sector of all 21 neutral-exit
survivors has NO tachyon at any radius sampled (Im T = 1.5 ... 0.2 on the diagonal of both Wilson-line tori, windings
scaled as 2.5/y); with them on, every model is tachyonic at small radius (data/w33_pass11146_no_holonomy_21.json, scan
analysis/w33_pass11146_scan_no_holonomy.py).  Lattice reason: level matching l^2 + (pL^2 - pR^2) = 1 with l in
E8xE8 + V0 needs l^2 = 1 for a pure winding (pL^2 = pR^2), and E8xE8 + V0 has no norm-1 vector with the beta phase;
winding-momentum states (pL^2 - pR^2 = 2mn != 0) have pR^2 >= 1.  The Wilson line W supplies l = pi + V0 + N W with
l^2 = 1 (Pass 11134).  Read thermally (a circle with a (-1)^F-type twist is a Euclidean time circle; the Wilson line is
an imaginary chemical potential), the neutral exit is a Hagedorn transition that EXISTS ONLY because of the chemical
potential; its critical curve is the Gamma_0(3) disk of Pass 11136, and beyond it the T-dual description is the U(16)
tachyonic string (Pass 11140).
(2) The thermofield double is the Choi state of Euclidean time.  For a qutrit Gibbs state rho = e^{-beta K}/Z,
|TFD> = sum_i sqrt(p_i)|ii>, and its partial transpose equals the Jamiolkowski matrix of the Euclidean half-evolution
A -> rho^(1/2) A rho^(1/2) -- verified to 1e-14 below.  So the spatial entanglement of the two copies IS the temporal
(pseudo-density) correlation of one qutrit across imaginary time beta/2, with negativity ((sum_i sqrt p_i)^2 - 1)/2:
1 at beta = 0 (the temporal Bell line of Pass 11143, the identity channel) falling to 0 as beta -> infinity (a product
line).  With the modular clock of the paper's Section 3 (K = -log rho), "thermal time" (Connes-Rovelli) is the flow
whose Euclidean continuation prepares this entanglement.
Interpretation (NOT computed here): at the Hagedorn point a winding condensate on the thermal circle is the
Horowitz-Polchinski string star, the string-scale end of a black hole, whose Lorentzian continuation is an eternal
black hole with an Einstein-Rosen bridge between the two TFD copies (ER = EPR).  In that reading the neutral condensate
would be a string-scale bridge built from entanglement across Euclidean time.  The internal Wilson-line torus is not
spacetime time; this is a formal dictionary (double Wick rotation), stated as such.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NOHOL = ROOT / "data" / "w33_pass11146_no_holonomy_21.json"
OUT = ROOT / "data" / "w33_pass11146_thermal_reading.json"


def tfd_checks(energies=(0.0, 1.0, 2.0), betas=(0.0, 0.3, 1.0, 3.0, 10.0)):
    rows = []
    for b in betas:
        p = np.exp(-b * np.array(energies)); p = p / p.sum()
        psi = np.zeros(9); psi[[0, 4, 8]] = np.sqrt(p)
        rho2 = np.outer(psi, psi)
        pt = rho2.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
        s = np.diag(np.sqrt(p))
        J = np.zeros((9, 9))
        for i in range(3):
            for j in range(3):
                Eji = np.zeros((3, 3)); Eji[j, i] = 1
                Eij = np.zeros((3, 3)); Eij[i, j] = 1
                J += np.kron(Eij, s @ Eji @ s)
        ev = np.linalg.eigvalsh(pt)
        rows.append(dict(beta=b, pt_equals_jamiolkowski_half_evolution=float(np.abs(pt - J).max()),
                         negativity=float((np.abs(ev).sum() - 1) / 2),
                         formula=float((np.sqrt(p).sum() ** 2 - 1) / 2)))
    return rows


def summarize():
    res = dict(pass_id=11146, tfd=tfd_checks())
    if NOHOL.exists():
        d = json.loads(NOHOL.read_text())
        res['no_tachyon_without_holonomy'] = all(v[y]['without_wl'] == 0 for v in d.values() for y in v)
        res['tachyonic_with_holonomy_at_small_radius'] = all(any(v[y]['with_wl'] > 0 for y in v) for v in d.values())
        res['n_models'] = len(d)
    OUT.write_text(json.dumps(res, indent=1, sort_keys=True))
    return res


if __name__ == "__main__":
    print(json.dumps(summarize(), indent=1))
