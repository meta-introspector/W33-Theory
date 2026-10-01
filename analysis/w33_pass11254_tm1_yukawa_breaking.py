"""Pass 11254: a Yukawa model on the line stabiliser's S4 -- what the sequestering coupling eps does to TM1.

Pass 11243: TM1 = the neutrino symmetry swaps the charged-lepton point X with a second point Q of the same W33 line;
Pass 11248: in a sequestered two-flavon potential the sign of eps in eps (phi . chi)^2 selects TM1 (eps < 0), and the
flavons tilt away from their symmetric directions at O(eps).  This pass adds the Yukawa layer and asks what is fixed.

SCOPE FIRST.  In residual-symmetry TM1, theta13 is NOT a function of eps: at eps = 0 the Z2 fixes one column, and
theta13, theta23 and delta inside the remaining 2x2 block are set by the Yukawa couplings.  What eps does is BREAK the
TM1 relations (|U_e1|^2 = 2/3 and the sum rule tying cos delta to theta13, theta23).  That is what is computed.

MODEL (lepton doublets L in the 3' of S4, flavons phi, chi in 3'; S4-equivariant mass matrices, checked numerically):
  M_e^dag M_e = H_e(phi) = alpha |phi|^2 + beta phi phi^T + gamma diag(phi_i^2) + i kappa1 [phi]_x + i kappa3 [phi^3]_x
      ([v]_x the cross-product matrix; the antisymmetric terms are what split m_mu from m_tau under the residual Z3)
  M_nu(chi) = a + b |chi|^2 + c chi chi^T + d diag(chi_i^2) + e K(chi)   (Weinberg operator, complex a..e)
      K = the unique equivariant symmetric cubic; WITHOUT it every quadratic M_nu at the chord commutes with the
      non-S4 reflection diag(-1,1,1), which freezes all of U_PMNS (theta13 fixed): an accidental symmetry
(not the most general equivariant forms -- a representative family at low degree).
At eps = 0 (symmetric vacua) H_e commutes with the Z3 fixing d_X and M_nu with the transposition (XQ): the first PMNS
column is exactly TM1 for EVERY choice of couplings (checked).  alpha, beta, kappa1 are solved so that the lightest
eigenvector is the Z3 eigenvector carrying the TM1 weight 2/3 and m_e : m_mu : m_tau are physical; gamma, kappa3 random.
The neutrino couplings are fitted to NuFIT's sin^2 theta13 and Delta m21^2/Delta m31^2 (normal ordering), with the
lightest mass m1^2/Delta m31^2 drawn in [0, 0.15] (Sum m_nu <~ 0.08 eV).
Then the exact eps-tilted vacua of Pass 11248 are inserted and the mixing is recomputed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import w33_pass11230_tm1_alignment as A  # noqa: E402
import w33_pass11248_sequestered_tm1 as SQ  # noqa: E402

OUT = ROOT / "data" / "w33_pass11254_tm1_yukawa_breaking.json"
ME, MMU, MTAU = 0.000511, 0.10566, 1.77686
R_DM = 7.49e-5 / 2.513e-3                                    # Delta m21^2 / Delta m31^2 (NuFIT 6.0, NO)
S13_RANGE = (0.02030, 0.02388)
S12_RANGE = (0.275, 0.345)
G = A.G3P


def cross(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def he_terms(phi):
    p2 = phi @ phi
    return [p2 * np.eye(3), np.outer(phi, phi), np.diag(phi ** 2), 1j * cross(phi), 1j * cross(phi ** 3)]


def kcubic(chi):
    """the unique S4-equivariant symmetric cubic (Reynolds projection of all 60 cubic monomial maps has rank 1)"""
    x, y, z = chi
    a, b, c = z * (y * y - x * x), y * (x * x - z * z), x * (z * z - y * y)
    return np.array([[0, a, b], [a, 0, c], [b, c, 0]])


def mnu_terms(chi):
    return [np.eye(3), (chi @ chi) * np.eye(3), np.outer(chi, chi), np.diag(chi ** 2), kcubic(chi)]


def equivariance_check(rng):
    p = rng.normal(size=3)
    ok = True
    for g in G:
        for A_, B_ in zip(he_terms(g @ p), he_terms(p)):
            ok &= np.allclose(A_, g @ B_ @ g.T)
        for A_, B_ in zip(mnu_terms(g @ p), mnu_terms(p)):
            ok &= np.allclose(A_, g @ B_ @ g.T)
    return bool(ok)


def pmns(He, Mnu):
    ev, Ue = np.linalg.eigh(He)                               # ascending: e, mu, tau
    m2, Un = np.linalg.eigh(Mnu.conj().T @ Mnu)               # ascending: nu1, nu2, nu3
    return Ue.conj().T @ Un, ev, m2


def angles(U):
    s13 = abs(U[0, 2]) ** 2
    c13 = 1 - s13
    s12 = abs(U[0, 1]) ** 2 / c13
    s23 = abs(U[1, 2]) ** 2 / c13
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    col = (abs(U[0, 0]) ** 2, abs(U[1, 0]) ** 2, abs(U[2, 0]) ** 2)
    cd = A.cos_delta(col, s13, s23)
    return dict(s12=s12, s13=s13, s23=s23, J=J, Ue1=col[0], cos_delta=None if cd is None else cd[1])


def tm1_sum_rule_cos_delta(s13, s23):
    out = A.cos_delta((2 / 3, 1 / 6, 1 / 6), s13, s23)
    return None if out is None else out[1]


def setup_charged(phi0, rng, perm):
    """solve alpha, beta, kappa1 so that H_e(phi0) has the physical masses on the right Z3 eigenvectors"""
    z3 = [g for g in G if np.allclose(g @ phi0, phi0) and not np.allclose(g, np.eye(3))
          and np.allclose(np.linalg.matrix_power(g, 3), np.eye(3))][0]
    w, V = np.linalg.eig(z3)
    V = np.linalg.qr(V)[0]
    gamma, kappa3 = rng.normal(), rng.normal()
    terms = he_terms(phi0)
    D = np.array([[np.real(V[:, i].conj() @ T @ V[:, i]) for T in terms] for i in range(3)])   # (3 eigvecs, 5 terms)
    return V, gamma, kappa3, D, perm


def run(n_models=60, seed=11254):
    rng = np.random.default_rng(seed)
    res = dict(pass_id=11254, equivariant=equivariance_check(rng))
    eps_list = [-0.0005, -0.001, -0.0025, -0.005, -0.01, -0.02, -0.04]
    raw = {e: SQ.global_min(e).x for e in eps_list}
    x_ref = raw[eps_list[0]]
    phi_s, chi_s, _, _ = SQ.nearest_symmetric(x_ref)
    phi0 = phi_s / np.linalg.norm(phi_s) * np.linalg.norm(x_ref[:3])
    chi0 = chi_s / np.linalg.norm(chi_s) * np.linalg.norm(x_ref[3:])

    def align(x):
        # the potential is S4-invariant AND even in phi, chi separately; the Yukawa sector is not even in phi,
        # so every eps-vacuum is mapped to the branch continuously connected to (phi0, chi0)
        best, bx = -np.inf, None
        for g in G:
            for sf in (1, -1):
                for sc in (1, -1):
                    f, c = sf * g @ x[:3], sc * g @ x[3:]
                    s = f @ phi0 / np.linalg.norm(f) / np.linalg.norm(phi0) + c @ chi0 / np.linalg.norm(c) / np.linalg.norm(chi0)
                    if s > best:
                        best, bx = s, np.concatenate([f, c])
        return bx
    vacua = {e: align(x) for e, x in raw.items()}
    res["vacuum_tilt_deg"] = {str(e): [float(np.degrees(np.arccos(min(1, abs(v[:3] @ phi0) / np.linalg.norm(v[:3]) / np.linalg.norm(phi0))))),
                                       float(np.degrees(np.arccos(min(1, abs(v[3:] @ chi0) / np.linalg.norm(v[3:]) / np.linalg.norm(chi0)))))]
                              for e, v in vacua.items()}
    # the Z3 eigenvector with TM1 weight 2/3 against the (XQ) isolated eigenvector
    invol = [g for g in G if np.allclose(g @ chi0, chi0) and np.allclose(g @ g, np.eye(3)) and not np.allclose(g, np.eye(3))]
    ev_i, Vi = np.linalg.eigh(invol[0])
    iso = Vi[:, 0] if np.sum(ev_i < 0) == 1 else Vi[:, 2]
    models = []
    tries = 0
    rej = {k: 0 for k in ('singular_charged_solve', 'He_not_positive', 'neutrino_fit_failed', 'not_TM1_at_eps0')}
    while len(models) < n_models and tries < 40 * n_models:
        tries += 1
        V, gamma, kappa3, Dm, _ = setup_charged(phi0, rng, None)
        w2 = np.abs(V.conj().T @ iso) ** 2
        e_idx = int(np.argmax(w2))                          # the 2/3 eigenvector must be the electron
        others = [i for i in range(3) if i != e_idx]
        if rng.integers(2):
            others = others[::-1]
        targets = np.zeros(3)
        targets[e_idx], targets[others[0]], targets[others[1]] = ME ** 2, MMU ** 2, MTAU ** 2
        # unknowns alpha, beta, kappa1 (terms 0, 1, 3); gamma, kappa3 fixed
        rhs = targets - Dm[:, 2] * gamma - Dm[:, 4] * kappa3
        try:
            sol = np.linalg.solve(Dm[:, [0, 1, 3]], rhs)
        except np.linalg.LinAlgError:
            rej["singular_charged_solve"] += 1
            continue
        ce = np.array([sol[0], sol[1], gamma, sol[2], kappa3])
        He0 = sum(c * T for c, T in zip(ce, he_terms(phi0)))
        if np.min(np.linalg.eigvalsh(He0)) < 0:
            rej["He_not_positive"] += 1
            continue

        def Mnu(cn, chi):
            z = cn[:5] + 1j * cn[5:]
            return sum(c * T for c, T in zip(z, mnu_terms(chi)))

        m1sq = rng.uniform(0.0, 0.15)                   # m1^2 / Delta m31^2: Sum m_nu <~ 0.08 eV (cosmology-safe)

        def resid(cn):
            U, _, m2 = pmns(He0, Mnu(cn, chi0))
            a = angles(U)
            d31 = m2[2] - m2[0]
            return [(a["s13"] - 0.02195) / 0.0006, ((m2[1] - m2[0]) / d31 - R_DM) / 0.001, d31 - 1.0,
                    (m2[0] / d31 - m1sq) / 0.01]

        fit = least_squares(resid, rng.normal(size=10), xtol=1e-12, ftol=1e-12)
        if np.max(np.abs(fit.fun)) > 1e-3:
            rej["neutrino_fit_failed"] += 1
            continue
        cn = fit.x
        U0, _, _ = pmns(He0, Mnu(cn, chi0))
        a0 = angles(U0)
        if abs(a0["Ue1"] - 2 / 3) > 1e-9:
            rej["not_TM1_at_eps0"] += 1
            continue
        rows = []
        for e in eps_list:
            x = vacua[e]
            phi_t, chi_t = x[:3], x[3:]
            phi_s = phi0 / np.linalg.norm(phi0) * np.linalg.norm(phi_t)       # symmetric direction, same size
            chi_s = chi0 / np.linalg.norm(chi0) * np.linalg.norm(chi_t)
            out = {}
            for tag, (f, c) in dict(both=(phi_t, chi_t), nu_only=(phi_s, chi_t), e_only=(phi_t, chi_s)).items():
                He = sum(cc * T for cc, T in zip(ce, he_terms(f)))
                U, _, _ = pmns(He, Mnu(cn, c))
                out[tag] = angles(U)
            a = out["both"]
            sr = tm1_sum_rule_cos_delta(a["s13"], a["s23"])
            rows.append(dict(eps=e, dUe1=a["Ue1"] - 2 / 3, dUe1_nu_only=out["nu_only"]["Ue1"] - 2 / 3,
                             dUe1_e_only=out["e_only"]["Ue1"] - 2 / 3, s12=a["s12"], s13=a["s13"], s23=a["s23"],
                             cos_delta=a["cos_delta"], sum_rule_cos_delta=sr,
                             sum_rule_violation=None if (sr is None or a["cos_delta"] is None) else a["cos_delta"] - sr))
        models.append(dict(eps0=dict(s12=a0["s12"], s13=a0["s13"], s23=a0["s23"], cos_delta=a0["cos_delta"],
                                     m1sq_over_dm31=m1sq),
                           scan=rows))
    res["models"] = len(models)
    res["rejections"] = rej
    if not models:
        return res
    res["tries"] = tries
    # eps0 spread: theta23 and delta are free at eps = 0 (set by the couplings)
    res["eps0_s23_range"] = [min(m["eps0"]["s23"] for m in models), max(m["eps0"]["s23"] for m in models)]
    res["eps0_cos_delta_range"] = [min(m["eps0"]["cos_delta"] for m in models if m["eps0"]["cos_delta"] is not None),
                                   max(m["eps0"]["cos_delta"] for m in models if m["eps0"]["cos_delta"] is not None)]
    res["eps0_s12_tm1"] = float(np.median([m["eps0"]["s12"] for m in models]))
    # slopes at small eps
    slope_u = [m["scan"][0]["dUe1"] / abs(m["scan"][0]["eps"]) for m in models]
    slope_sr = [m["scan"][0]["sum_rule_violation"] / abs(m["scan"][0]["eps"]) for m in models
                if m["scan"][0]["sum_rule_violation"] is not None]
    res["dUe1_per_unit_eps"] = dict(median=float(np.median(slope_u)), p10=float(np.percentile(slope_u, 10)),
                                    p90=float(np.percentile(slope_u, 90)))
    th = {e: np.radians(res["vacuum_tilt_deg"][str(e)]) for e in eps_list}      # (phi tilt, chi tilt) in radians
    for tag, k in (("nu_only", 1), ("e_only", 0)):
        s = [abs(m["scan"][0]["dUe1_" + tag]) / th[eps_list[0]][k] for m in models]
        s2 = [abs(m["scan"][1]["dUe1_" + tag]) / th[eps_list[1]][k] for m in models]
        res[f"abs_dUe1_per_radian_tilt_{tag}"] = dict(median=float(np.median(s)), p10=float(np.percentile(s, 10)),
                                                      p90=float(np.percentile(s, 90)))
        res[f"linear_regime_{tag}"] = float(np.mean([abs(a / b - 1) < 0.1 for a, b in zip(s, s2) if b > 1e-12]))
    res["amplification_mtau2_over_mmu2"] = MTAU ** 2 / MMU ** 2
    res["sum_rule_violation_per_unit_eps"] = dict(median=float(np.median(slope_sr)),
                                                  p10=float(np.percentile(slope_sr, 10)),
                                                  p90=float(np.percentile(slope_sr, 90))) if slope_sr else None
    lin = []
    for m in models:
        s = [r["dUe1"] / r["eps"] for r in m["scan"][:3]]
        lin.append(all(abs(x / s[0] - 1) < 0.1 for x in s[1:]) if abs(s[0]) > 1e-12 else True)
    res["dUe1_linear_in_eps_fraction"] = float(np.mean(lin))
    # how large may |eps| be before s12 leaves the NuFIT 3 sigma window?
    within = {}
    for e_i, e in enumerate(eps_list):
        within[str(e)] = float(np.mean([S12_RANGE[0] <= m["scan"][e_i]["s12"] <= S12_RANGE[1] for m in models]))
    res["fraction_s12_in_3sigma"] = within
    res["examples"] = models[:3]
    return res


def main():
    res = run()
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1, default=float))


if __name__ == "__main__":
    main()
