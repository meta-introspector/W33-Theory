"""Pass 10980: can flavour structure suppress the FI-regenerated RPV below the proton bound at the Higgs-quartic ceiling?
Z6-I, 23 D-flat models, vacuum = Pass 10960 D-flat witness (contains the forced nu^c-like singlet).
Rules: gauge + space group (Z6-I defines no R-rules), which OVER-permit couplings, so every coupling below is an
UPPER bound on its true size.  Couplings ~ eps^(order-3) with random O(1) coefficients (eps = 0.3).
Light states: exotic vector-like pairs removed via their mass matrices; H_u = light bl; u1, q1 from the up-Yukawa
SVD.  Proton-relevant P = max over d^c flavours j,k and light lepton directions of |lambda''(u1 dj dk)| *
|lambda'(q1 L dk)|  (conservative: the down basis is undetermined when H_d is absent).  Bound at m_sq = 3e10 GeV:
|lambda' lambda''| < 1e-27 (m/100 GeV)^2 ~ 1e-10."""
import importlib.util, json, re, sys
from fractions import Fraction as F
from math import lcm
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
ROOT = __import__("pathlib").Path(r"C:\Repos\Theory of Everything")
SPEC = importlib.util.spec_from_file_location("m", ROOT / "analysis" / "w33_pass10967_fi_vacuum_regenerates_rpv.py")
M = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(M)
L, _ = M.P.load_ledger(); D = M.load(M.DISC)
C10960 = json.loads((ROOT / "data" / "w33_pass10960_matter_even_dflat_closure.json").read_text())["models"]
EPS, TRIALS = 0.3, 5
models = [n for n, r in C10960.items() if r.get("dflat_any") and r.get("bl_direction") and n.startswith("Z6-I|")]
out = {}
for name in sorted(models):
    d = D[name.split("|")[1]]
    orders, left = M.fields(L[name], d)
    nq = len(left[0]["q"])
    vec = {f["name"]: f["q"] + f["ext"] for f in left}
    lab = lambda p: sorted((f["name"] for f in left if f["base"] == p), key=lambda s: int(s.split("_")[1]))
    _w = C10960[name]["dflat_witness"]; S = list((eval(_w) if isinstance(_w, str) else _w).keys())
    dim = len(vec[S[0]])
    den = lcm(*[x.denominator for v in vec.values() for x in v])
    A = np.array([[int(x * den) for x in vec[s]] for s in S], dtype=float).T
    n, ncong = len(S), dim - nq
    cache = {}

    def order(fields):
        key = tuple(sorted(fields))
        if key in cache:
            return cache[key]
        t = np.array([int(sum(vec[f][i] for f in fields) * den) for i in range(dim)], dtype=float)
        Aeq = np.zeros((dim, n + ncong)); Aeq[:, :n] = A
        for j in range(ncong):
            Aeq[nq + j, n + j] = -den
        r = milp(np.concatenate([np.ones(n), np.zeros(ncong)]), constraints=[LinearConstraint(Aeq, -t, -t)],
                 integrality=np.ones(n + ncong), bounds=Bounds(np.concatenate([np.zeros(n), -1e7 * np.ones(ncong)]),
                                                             np.concatenate([40 * np.ones(n), 1e7 * np.ones(ncong)])),
                 options={"time_limit": 20})
        o = int(round(r.fun)) + len(fields) if r.status == 0 else None
        cache[key] = o
        return o

    rng = np.random.default_rng(1)
    coef = lambda: rng.normal() + 1j * rng.normal()
    mag = lambda o: 0.0 if o is None else EPS ** max(o - 3, 0)

    def mass(Xb, X):
        R_, C_ = lab(Xb), lab(X)
        Mx = np.array([[coef() * (0 if order([a, b]) is None else EPS ** (order([a, b]) - 2)) for b in C_] for a in R_])
        return R_, C_, Mx

    def light(Xb, X, left_side=False):
        R_, C_, Mx = mass(Xb, X)
        if not R_:
            return (C_, np.eye(len(C_))) if not left_side else (R_, np.eye(len(R_)))
        if left_side:
            u, s, vh = np.linalg.svd(Mx.T)
            r = int((s > 1e-9 * max(s.max(), 1e-300)).sum())
            return R_, vh[r:].conj().T
        u, s, vh = np.linalg.svd(Mx)
        r = int((s > 1e-9 * max(s.max(), 1e-300)).sum())
        return C_, vh[r:].conj().T

    try:
        Qn, Qv = light("bq", "q")
        Un, Uv = light("u", "bu")
        Dn, Dv = light("d", "bd")
        Ln, Lv = light("bl", "l")
        Hn, Hv = light("bl", "l", left_side=True)
    except Exception as e:
        out[name] = dict(error=str(e)); continue
    if Hv.shape[1] == 0:
        out[name] = dict(note="no light H_u"); continue
    Hu = Hv[:, 0]
    # up Yukawa in light basis
    Yu = np.zeros((len(Qn), len(Un)), complex)
    for i, a in enumerate(Qn):
        for j, b in enumerate(Un):
            Yu[i, j] = sum(coef() * mag(order([a, b, h])) * Hu[k] for k, h in enumerate(Hn))
    Yl = Qv.T @ Yu @ Uv
    u_, s_, vh_ = np.linalg.svd(Yl)
    q1 = Qv @ u_[:, -1].conj()          # lightest left-handed up combination in field basis
    u1 = Uv @ vh_[-1].conj()            # lightest right-handed up combination
    # lambda'' (u1 dj dk) and lambda' (q1 L dk) in the light d^c basis (max over flavours, conservative)
    Dl = Dv
    lam2 = np.zeros((Dl.shape[1], Dl.shape[1]), complex)
    T2 = {}
    for a in Un:
        for b in Dn:
            for c in Dn:
                T2[(a, b, c)] = coef() * mag(order([a, b, c]))
    for j in range(Dl.shape[1]):
        for k in range(Dl.shape[1]):
            lam2[j, k] = sum(u1[ia] * Dl[ib, j] * Dl[ic, k] * T2[(a, b, c)]
                             for ia, a in enumerate(Un) for ib, b in enumerate(Dn) for ic, c in enumerate(Dn))
    lam1 = np.zeros((Lv.shape[1], Dl.shape[1]), complex)
    for li in range(Lv.shape[1]):
        for k in range(Dl.shape[1]):
            lam1[li, k] = sum(q1[iq] * Lv[il, li] * Dl[ik, k] * coef() * mag(order([a, l_, c]))
                              for iq, a in enumerate(Qn) for il, l_ in enumerate(Ln) for ik, c in enumerate(Dn))
    P = max(abs(lam2[j, k]) * abs(lam1[li, k]) for j in range(lam2.shape[0]) for k in range(lam2.shape[1]) for li in range(lam1.shape[0]))
    out[name] = dict(P=float(P), yukawa_sv=[float(x) for x in s_], light=(Qv.shape[1], Uv.shape[1], Dv.shape[1], Lv.shape[1], Hv.shape[1]),
                     max_lam2=float(np.abs(lam2).max()), max_lam1=float(np.abs(lam1).max()))
    print(name.split("|")[1], f"P = {P:.2e}  max|l''| {np.abs(lam2).max():.2e}  max|l'| {np.abs(lam1).max():.2e}  "
          f"light (q,u,d,l,Hu) {out[name]['light']}", flush=True)
OUT = ROOT / "data" / "w33_pass10980_rpv_flavour_upper_bound.json"
vals0 = [v["P"] for v in out.values() if "P" in v]
summ = {"models": len(out), "with_P": len(vals0), "no_light_Hu": sum(1 for v in out.values() if v.get("note") == "no light H_u"),
        "five_light_d": sum(1 for v in out.values() if tuple(v.get("light", (0, 0, 0)))[2] == 5),
        "min_P": min(vals0) if vals0 else None, "below_1e-10": sum(x < 1e-10 for x in vals0)}
OUT.write_text(json.dumps(dict(pass_id=10980, eps=EPS, summary=summ, models=out), indent=1, sort_keys=True))
vals = [v["P"] for v in out.values() if "P" in v]
print("models with P:", len(vals), " min P:", min(vals) if vals else None, " below 1e-10:", sum(v < 1e-10 for v in vals))
