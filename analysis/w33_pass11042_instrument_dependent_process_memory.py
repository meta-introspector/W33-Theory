#!/usr/bin/env python3
"""Pass 11042: causal one-slot qutrit comb with instrument-dependent memory visibility."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11042_instrument_dependent_process_memory.json"
D=3


def idx(f,o,i):
    return (f*D+o)*D+i


def classical_comb():
    Y=sp.zeros(D**3)
    for h in range(D):
        for o in range(D):
            Y[idx(h,o,h),idx(h,o,h)]=sp.Rational(1,D)
    return Y


def coherent_comb():
    Y=sp.zeros(D**3)
    for o in range(D):
        for a in range(D):
            for b in range(D):
                Y[idx(a,o,a),idx(b,o,b)]=sp.Rational(1,D)
    return Y


def trace_F(Y):
    T=sp.zeros(D*D)
    for o in range(D):
        for i in range(D):
            for op in range(D):
                for ip in range(D):
                    T[o*D+i,op*D+ip]=sum(
                        Y[idx(f,o,i),idx(f,op,ip)] for f in range(D))
    return T
def partial_transpose_F(A):
    B=np.zeros_like(A,dtype=complex)
    for f in range(D):
      for o in range(D):
       for i in range(D):
        for fp in range(D):
         for op in range(D):
          for ip in range(D):
           B[idx(fp,o,i),idx(f,op,ip)]=A[idx(f,o,i),idx(fp,op,ip)]
    return B


def negativity_process(Y):
    # Normalize the comb to trace one only for this entanglement diagnostic.
    A=np.array(Y.tolist(),dtype=complex)
    A=A/np.trace(A)
    ev=np.linalg.eigvalsh(partial_transpose_F(A))
    return float(-ev[ev< -1e-12].sum())


def instrument_distribution(kind):
    # Hidden h is copied to past P and future F. The exposed middle system is |h>.
    # Instrument is a projective measurement followed by a fixed replacement.
    out=np.zeros((D,D,D),float) # p,m,f
    w=np.exp(2j*np.pi/D)
    for h in range(D):
      for m in range(D):
        if kind=="Z":
            prob=1.0 if h==m else 0.0
        elif kind=="X":
            amp=w**(-m*h)/math.sqrt(D)
            prob=abs(amp)**2
        else:
            raise ValueError(kind)
        out[h,m,h]+=prob/D
    assert abs(out.sum()-1)<1e-12
    return out


def conditional_mutual_information(PMF):
    val=0.0
    pm=PMF.sum(axis=(0,2))
    for m in range(D):
        if pm[m]==0: continue
        q=PMF[:,m,:]/pm[m]
        qp=q.sum(axis=1); qf=q.sum(axis=0)
        I=0.0
        for p in range(D):
            for f in range(D):
                if q[p,f]>0:
                    I+=q[p,f]*math.log2(q[p,f]/(qp[p]*qf[f]))
        val+=pm[m]*I
    return val
def comb_outcome_prob(Y,kind):
    # J_m = |0><0|_O tensor E_m^T_I.
    w=np.exp(2j*np.pi/D)
    probs=[]
    for m in range(D):
        E=np.zeros((D,D),complex)
        if kind=="Z":
            E[m,m]=1
        else:
            v=np.array([w**(m*j)/math.sqrt(D) for j in range(D)],complex)
            E=np.outer(v,v.conj())
        J=np.zeros((D*D,D*D),complex)
        # fixed output |0>
        for i in range(D):
            for j in range(D):
                J[i,j]=E.T[i,j]  # rows/cols correspond (o=0,i)
        A=np.array(Y.tolist(),dtype=complex)
        total=0j
        for f in range(D):
            sl=[idx(f,o,i) for o in range(D) for i in range(D)]
            total+=np.trace(A[np.ix_(sl,sl)]@J)
        probs.append(float(total.real))
    return probs


def payload():
    Yc=classical_comb(); Yq=coherent_comb()
    target=sp.eye(D*D)/D
    assert trace_F(Yc)==target
    assert trace_F(Yq)==target
    assert all(x>=0 for x in Yc.eigenvals())
    assert all(x>=0 for x in Yq.eigenvals())

    nc=negativity_process(Yc); nq=negativity_process(Yq)
    assert abs(nc)<1e-12 and abs(nq-1)<1e-12

    Pz=instrument_distribution("Z")
    Px=instrument_distribution("X")
    Iz=conditional_mutual_information(Pz)
    Ix=conditional_mutual_information(Px)
    assert abs(Iz)<1e-12
    assert abs(Ix-math.log2(3))<1e-12

    pz=comb_outcome_prob(Yc,"Z")
    px=comb_outcome_prob(Yc,"X")
    assert np.max(np.abs(np.array(pz)-1/3))<1e-12
    assert np.max(np.abs(np.array(px)-1/3))<1e-12
    checks={
      "classical_comb_positive":True,
      "coherent_comb_positive":True,
      "classical_causal_trace":trace_F(Yc)==target,
      "coherent_causal_trace":trace_F(Yq)==target,
      "comb_trace_is_output_dimension":sp.trace(Yc)==sp.trace(Yq)==D,
      "classical_process_negativity_zero":abs(nc)<1e-12,
      "coherent_process_negativity_one":abs(nq-1)<1e-12,
      "Z_instrument_screens_past_future":abs(Iz)<1e-12,
      "X_instrument_exposes_log2_3_memory":abs(Ix-math.log2(3))<1e-12,
      "both_instruments_normalized_uniform_outcomes":True,
    }
    assert all(checks.values())
    checks={k:bool(v) for k,v in checks.items()}
    return {
      "schema":"w33.pass11042.instrument-dependent-process-memory.v1",
      "status":"PASS",
      "headline":(
        "A fully positive, causally normalized qutrit one-slot process already "
        "shows instrument-dependent memory visibility. For a hidden classical "
        "trit copied from past to future, a Z-basis middle instrument screens "
        "past from future exactly, I(P:F|M)=0, while a Fourier-basis instrument "
        "reveals no hidden label and leaves I(P:F|M)=log2(3). Causal order is "
        "unchanged in both cases."
      ),
      "process_combs":{
        "spaces":"F(final qutrit) x O(intervention output) x I(intervention input)",
        "dimension":27,
        "causal_constraint":"Tr_F Upsilon = I_O tensor (I_I/3)",
        "trace":3,
        "classical_hidden_memory":{
          "formula":"(1/3) sum_h |h><h|_F tensor I_O tensor |h><h|_I",
          "normalized_F_vs_OI_negativity":nc,
        },
        "coherent_hidden_memory":{
          "formula":"(1/3) |Omega_un><Omega_un|_(F,I) tensor I_O",
          "normalized_F_vs_OI_negativity":nq,
        },
      },
      "instrument_dependence":{
        "hidden_process":"P=H, middle exposed state |H>, F=H with H uniform",
        "Z_measure_and_reprepare":{
          "outcome_probabilities":pz,
          "conditional_mutual_information_bits":Iz,
          "reading":"M reveals H, so conditioning on M screens P from F",
        },
        "Fourier_measure_and_reprepare":{
          "outcome_probabilities":px,
          "conditional_mutual_information_bits":Ix,
          "reading":"M is independent of H, so P and F remain perfectly correlated given M",
        },
        "difference_bits":Ix-Iz,
      },
      "conceptual_result":(
        "This is an exact finite example of instrument-specific Markov order. "
        "Memory disappearing under one intervention family does not mean the "
        "underlying process is globally memoryless, and memory reappearing under "
        "another family does not imply retrocausality."
      ),
      "connection_to_W33":(
        "Passes 11038-11039 provide candidate W33 operator frames for such instrument "
        "comparisons, but the raw off-Bell rays need CP completion before they can "
        "replace the physical Z/Fourier instruments used in this control."
      ),
      "boundary":(
        "This is the minimal one-slot comb, not a full many-slot process-tensor "
        "reconstruction. It proves the causal/instrument-dependence principle; the "
        "next hardware packet should construct CP W33-derived instruments and repeat "
        "the conditional-memory test on the signed 24-mode carrier."
      ),
      "checks":checks,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    a=ap.parse_args(); p=payload()
    text=json.dumps(p,indent=2,sort_keys=True)+"\n"
    if a.check:
        if not a.output.exists() or a.output.read_text()!=text: raise SystemExit("certificate drift")
    else:
        a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],
      "I_Z":p["instrument_dependence"]["Z_measure_and_reprepare"]["conditional_mutual_information_bits"],
      "I_X":p["instrument_dependence"]["Fourier_measure_and_reprepare"]["conditional_mutual_information_bits"],
      "quantum_neg":p["process_combs"]["coherent_hidden_memory"]["normalized_F_vs_OI_negativity"]},
      sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
