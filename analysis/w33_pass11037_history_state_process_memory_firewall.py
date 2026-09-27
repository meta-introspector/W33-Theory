#!/usr/bin/env python3
"""Pass 11037: Page-Wootters cubic history state versus genuine process memory."""
from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/w33_pass11037_history_state_process_memory_firewall.json"


def entropy(vals):
    return -sum(x*math.log2(x) for x in vals if x>0)


def partial_transpose(rho,dA,dB):
    return rho.reshape(dA,dB,dA,dB).transpose(2,1,0,3).reshape(dA*dB,dA*dB)


def negativity(rho,dA,dB):
    ev=np.linalg.eigvalsh(partial_transpose(rho,dA,dB))
    return float(-ev[ev<0].sum())


def bell(d=3):
    v=np.zeros(d*d,complex)
    for j in range(d): v[j*d+j]=1/math.sqrt(d)
    return v


def isotropic(eta,d=3):
    v=bell(d); pure=np.outer(v,v.conj())
    return eta*pure+(1-eta)*np.eye(d*d)/(d*d)
def payload():
    # CCZ eigenvalue census for f(a,b,c)=abc mod 3.
    counts=[0,0,0]
    for a,b,c in itertools.product(range(3),repeat=3):
        counts[(a*b*c)%3]+=1
    assert counts==[19,4,4]
    overlap=(counts[0]-(counts[1]+counts[2])/2)/27
    assert abs(overlap-5/9)<1e-15

    rho_clock=np.full((3,3),5/27,dtype=float)
    np.fill_diagonal(rho_clock,1/3)
    eig=np.linalg.eigvalsh(rho_clock)[::-1]
    expected=np.array([19/27,4/27,4/27])
    assert np.max(np.abs(eig-expected))<1e-14
    S=entropy(eig)
    assert abs(S-1.1730125852957602)<1e-12

    # Maximal history input: equal amplitudes from eigenvalue sectors 1,w,w^2.
    # Its three orbit states have Gram matrix I_3.
    w=np.exp(2j*np.pi/3)
    amps=np.array([1,1,1],complex)/math.sqrt(3)
    gram=np.empty((3,3),complex)
    for s in range(3):
      for t in range(3):
        gram[s,t]=sum(abs(amps[j])**2*w**(j*(t-s)) for j in range(3))
    assert np.max(np.abs(gram-np.eye(3)))<1e-12
    # Causal-break memory witness.
    phi=bell(3); bellrho=np.outer(phi,phi.conj())
    quantum_neg=negativity(bellrho,3,3)
    classical=np.zeros((9,9),complex)
    for j in range(3):
        v=np.zeros(9); v[j*3+j]=1
        classical += np.outer(v,v)/3
    classical_neg=negativity(classical,3,3)
    assert abs(quantum_neg-1)<1e-12 and abs(classical_neg)<1e-12

    # Exact depolarizing threshold for a qutrit memory channel.
    samples={}
    for eta in (0,0.2,0.25,0.3,0.5,1):
        n=negativity(isotropic(eta),3,3)
        formula=max(0,(4*eta-1)/3)
        assert abs(n-formula)<1e-12
        samples[str(eta)]=n

    # Phase-damping on the Bell coherences gives negativity exactly lambda.
    deph={}
    for lam in (0,0.1,0.5,1):
        rho=np.zeros((9,9),complex)
        for i in range(3):
          for j in range(3):
            rho[i*3+i,j*3+j]=(1 if i==j else lam)/3
        n=negativity(rho,3,3)
        assert abs(n-lam)<1e-12
        deph[str(lam)]=n
    checks={
      "ccz_eigenspace_counts_19_4_4":counts==[19,4,4],
      "symmetric_input_overlap_5_over_9":abs(overlap-5/9)<1e-15,
      "clock_spectrum_19_4_4_over_27":np.max(np.abs(eig-expected))<1e-14,
      "maximal_history_input_gram_identity":np.max(np.abs(gram-np.eye(3)))<1e-12,
      "quantum_hidden_memory_negativity_one":abs(quantum_neg-1)<1e-12,
      "classical_measure_prepare_negativity_zero":abs(classical_neg)<1e-12,
      "depolarizing_entanglement_threshold_one_quarter":samples["0.25"]<1e-12 and samples["0.3"]>0,
      "dephasing_negativity_equals_coherence":all(abs(float(k)-v)<1e-12 for k,v in deph.items()),
    }
    assert all(checks.values())
    checks={k:bool(v) for k,v in checks.items()}
    return {
      "schema":"w33.pass11037.history-state-process-memory-firewall.v1",
      "status":"PASS",
      "headline":(
        "The cubic CCZ history state has exact relational-time observables but "
        "those are not quantum-memory witnesses. For |+>^3 the qutrit clock "
        "spectrum is (19,4,4)/27 and entropy is about 1.173 bits; an equal "
        "superposition of the three CCZ eigensectors makes the clock maximally "
        "mixed. Across an operational causal break, however, a preprogrammed "
        "history with no hidden memory has zero reference negativity, while a "
        "hidden qutrit memory can retain negativity 1."
      ),
      "relational_history":{
        "ccz_eigenspace_dimensions":counts,
        "symmetric_input_overlap":"5/9",
        "clock_reduced_spectrum":["19/27","4/27","4/27"],
        "clock_entropy_bits":S,
        "maximal_input_clock_spectrum":["1/3","1/3","1/3"],
        "maximal_entropy_bits":math.log2(3),
        "stationarity_convention":"(X_clock tensor U_CCZ)|Psi> = |Psi>",
      },
      "causal_break_firewall":{
        "no_hidden_memory":(
          "Trace-and-prepare on the exposed qutrit is entanglement breaking; "
          "subsequent preprogrammed local gates cannot restore reference-system "
          "entanglement."
        ),
        "quantum_store_retrieve_negativity":quantum_neg,
        "classical_measure_prepare_negativity":classical_neg,
        "witness":"reference negativity surviving the causal break",
      },
      "noise_thresholds":{
        "qutrit_depolarizing_model":"N(eta)=max(0,(4*eta-1)/3)",
        "entanglement_survives_iff":"eta > 1/4",
        "sample_negativities":samples,
        "phase_damping_model":"N(lambda)=lambda",
        "phase_damping_samples":deph,
      },
      "interpretation":(
        "Clock-system entanglement certifies a relational history encoding; the "
        "nested cubic echo certifies holonomy; reference entanglement across a "
        "causal break certifies quantum temporal memory. These are three distinct "
        "resources and should not be conflated."
      ),
      "boundary":(
        "The hidden-memory and classical controls are exact finite protocols, "
        "not a derivation that the universe implements such a memory. A full "
        "multi-time Choi/process-tensor implementation should additionally freeze "
        "all causal partial-trace constraints for a chosen laboratory circuit."
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
      if not a.output.exists() or a.output.read_text()!=text:
        raise SystemExit("certificate drift")
    else:
      a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text)
    print(json.dumps({"status":p["status"],
      "clock_entropy":p["relational_history"]["clock_entropy_bits"],
      "quantum_neg":p["causal_break_firewall"]["quantum_store_retrieve_negativity"],
      "depol_threshold":p["noise_thresholds"]["entanglement_survives_iff"]},sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
