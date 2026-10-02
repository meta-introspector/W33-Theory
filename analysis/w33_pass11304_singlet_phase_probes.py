"""E6-singlet vectorlike family probes detect the old four Takagi flats."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11294_hierarchical_cp_phase_vacuum as F

def fields(t,e=.22):
    U=np.diag(np.array([1.,2.,3.])*np.exp(1j*np.r_[t[:2],-sum(t[:2])]))
    Q=F.mixing(e);D=Q@np.diag(np.array([1.,3.,5.])*np.exp(1j*np.r_[t[2:],-sum(t[2:])]))@Q.T
    return U,D

def moments(t,e=.22):
    U,D=fields(t,e);out=[]
    for k in (1.,2.):
        A=U+k*D;S=A@A.conj().T
        out.extend(np.trace(np.linalg.matrix_power(S,p)).real for p in (1,2,3))
    return np.array(out)

def cw(t,e=.22):
    U,D=fields(t,e);value=0.;spectra=[]
    for k in (1.,2.):
        A=.2*(U+k*D).conj();M=np.block([[A,3*np.eye(3)],[3*np.eye(3),np.zeros((3,3))]])
        w=np.linalg.eigvalsh(M.conj().T@M);assert min(w)>0
        value-=2*np.sum(w*w*(np.log(w/9)-1.5))/(64*np.pi**2);spectra.append(w.tolist())
    return float(value),spectra

def payload():
    controls=[];h=2e-5
    for e in (.18,.22,.3,.4):
        J=np.column_stack([(moments(np.eye(4)[i]*h,e)-moments(-np.eye(4)[i]*h,e))/(2*h) for i in range(4)])
        J/=np.linalg.norm(J,axis=1)[:,None];sv=np.linalg.svd(J,compute_uv=False);assert sv[-1]>1e-4
        controls.append({'epsilon':e,'scaled_phase_singular_values':sv.tolist(),'rank':int(sum(sv>1e-7))})
    a,spectra=cw(np.zeros(4));b,_=cw(np.array([.1,-.07,.13,.04]));assert abs(a-b)>1e-7
    return {'status':'PASS','result_scope':'PASS_SM_PRESERVING_PHASE_SENSITIVE_SINGLET_PROBE_CONSTRUCTION',
      'action':'For k=1,2 add E6-singlet vectorlike Weyls chi_k(3)+chibar_k(bar3) under familySU3. L has M chi chibar+(y/2) chi^T(Fu†+k Fd†)chi+h.c. All y,k,M are real; the action is CP even and leaves the E6/SM stabilizer intact.',
      'mass_map':'The symmetric6x6 mass is [[y(Fu†+kFd†),M I],[M I,0]]. Its full determinant detects mixed phases. M=3,y=.2 are inputs.',
      'invariant_map':'Six moments Tr[((Fu+kFd)(Fu+kFd)†)^p], k=1,2 and p=1,2,3, have rank4 on the old four determinant-preserving Takagi phase directions at all four reference vacua.',
      'phase_controls':controls,'one_loop_potential_control':{'reference':a,'phase_displaced':b,'difference':b-a,'reference_squared_masses':spectra},
      'gauge_cost':'Two vectorlike family3+bar3 pairs have zero anomaly and cost4/3 family beta units, but zero E6/SO10 beta units. Family gauging remains a UV issue.',
      'boundaries':['A named renormalizable interaction and actual phase-sensitive fermion loop, not a derivation of the observed CP phase or hierarchy.','Local identifiability is not a positive Hessian at a stationary quantum vacuum. Probe masses, couplings and matching coefficients are inputs.','The sextet Gram targets and epsilon used in controls are prior imposed data; no measured CKM prediction.'],
      'prior_owners':['analysis/w33_pass11299_e6_flavor_spectral_obstruction.py','analysis/w33_pass11294_hierarchical_cp_phase_vacuum.py'],
      'primary_sources':['https://arxiv.org/abs/1103.2915','https://arxiv.org/abs/2008.08606']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11304_singlet_phase_probes.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])
