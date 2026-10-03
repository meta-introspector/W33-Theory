"""Renormalizable healthy Gaussian-adjoint mediation: exact flavor alignment obstruction."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1]

def eliminate(A,B,M,alpha,beta):
 n=len(A);A=A-np.trace(A)*np.eye(n)/n;B=B-np.trace(B)*np.eye(n)/n;source=alpha[:,None,None]*A+beta[:,None,None]*B;X=-np.einsum('ab,bij->aij',np.linalg.inv(M),source);energy=.5*np.einsum('ab,aij,bji->',M,X,X).real+np.einsum('aij,aji->',X,source).real
 return X,float(energy)
def payload():
 rng=np.random.default_rng(11336);n=3;Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));A=Z@Z.conj().T;Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));B=Z@Z.conj().T;M=np.array([[2.,.3],[.3,3.]]);alpha=np.array([.4,-.2]);beta=np.array([.1,.5]);X,v=eliminate(A,B,M,alpha,beta);chi=float(alpha@np.linalg.solve(M,beta));A0=A-np.trace(A)*np.eye(3)/3;B0=B-np.trace(B)*np.eye(3)/3
 exact=-.5*(alpha@np.linalg.solve(M,alpha))*np.trace(A0@A0).real-.5*(beta@np.linalg.solve(M,beta))*np.trace(B0@B0).real-chi*np.trace(A0@B0).real;assert abs(v-exact)<1e-12
 T=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));T=(T+T.conj().T)/2;h=1e-5;Up=expm(1j*h*T);Um=expm(-1j*h*T);fd=(eliminate(A,Up@B@Up.conj().T,M,alpha,beta)[1]-eliminate(A,Um@B@Um.conj().T,M,alpha,beta)[1])/(2*h);analytic=(-1j*chi*np.trace((B@A-A@B)@T)).real;assert abs(fd-analytic)<1e-7
 return {'status':'PASS','scope':'Exact tree-level no-go for any number of positive-mass Gaussian family-adjoint mediators coupled only to traceless UUdag and DDdag; not every renormalizable mediator theory.','action':'V=Vsep(spectraA,spectraB)+.5 M2_ab TrXaXb+TrXa(alpha_a A0+beta_a B0), M2 positive definite; canonical positive kinetic terms. Cubic Xa UUdag couplings are power-counting renormalizable.','elimination':'Xa=-(M2^-1)_ab(alpha_b A0+beta_b B0); Veff=Vsep-.5 aa TrA0²-.5 bb TrB0²-chi TrA0B0, chi=alphaT M2^-1 beta.','alignment_proof':'For relative orientation B->exp(i tT)Bexp(-i tT), dV=-i chi Tr([B,A]T). Stationarity for all Hermitian tracelessT and chi!=0 forces[A,B]=0. Jarlskog commutator vanishes, even with many mediators. Ifchi0 every relative orientation is flat and cannot be uniquely selected by this action.','positive_mass_eigenvalues':np.linalg.eigvalsh(M).tolist(),'chi':chi,'energy_replay_error':abs(v-exact),'orientation_derivative_error':abs(fd-analytic),'heavy_stationarity_residual':float(max(abs(np.einsum('ab,bij->aij',M,X)+alpha[:,None,None]*A0+beta[:,None,None]*B0).flat)),'implication':'This healthy renormalizable subclass cannot derive the11326 chiral squared invariant. Non-Gaussian interactions, extra sources or controlled loop-generated orientation terms must be built explicitly. A sufficiently strong spectral quartic can bound the combined potential, but no UV-complete AF model is asserted.','prior_owners':['analysis/w33_pass11326_noncommuting_flavor_operator.py','analysis/w33_pass11280_factor_higgs_mediator.py'],'primary_sources':['https://arxiv.org/abs/hep-ph/0211440']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11336_gaussian_mediator_alignment.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['chi'])
