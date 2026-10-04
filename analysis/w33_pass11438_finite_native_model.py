"""Analytic gradient of the explicitly supplied324-real-field bounded model.

Native invariant functions and compact generators are owned by11384 and the
20261001 native Cartan circuit. Uses supplied kappa=1e4,rho=1e-6,sigma=.1.
"""
from contextlib import contextmanager
import numpy as np
import w33_pass11428_11432_native_alignment_measure_coarse_flags as M
import w33_20261001_degree18_phase_completion as D
TENSORS=M.native_tensors()[1]
H=M.compact_generators(M.native_tensors()[0])
kappa=1e4
rho=1e-6
sigma=.1

@contextmanager
def native_context():
    original=D.I.tensors
    D.I.tensors=lambda:TENSORS
    try:
        yield
    finally:
        D.I.tensors=original

pack=lambda p,s:np.r_[p.real,p.imag,s.real,s.imag]
def native(v):
 a,b,ga,gb=D.I.evaluate(v);c,gc=D.evaluate(v);r=b-a*a;t=c-a**3
 acts=H@v;mu=(acts.conj()@v).real
 pot=.5*mu@mu+abs(a)**4-abs(a)**2+a.real**2+a.real/2+abs(r)**2+abs(t)**2
 hol=(2*a.conjugate()*(2*abs(a)**2-1)+2*a.real+.5)*ga+2*r.conjugate()*(gb-2*a*ga)+2*t.conjugate()*(gc-3*a*a*ga)
 md=2*np.einsum('a,ai->i',mu,acts)
 grad=np.r_[hol.real+md.real,-hol.imag+md.imag]
 return pot,grad,a,ga

def fun(x):
 p=(x[:81]+1j*x[81:162]);s=(x[162:243]+1j*x[243:]);v,gp,a,ga=native(p);w,gs,b,gb=native(s)
 pp=p.reshape(27,3);ss=s.reshape(27,3);norm=np.vdot(p,p).real+np.vdot(s,s).real;den=norm/2+sigma*sigma
 C=pp.T@ss.conj()/den;U=2*np.eye(3)+(C+C.conj().T)/2;B=np.eye(3)+C@C.conj().T;W=U@B-B@U;q=np.imag(np.trace(W@W@W))
 dcs=[]
 for fld,imag in [(0,False),(0,True),(1,False),(1,True)]:
  base=pp if fld==0 else ss
  for aa in range(27):
   for ii in range(3):
    dc=np.zeros((3,3),complex)
    if fld==0:dc[ii,:]=(1j if imag else 1)*ss[aa,:].conj()/den
    else:dc[:,ii]=(-1j if imag else 1)*pp[aa,:]/den
    dg=base[aa,ii].imag if imag else base[aa,ii].real
    dcs.append(dc-C*dg/den)
 dc=np.array(dcs);du=(dc+dc.conj().transpose(0,2,1))/2;db=dc@C.conj().T+C@dc.conj().transpose(0,2,1)
 dw=du@B-B@du+U@db-db@U;dq=3*np.imag(np.einsum('ij,aji->a',W@W,dw))
 chi=a.imag;sat=chi/(1+chi*chi);dchi=np.r_[ga.imag,ga.real,np.zeros(162)]
 grad=np.r_[gp,gs]-kappa*(sat*dq+q*(1-chi*chi)/(1+chi*chi)**2*dchi)+4*rho*norm*x
 value=v+w-kappa*sat*q+rho*norm*norm
 return float(value),grad
