import sys,itertools,numpy as np
sys.path.insert(0,'analysis')
import w33_pass11557_pin_equivariant_event_dirac as P
gam,grade=P.small_clifford(); gam=[np.array(g.evalf(30),complex) for g in gam]; G0=np.array(grade.evalf(30),complex)
def shifts(L=3,flux=0):
 sites=list(itertools.product(range(L),repeat=3));idx={x:i for i,x in enumerate(sites)};Ts=[]
 for mu in range(3):
  T=np.zeros((L**3,L**3),complex)
  for i,x in enumerate(sites):
   phase=1
   if mu==0: phase=np.exp(-2j*np.pi*flux*x[1]/L**2)
   if mu==1 and x[1]==L-1: phase=np.exp(2j*np.pi*flux*x[0]/L)
   y=list(x);y[mu]=(y[mu]+1)%L;T[i,idx[tuple(y)]]=phase
  Ts.append(T)
 return Ts
def ov(L=3,flux=0,m0=1.,r=1.):
 T=shifts(L,flux);n=L**3;I=np.eye(n);Dn=np.zeros((4*n,4*n),complex);W=np.zeros((n,n),complex)
 for mu in range(3):
  Dn+=np.kron((T[mu]-T[mu].conj().T)/2,gam[mu])
  W+=I-(T[mu]+T[mu].conj().T)/2
 K=Dn+np.kron(r*W-m0*I,np.eye(4));G=np.kron(np.eye(n),G0);H=G@K
 herm=np.linalg.norm(H-H.conj().T);w,V=np.linalg.eigh(H);S=(V*np.sign(w))@V.conj().T;D=np.eye(4*n)+G@S
 gw=np.linalg.norm(G@D+D@G-D@G@D);index=-np.trace(S).real/2
 return herm,gw,index,min(abs(w)),np.count_nonzero(abs(np.linalg.eigvalsh(D.conj().T@D))<1e-9)
for f in [0,1,-1]:
 for m in [.5,1.,1.5,2.5]:
  print(f,m,ov(flux=f,m0=m))
