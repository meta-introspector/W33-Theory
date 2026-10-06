import numpy as np,itertools
def gammas():
 sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1,-1]).astype(complex);I=np.eye(2)
 g=[np.kron(sx,s) for s in [sx,sy,sz]]+[np.kron(sy,I)]
 return g,g[0]@g[1]@g[2]@g[3]
def links(L,q,flux=(1,1)):
 sites=list(itertools.product(range(L),repeat=4));look={x:i for i,x in enumerate(sites)};Ts=[]
 for mu in range(4):
  T=np.zeros((L**4,L**4),complex)
  for i,x in enumerate(sites):
   phase=1+0j
   for first,m in [(0,flux[0]),(2,flux[1])]:
    if mu==first: phase*=np.exp(-2j*np.pi*q*m*x[first+1]/L**2)
    if mu==first+1 and x[first+1]==L-1: phase*=np.exp(2j*np.pi*q*m*x[first]/L)
   y=list(x);y[mu]=(y[mu]+1)%L;T[i,look[tuple(y)]]=phase
  Ts.append(T)
 return Ts
def ov(q,L=4,mass=1.):
 g,g5=gammas();T=links(L,q);n=L**4;G=np.kron(np.eye(n),g5)
 W=(4-mass)*np.eye(4*n,dtype=complex)
 for mu in range(4):
  W-=(np.kron(T[mu],np.eye(4)-g[mu])+np.kron(T[mu].conj().T,np.eye(4)+g[mu]))/2
 H=G@W;herm=np.linalg.norm(H-H.conj().T);w,V=np.linalg.eigh(H);S=(V*np.sign(w))@V.conj().T;D=np.eye(4*n)+G@S
 return dict(q=q,index=int(round(-np.sign(w).sum()/2)),gap=float(min(abs(w))),GW=float(np.linalg.norm(G@D+D@G-D@G@D)),herm=float(herm))
for L,q in [(4,0),(4,1),(4,2),(4,3),(4,4),(4,6),(5,3),(5,4)]:
 r=ov(q,L=L);r['L']=L;print(r,flush=True)
print("PASS")
