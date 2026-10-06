import json,sympy as sp,numpy as np,itertools
gdata=json.load(open('data/w33_pass10961_albert_clifford9_gammas.json',encoding='utf-8'))
g9=[sp.Matrix([[sp.Rational(x) for x in row] for row in M]) for M in gdata['gamma9']]
I16=sp.eye(16);Z16=sp.zeros(16)
g=[sp.Matrix.vstack(sp.Matrix.hstack(Z16,x),sp.Matrix.hstack(x,Z16)) for x in g9]+[sp.diag(I16,-I16)]
prod=sp.eye(32)
for x in g:prod=prod*x
Chi=sp.I*prod; W=sp.Matrix.hstack(*(Chi-sp.eye(32)).nullspace());p=list(W.T.rref()[1]);L=W[p,:].inv()
def rw(M):
 im=M*W; return L*im[p,:]
pairs=list(itertools.combinations(range(10),2))
R=[np.array(rw(g[i]*g[j]).evalf(),complex) for i,j in pairs]
# vector rep corresponding to B_ij=gamma_i gamma_j has 2*(Eij-Eji)
V=[]
for i,j in pairs:
 M=np.zeros((10,10),complex);M[i,j]=2;M[j,i]=-2;V.append(M)
Gram=np.array([[-np.trace(A@B).real for B in V] for A in V])
print("Gram diag/off",np.unique(np.round(Gram,10)))
Gi=np.linalg.inv(Gram)
# symmetric embedding S: C^136 -> C^(16x16)
sympairs=[(i,j) for i in range(16) for j in range(i,16)]
S=np.zeros((256,len(sympairs)),complex)
for k,(i,j) in enumerate(sympairs):
 if i==j:S[i*16+j,k]=1
 else:
  S[i*16+j,k]=1/np.sqrt(2);S[j*16+i,k]=1/np.sqrt(2)
I=np.eye(16)
Rs=[]
for A in R:
 K=np.kron(A,I)+np.kron(I,A)
 Rs.append(S.conj().T@K@S)
C=np.zeros((136,136),complex)
for a in range(45):
 for b in range(45):
  if abs(Gi[a,b])>1e-12:C-=Gi[a,b]*(Rs[a]@Rs[b])
C=(C+C.conj().T)/2
ev=np.linalg.eigvalsh(C)
# cluster
clusters=[]
for x in ev:
 if not clusters or abs(x-clusters[-1][0])>1e-7:clusters.append([x,1])
 else:clusters[-1][1]+=1
print("Sym2 Casimir clusters",clusters)
# vector Casimir
Cv=np.zeros((10,10),complex)
for a in range(45):
 for b in range(45):
  if abs(Gi[a,b])>1e-12:Cv-=Gi[a,b]*(V[a]@V[b])
print("vector Casimir",np.unique(np.round(np.linalg.eigvalsh((Cv+Cv.conj().T)/2),10)))
assert sorted(m for _,m in clusters)==[10,126]
# exterior square dimension120 and should be irreducible
apairs=list(itertools.combinations(range(16),2));Aemb=np.zeros((256,120),complex)
for k,(i,j) in enumerate(apairs):
 Aemb[i*16+j,k]=1/np.sqrt(2);Aemb[j*16+i,k]=-1/np.sqrt(2)
Ra=[]
for X in R:
 K=np.kron(X,I)+np.kron(I,X);Ra.append(Aemb.conj().T@K@Aemb)
Ca=np.zeros((120,120),complex)
for a in range(45):
 for b in range(45):
  if abs(Gi[a,b])>1e-12:Ca-=Gi[a,b]*(Ra[a]@Ra[b])
ae=np.linalg.eigvalsh((Ca+Ca.conj().T)/2)
print("Alt2 spread",float(ae.max()-ae.min()),"value",float(ae.mean()))
assert ae.max()-ae.min()<1e-7
print("PASS")
