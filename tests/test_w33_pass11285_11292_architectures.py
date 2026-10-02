"""Independent covariance, quotient vertices, topology and physical-scope controls."""
from pathlib import Path
import sys,json,math
import numpy as np
import sympy as s
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11285_semisimple_factor_higgs as H
import w33_pass11286_misaligned_family_vacua as F
import w33_pass11287_quotient_self_energy_matrix as Q
import w33_pass11288_metric_constraint_audit as G
import w33_pass11289_cycle_gram_gluing_flux as T
import w33_pass11290_gauged_family_decay as V

def cert(stem):return json.loads((ROOT/'data'/('w33_pass'+stem+'.json')).read_text())

def test_so10_covariance_and_projection_hessian_repair():
 v,A,U,W,Z,K=H.reference();b=np.zeros((10,10));b[0,6]=1;b[6,0]=-1;R=expm(.13*b)
 e=expm(.02j*H.F.H.setup()[0][0]);assert np.max(abs(e.conj().T@e-np.eye(27)))<1e-10
 vals=H.constraints(e@v,e@A@e.conj().T,e@U@R.T,e@W@R.T,e@Z@R.T,R@K@R.T)
 assert max(np.max(abs(a)) for a in vals)<1e-9
 Ep=np.vstack([np.eye(5),1j*np.eye(5)])/np.sqrt(2);du=np.outer(np.eye(27)[:,0],Ep[:,0].conj())
 assert np.max(abs(U.conj().T@du+du.conj().T@U))<1e-12
 assert np.linalg.norm(du@K+1j*du)>1
 # This is an actual normal displacement invisible to the linearized Gram:
 # explicit projection makes its positive quadratic energy visible.
 a=1e-5;vals=H.constraints(v,A,U+a*du,W,Z,K)
 assert np.linalg.norm(vals[1])/a>1

def test_mixed_moments_fix_fourier_mixing_and_cp_pair():
 p=s.symbols('p0:9');P=s.Matrix(3,3,p);u=[1,4,9];d=[1,9,25]
 eq=[sum(P[i,j] for j in range(3))-1 for i in range(3)]+[sum(P[i,j] for i in range(3))-1 for j in range(3)]
 for a in (1,2):
  for b in (1,2):eq.append(sum(u[i]**a*d[j]**b*P[i,j] for i in range(3) for j in range(3))-s.Rational(sum(x**a for x in u)*sum(x**b for x in d),3))
 solution=s.linsolve(eq,p);assert solution==s.FiniteSet((s.Rational(1,3),)*9)
 Fu,Fd=F.vacuum();assert F.potential(Fu,Fd)<1e-15 and F.potential(Fu,Fd.conj())<1e-15
 a=np.array([[0,1j,.2], [1j,0,0],[-.2,0,0]],complex);R=expm(a)
 assert F.potential(R@Fu@R.T,R@Fd@R.T)<1e-12
 omega=(-1+s.sqrt(3)*s.I)/2;j=s.im(omega/9);assert s.simplify(j*j)==s.Rational(1,108)
 assert abs(F.jarlskog(F.fourier()))>0.09

def test_full_quotient_vertices_and_basis_independent_cuts():
 d=cert('11287_quotient_self_energy_matrix');h=np.array(d['Hessian']);c=np.array(d['cubic_tensor']);m=np.array(d['Weyl_mass_real'])+1j*np.array(d['Weyl_mass_imag']);y=np.array(d['Yukawa_real'])+1j*np.array(d['Yukawa_imag'])
 assert max(abs(c-c.transpose(1,0,2)).ravel())<1e-8
 p,N,*_=Q.geometry();a=N.conj().T@p;r=np.r_[a.real,a.imag]
 z=s.symbols('z');vr=s.Rational(81,49)*z**-14-s.Rational(27,7)*z**-6
 expected=float(s.diff(vr,z,3).subs(z,1)/(2*s.sqrt(2)));actual=np.einsum('ijk,i,j,k',c,r,r,r)
 assert abs(actual-expected)<1e-7
 sm,fm,O,groups=Q.cut_groups(h,c,m,y);spectral=sum((Q.rho(g,sm[-1]) for g in groups),start=np.zeros((22,22)))
 assert min(np.linalg.eigvalsh(spectral))>-1e-11
 rng=np.random.default_rng(11287);R=np.linalg.qr(rng.normal(size=(22,22)))[0]
 h2=R.T@h@R;c2=np.einsum('ia,jb,kd,ijk->abd',R,R,R,c);y2=np.einsum('ia,ijk->ajk',R,y)
 sm2,_,_,gs=Q.cut_groups(h2,c2,m,y2);sp2=sum((Q.rho(g,sm2[-1]) for g in gs),start=np.zeros((22,22)))
 assert max(abs(np.linalg.eigvalsh(spectral)-np.linalg.eigvalsh(sp2)))<1e-10
 # The massless integration endpoint is explicit and finite.
 zero={'kind':'scalar','a':0.,'b':0.,'R':np.eye(22),'I':np.zeros((22,22))}
 assert np.isfinite(Q.dispersion(zero,.01,-.0001)).all()

def test_metric_pullback_has_fiber_and_lapse_shift_nulls():
 from w33_pass11283_symplectic_metric_field import chart_Jacobian
 A=s.Matrix(chart_Jacobian());null=A.nullspace();assert len(null)==1 and A*null[0]==s.zeros(16,1)
 spatial=s.Matrix.vstack(*[A[4*i+j,:] for i,j in [(1,1),(2,2),(3,3),(1,2),(1,3),(2,3)]])
 D=s.diag(1,1,1,2,2,2)-s.Matrix([1,1,1,0,0,0])*s.Matrix([[1,1,1,0,0,0]])
 M=spatial.T*D*spatial;assert M.rank()==6 and len(M.nullspace())==5
 ev=np.linalg.eigvalsh(np.array(M,float));assert sum(ev>1e-8)==5 and sum(ev<-1e-8)==1
 assert (22-2*9)//2==2

def test_cycle_gluing_uses_prior_critical_group_and_integral_basis():
 d=cert('11289_cycle_gram_gluing_flux');inc,C=T.cycle_basis();B=np.array(d['cycle_Gram'],dtype=np.int64)
 assert not np.any(inc@C) and np.array_equal(C.T@C,B)
 prior=json.loads((ROOT/'data/PART_W33_PASS5031_CRITICAL_GROUPS.json').read_text())['levi']['critical_group']
 own={str(x):d['Smith_diagonal'].count(x) for x in set(d['Smith_diagonal']) if x!=1};assert own==prior
 U=np.eye(81,dtype=np.int64);U[0,1]=1;C2=C@U;assert np.array_equal(C2.T@C2,U.T@B@U)
 assert all(B.diagonal()%2==0) and d['glued_product_rational_Betti']==[1,1,0,1,1]


def test_gauging_moves_goldstone_width_to_vector_thresholds():
 d=cert('11290_gauged_family_decay');assert d['rows'][0]['open_vector_channels']==8 and d['rows'][-1]['open_vector_channels']==0
 M2=.009;assert abs(V.width(1e-12,M2)/(M2**1.5/(64*np.pi))-1)<1e-8
 assert V.width(M2/4,M2)==0 and 0.99999<d['small_gF_Goldstone_width_ratio']<1.00001


def test_rank_five_window_is_inventory_conditional():
 good=[]
 for n in range(3,20):
  e6=s.Integer(34)-s.Rational(1,3)*(3*2*n)*3
  so=s.Rational(11,3)*(2*n-2)-s.Rational(1,3)*81-s.Rational(1,6)*(2*n-2)
  if e6>0 and so>0:good.append(n)
 assert good==[5]
 assert cert('11291_rank_five_AF_window')['at_r5']=={'E6':4,'SO10':1}


def test_cycle_cs_quadratic_form_and_doubled_boundary():
 d=cert('11289_cycle_gram_gluing_flux');B=s.Matrix(d['cycle_Gram']);inv=B.inv(method='DM');x=s.eye(81)[:,0];a=s.eye(81)[:,1]
 delta=((x+B*a).T*inv*(x+B*a)-x.T*inv*x)[0]/2;assert delta.is_Integer
 assert math.isqrt(d['exact_determinant'])**2!=d['exact_determinant']
 q=(x.T*inv*x)[0]/2;assert q==s.Rational(81,320)
 doubled_inverse=s.diag(inv,-inv);v=x.col_join(x);assert (v.T*doubled_inverse*v)[0]==0


def test_nonradial_cubic_against_original_d_flat_potential():
 d=cert('11287_quotient_self_energy_matrix');control=Q.finite_slice_control(np.array(d['cubic_tensor']))
 assert abs(control['Richardson_error'])<2e-5
 assert max(row['moment_residual'] for row in control['rows'])<1e-12
 assert abs(control['rows'][2]['error'])<abs(control['rows'][0]['error'])/12
