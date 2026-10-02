"""Finite gauge covariance, analytic cuts and independent field-map controls."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11280_factor_higgs_mediator as H
import w33_pass11281_family_sextet_selection as F
import w33_pass11282_all_modulus_radial_self_energy as Q
import w33_pass11283_symplectic_metric_field as G
import w33_pass11284_closed_history_flux as C

def certificate(n):return json.loads(next((ROOT/'data').glob(f'w33_pass{n}_*.json')).read_text())

def test_fundamental_factor_lift_under_both_gauge_groups():
    v,A,U,V,W=H.reference();b=H.H.setup()[0];E=expm(.07j*b[19]+.11j*b[22])
    rng=np.random.default_rng(11280);a=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));R=expm(.1j*(a+a.conj().T))
    vals=H.constraints(E@v,E@A@E.conj().T,E@U@R,E@V@R,E@W@R)
    assert max(np.max(abs(z)) for z in vals)<1e-12
    d=certificate(11280);assert d['E6_beta_combined']=='13' and d['aux_SU5_beta']=='29/6'
    assert d['real_scalar_dimensions']-d['gauge_orbit_dimension']==905
    # Common noncompact rescaling preserves factor product but is removed by Gram constraints.
    vals=H.constraints(v,A,2*U,V/2,2*W);assert np.linalg.norm(vals[0])==0 and np.linalg.norm(vals[1])>1

def test_family_selector_covariance_and_independent_normal_displacements():
    D=np.diag([1,2,3]).astype(complex);a=np.array([[.2,.1j,.15],[-.1j,-.1,.05],[.15,.05,-.1]])
    U=expm(1j*a);assert abs(np.linalg.det(U)-1)<1e-12
    assert F.potential(U@D@U.T)<1e-20
    for change in (np.diag([.001,0,0]),1j*np.diag([.001,0,0])):assert F.potential(D+change)>1e-6
    old=certificate(11281);assert old['exact_normal_rank']==4 and old['compact_SU3_orbit_dimension']==8
    # Bilinear covariant under separate E6 and family frames.
    v,A,_,_,_=H.reference();field=np.zeros((27,3,3),complex);field[0]=D
    E=expm(.04j*H.H.setup()[0][24]);v2=E@v[:,0];h2=np.einsum('ab,bij->aij',E,field)
    assert np.max(abs(np.einsum('a,aij->ij',v2,h2.conj())-D))<1e-12

def test_radial_cuts_ward_width_and_closed_thresholds():
    d=certificate(11282);M2=d['tree_radial_mass_squared'];M=np.sqrt(M2);cs=d['channels']
    gold=next(c for c in cs if c['kind']=='scalar' and c['m1_squared']==c['m2_squared']==0)
    assert abs(Q.rho(gold,M2)/M-M**3/(8*np.pi))<1e-13
    assert all(Q.rho(c,Q.threshold(c)*.9)==0 for c in cs if Q.threshold(c)>0)
    assert abs(sum(Q.rho(c,M2) for c in cs)-d['absorptive_total'])<1e-14
    assert d['strict_one_loop_pole_squared']['imag']<0 and d['width_over_mass']<.01
    # Exact soft-profile derivative reproduces Ward vertex; fermion homogeneity too.
    h=1e-6;m=d['m'];g=lambda r:m*m*81/7*(r**-8-r**-16)
    assert abs((g(1+h)-g(1-h))/(2*h*np.sqrt(2))-M2/np.sqrt(2))<1e-10
    c=next(c for c in cs if c['kind']=='Majorana' and c['m1_squared']==81*m*m)
    assert Q.rho(c,M2)==0 # heaviest Weyl pair is not kinematically accessible.

def test_metric_chart_signature_frame_covariance_and_full_rank():
    rng=np.random.default_rng(11283);D=sum(a*b for a,b in zip(rng.normal(size=6)*.15,G.tangents()));Gp=expm(D)
    u=rng.normal(size=4);u/=np.sqrt(u@Gp@u)
    g=G.metric(Gp,u,.31,.12)
    ev=np.linalg.eigvalsh(g);assert sum(ev<0)==1 and sum(ev>0)==3
    assert abs(np.linalg.det(g)+np.exp(.96))<1e-10
    K=-G.J@np.diag([.2,-.3,.1,.4]);S=expm(K);inv=np.linalg.inv(S)
    assert np.max(abs(S.T@G.J@S-G.J))<1e-12
    transformed=G.metric(inv.T@Gp@inv,S@u,.31,.12)
    assert np.max(abs(transformed-inv.T@g@inv))<1e-11
    a=G.chart_Jacobian();assert H.H.rank_mod(a,103)==10 and H.H.rank_mod(a[:,:10],103)==9 and H.H.rank_mod(a[:,:9],103)==8

def test_closed_history_rigidity_and_flux_are_not_clique_history():
    for n in (2,5,8):
        D=C.top_boundary(n);assert H.H.rank_mod(D,103)==n-1
        assert np.array_equal(D@np.ones(n,dtype=int),np.zeros(len(D),int))
        rng=np.random.default_rng(n);A=rng.normal(size=len(D));Q=7.3;F=D.T@A+Q/n
        assert abs(sum(F)-Q)<1e-12
    d=certificate(11284);assert d['top_boundary_rank_mod101']==2
    # Fixed linear-sigma sectors cancel arbitrary constant stress shifts.
    trace=np.array([1.,-2.,4.]);vol=np.array([2.,1.,3.]);shift=.7
    avg=np.dot(vol,trace)/sum(vol);avg2=np.dot(vol,trace-4*shift)/sum(vol)
    assert np.max(abs((trace/4-avg/4)-((trace-4*shift)/4-avg2/4)))<1e-14
