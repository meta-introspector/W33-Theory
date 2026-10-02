"""A reduced fixed-flux saddle, quantized membrane period and constrained radial mode."""
from pathlib import Path
import sys,json
import numpy as np
import mpmath as mp
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11302_sequestered_membrane_junction as B
from w33_pass11284_closed_history_flux import top_boundary
mp.mp.dps=60

def action(x,Qmem,Q,Qhat,mu,bare=0):
    R,L,G,f=x;q=mp.mpf('.4');T=mp.mpf('.4');area=2*mp.pi**2*R**3
    out=T*area-L*Q/mu-G*Qhat-f*Qmem
    for sign,ff in [(-1,f+q),(1,f)]:
        rho=L+bare+ff*ff/2;h=rho/(3*G);c=sign*mp.sqrt(1-h*R*R)
        vol=2*mp.pi**2/h**2*(mp.mpf(2)/3-c+c**3/3)
        out-=rho*vol+3*G*c*area/R
    return out

def gradient(x,Qmem,Q,Qhat,mu,bare=0):
    vals=[];x=list(x)
    for i in range(4):
        def f(t):
            z=x.copy();z[i]=t;return action(z,Qmem,Q,Qhat,mu,bare)
        vals.append(mp.diff(f,x[i]))
    return vals

def saddle():
    r,l,g,mu,Q,Qhat=B.inputs();vol,_,_=B.caps(r,l,g);Qold=float(1.2*vol[0]+.8*vol[1]);N=round(.4*Qold/(2*np.pi));Qmem=2*mp.pi*N/mp.mpf('.4')
    scale=np.array([1000.,Q/mu,abs(Qhat),float(Qmem)])
    def f(x):
        try:return np.array([float(v) for v in gradient(list(map(mp.mpf,x)),Qmem,Q,Qhat,mu)])/scale
        except (ValueError,TypeError,ZeroDivisionError):return np.ones(4)*1e5
    sol=root(f,[r,l,g,.8],tol=1e-10);assert sol.success and max(abs(f(sol.x)))<1e-9
    x=mp.findroot(lambda R,L,G,f:tuple(gradient([R,L,G,f],Qmem,Q,Qhat,mu)),tuple(sol.x),tol=mp.mpf('1e-45'),maxsteps=30)
    return list(x),N,Qmem,mp.mpf(Q),mp.mpf(Qhat),mp.mpf(mu),Qold

def hessian(x,Qmem,Q,Qhat,mu):
    H=mp.matrix(4)
    for i in range(4):
        for j in range(4):
            def f(t,u):
                z=x.copy();z[i]=t;z[j]=u;return action(z,Qmem,Q,Qhat,mu)
            if i==j:
                def a(t):
                    z=x.copy();z[i]=t;return action(z,Qmem,Q,Qhat,mu)
                H[i,j]=mp.diff(a,x[i],2)
            else:H[i,j]=mp.diff(f,(x[i],x[j]),(1,1))
    return H

def payload():
    x,N,Qmem,Q,Qhat,mu,Qold=saddle();H=hessian(x,Qmem,Q,Qhat,mu);aux=H[1:4,1:4]
    curvature=H[0,0]-(H[0,1:4]*aux**-1*H[1:4,0])[0];assert curvature<0
    residual=max(abs(v) for v in gradient(x,Qmem,Q,Qhat,mu));assert residual<mp.mpf('1e-35')
    shifted=x.copy();shifted[1]-=mp.mpf('.03');shift_error=max(abs(v) for v in gradient(shifted,Qmem,Q,Qhat,mu,bare=mp.mpf('.03')));assert shift_error<mp.mpf('1e-35')
    D=top_boundary();assert not np.any(D@np.ones(3,dtype=int));flux=np.full(3,float(Qmem)/3);assert abs(.4*sum(flux)-2*np.pi*N)<1e-10
    return {'status':'PASS','result_scope':'PASS_QUANTIZED_FIXED_FLUX_REDUCED_SADDLE_WITH_CONSTRAINED_RADIAL_NEGATIVE_MODE',
      'reduced_action':'S=sum_caps[-rho V-3kappa² c Area/R]+T Area-Lambda Qseq/mu4-kappa² Qhat-fminus Qmem. The cap terms include GHY; fplus=fminus+q. This declares a fixed-flux ensemble. It is a reducedO4 action, not the complete field fluctuation determinant.',
      'stationarity':'dR S gives the Israel junction; dLambda S=Vol-Qseq/mu4; dkappa² S fixes bulk plus geometric-junction integrated curvature; dfminus S=fplus Vplus+fminus Vminus-Qmem. Away from the Israel saddle, the geometric junction6*c_sum*Area/R must replace the on-shell wall3T*Area/kappa².',
      'membrane_flux_quantization':{'compact_form_assumption':'q integral Fmem=2piN','charge':.4,'integer_sector':N,'unquantized_reference_period':Qold,'quantized_period':float(Qmem)},
      'saddle':{'radius':float(x[0]),'Lambda':float(x[1]),'kappa_squared':float(x[2]),'fminus':float(x[3])},'gradient_residual':float(residual),
      'full_reduced_Hessian':[[float(v) for v in row] for row in H.tolist()],'fixed_auxiliary_radial_curvature':float(H[0,0]),'fixed_flux_constrained_radial_curvature':float(curvature),
      'negative_mode_scope':'EliminateLambda,kappa²,fminus by their stationarity equations at fixedQseq,Qhat,Qmem. The Schur complement is the second derivative of the reduced radius action and is negative. This is oneO4 collective-coordinate negative direction; conformal-factor contours, nonradial modes and a full bounce count/rate are not proved.',
      'vacuum_shift_gradient_error':float(shift_error),
      'history_topology_map':'For the prior oriented closedM4=#81(S1xS2) x S1, collapse the complement of an embedded oriented4-ball to a point. M4->B4/boundaryB4=S4 has degree1 and pulls the integralH4 generator back to the history generator. On the existing3-tick chain, Fmem=(Qmem/3)(1,1,1)+D4^T A has totalperiodQmem. This transfers a flux class, not the two-cap metric or instanton.',
      'boundaries':['Compact membrane gauge group, q,T and the linear sequestering flux sectors are additional inputs.','No observed CC selection, nucleation rate or complete genus81 instanton geometry.','The degree1 collapse is a topological map; pulling back a spherical cap metric would be degenerate and is not performed.'],
      'prior_owners':['analysis/w33_pass11302_sequestered_membrane_junction.py','analysis/w33_pass11284_closed_history_flux.py'],
      'primary_sources':['https://arxiv.org/abs/1604.04000','https://arxiv.org/abs/1210.4740']}
if __name__=='__main__':
    out=payload();(ROOT/'data/w33_pass11307_quantized_membrane_radial_mode.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'],out['fixed_flux_constrained_radial_curvature'])
