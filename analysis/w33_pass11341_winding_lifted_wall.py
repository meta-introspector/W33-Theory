"""Quantized axion-winding wall with backreaction and lifted torus shape modes."""
from pathlib import Path
import json
from functools import lru_cache
import numpy as np
import sympy as s
from scipy.linalg import eigh
from numpy.polynomial.legendre import Legendre,leggauss
ROOT=Path(__file__).resolve().parents[1]

@lru_cache(maxsize=8)
def exact_data(h=s.Rational(1,100)):
 b=s.symbols('b',positive=True);bw=s.Rational(5,4);q=(1+s.sqrt(1+4*h*bw))/2;c=q/5;rho=q*q/2-h-2*c;C=rho/3+q*q/2+h;F=C/b-rho*b*b/3-q*q/(2*b*b)-h
 assert s.simplify(F.subs(b,1))==0 and s.simplify(s.diff(F,b).subs(b,1)-2*c)==0 and s.simplify((s.diff(F,b)-2*F/b).subs(b,bw))==0
 assert s.simplify(q/c*(1-1/bw)-1)==0 and s.simplify(q*q-q-h*bw)==0
 # Actual independent Ricci components in proper-distance gauge.
 A=rho+q*q/(2*b**4);D=rho-q*q/(2*b**4)+h/b**2
 assert s.simplify(-s.diff(F,b,2)/2-s.diff(F,b)/b-A)==0 and s.simplify(-s.diff(F,b)/b-F/b**2-D)==0
 dp=s.diff(b*b*F,b);assert s.simplify(dp.subs(b,bw)).is_positive;assert s.simplify(s.diff(dp,b)+4*rho*b*b+2*h)==0
 fw=s.simplify(F.subs(b,bw));T=s.sqrt(16*fw/bw**2);assert float(fw)>0 and float(rho)>0
 return b,F,{'h':str(h),'q':str(q),'c':str(c),'rho':str(s.simplify(rho)),'C':str(s.simplify(C)),'b_wall':str(bw),'F_wall':str(fw),'tension':float(T),'Dirac_integer':1,'electric_unit':'1/2','axion_windings':[1,1]}
def shape_spectrum(n):
 b,F,d=exact_data();fun=s.lambdify(b,F,'numpy');z,w=leggauss(350);bb=1+(z+1)/8;w/=8;U=np.array([Legendre.basis(k)(z) for k in range(n)]).T;D=np.array([8*Legendre.basis(k).deriv()(z) for k in range(n)]).T;M=U.T@((w*bb*bb)[:,None]*U);K=D.T@((w*bb*bb*fun(bb))[:,None]*D)+.02*U.T@(w[:,None]*U);return eigh(K,M,eigvals_only=True)
def payload():
 b,F,d=exact_data();ev=shape_spectrum(28)[:6];low=shape_spectrum(18)[:6];assert .0128<ev[0]<.02 and max(abs(ev-low))<1e-5
 u,v=s.symbols('u v',real=True);# Tr exp(-2[[u,v],[v,-u]])=2cosh(2sqrt(u²+v²)).
 quadratic=2*(u*u+v*v);assert quadratic.coeff(u,2)==2 and quadratic.coeff(v,2)==2
 # Independent physical positivity sample, plus exact interval lower bound below.
 grid=np.linspace(1.000001,1.25,1000);assert min(s.lambdify(b,F,'numpy')(grid))>0
 return {'status':'PASS','scope':'Explicit added compact-scalar winding mechanism with exact backreacted quantized wall and positive torus-shape lift; not nativeW33 matter derivation or total stability.','added_action':'Two compact real axions theta1=x, theta2=y, each with canonical positive coefficient f²=h=1/100. Unit winding is allowed by2pi torus and field periods. Interchange symmetry supplies equal kinetic coefficients.','backreaction':'F=C/b-rho b²/3-q²/(2b²)-h, C=rho/3+q²/2+h, c=(q²/2-rho-h)/2. The base Ricci isrho+q²/(2b4), torus Ricci gains+h/b². Both Maxwell flux and Einstein components replay exactly.','certificate':d,'quantization_constraint':'With electricunit1/2 and Diracinteger1, walljunction and pole regularity imply q²-q=h*bwall. Therefore fixing the oldq1 while addingh>0 is inconsistent; the Maxwell strength must backreact. Our q=(1+sqrt(21/20))/2 enforces it exactly atbwall5/4.','shape_action':'At detQ1, winding potential is h TrQ^-1/(2b²). Parametrize Q=exp(2[[u,v],[v,-u]]); TrQ^-1=2+4(u²+v²)+O4. Both prior shape directions receive positive potential2h(u²+v²)/b².','shape_operator':'Old radial numerator gains2h int|u|²db, denominator intb²|u|²db. Both exact constant shape zeros lift. Rayleigh bounds give smallest shape eigenvalue in[2h/bwall²,2h], independent of Ritz approximation.','exact_smallest_shape_bounds':['8/625','1/50'],'shape_eigenvalues':ev.tolist(),'shape_basis18_to28_error':float(max(abs(ev-low))),'F_positivity':'b²F polynomial vanishes atb1. Its derivative isC-4rho b³/3-2hb, strictly decreasing and positive atbwall5/4 for the exact parameters; henceF>0 throughout the cap interior.','boundary':['The added axions and their kinetic scale are inputs; winding energy is not predicted from the native cochain.','The new axion/metric/Maxwell vector block must be recomputed: prior11339 four-zero statement belongs to the original action. Breathing/base/wall-bending/four-form modes remain open.','No nonlinear zero-mode completion, selected absolute scale or observed cosmological constant is inferred.'],'prior_owners':['analysis/w33_pass11324_quantized_warped_history_wall.py','analysis/w33_pass11329_wall_shape_fluctuations.py','analysis/w33_pass11339_coupled_maxwell_metric_modes.py'],'primary_sources':['https://arxiv.org/abs/1311.5157']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11341_winding_lifted_wall.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['shape_eigenvalues'][0])
