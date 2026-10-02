"""Exact positive-tension Maxwell-supported closed one-handle history wall."""
from pathlib import Path
import json,sys
import sympy as s
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))

def exact_data():
 b=s.symbols('b',positive=True);rho=s.Rational(1,10);q=s.S(1);k=s.S(1);b0=s.S(1);c=s.Rational(1,5);C=s.Rational(8,15)
 F=C/b-rho*b*b/(3*k)-q*q/(2*k*b*b);bw=s.Rational(5,4);Fw=s.factor(F.subs(b,bw));T=4*k*s.sqrt(Fw)/bw;rwall=s.sqrt(Fw)/c
 assert F.subs(b,b0)==0 and s.diff(F,b).subs(b,b0)==2*c
 assert s.factor(s.diff(F,b)-2*F/b).subs(b,bw)==0
 # Independent Einstein equations in proper-distance gauge, b'=sqrtF, r=sqrtF/c.
 bp=s.sqrt(F);r=s.sqrt(F)/c;rp=s.diff(F,b)/(2*c);rpp=s.diff(rp,b)*bp;bpp=s.diff(F,b)/2
 A=(rho+q*q/(2*b**4))/k;D=(rho-q*q/(2*b**4))/k
 assert s.simplify(-rpp/r-2*bpp/b-A)==0
 assert s.simplify(-rpp/r-2*rp/r*bp/b-A)==0
 assert s.simplify(-bpp/b-bp*bp/b**2-rp/r*bp/b-D)==0
 flux=s.factor(4*s.pi*q/c*(1/b0-1/bw));vol=s.factor(2*(2*s.pi)**3*(bw**3-b0**3)/(3*c));area=(2*s.pi)**3*rwall*bw*bw;average=s.factor(4*rho/k+3*T*area/(k*vol))
 e=s.Rational(1,2);assert e*flux/(2*s.pi)==1
 # The naive one-junction radial potential is NOT a full admissible mode.
 radial=s.factor(s.diff(F-T*T*b*b/(16*k*k),b,2).subs(b,bw));mismatch=s.factor(s.diff(s.diff(F,b)/(2*F)-1/b,b).subs(b,bw));assert radial<0 and mismatch!=0
 return {'b0':str(b0),'c':str(c),'C':str(C),'rho':str(rho),'q_Maxwell':str(q),'kappa_squared':str(k),'b_wall':str(bw),'F_wall':str(Fw),'r_wall_squared':str(s.factor(rwall*rwall)),'tension_squared':str(s.factor(T*T)),'tension':float(T),'magnetic_flux':str(flux),'electric_charge_unit':str(e),'Dirac_integer':1,'four_volume':str(vol),'average_scalar_curvature_with_wall':str(average),'wall_area_times_tension':str(s.factor(area*T)),'naive_radial_potential_second_derivative':str(radial),'other_spatial_junction_linear_mismatch':str(mismatch)}

def numerical_replay(h):
 c=.2;rho=.1;q=1.
 def flow(x,v):
  r,rp,b,I,J=v;bp=c*r;A=rho+q*q/(2*b**4)
  return [rp,-A*r-2*rp*bp/b,bp,r*b*b,q*r/b**2]
 a0=rho+q*q/2+4*c;init=[h-a0*h**3/6,1-a0*h*h/2,1+c*h*h/2,h*h/2,q*h*h/2]
 def end(x,v):return v[2]-1.25
 end.terminal=True;z=solve_ivp(flow,[h,3],init,events=end,rtol=2e-12,atol=2e-13,max_step=.01);assert z.success and len(z.t_events[0])==1
 r,rp,b,I,J=z.y_events[0][0];junction=abs(rp/r-c*r/b);constraint=abs(2*rp/r*c*r/b+(c*r/b)**2+rho-q*q/(2*b**4));return {'pole_offset':h,'proper_cap_length':float(z.t_events[0][0]),'junction_error':float(junction),'bulk_constraint_error':float(constraint),'volume_error':float(abs(2*(2*np.pi)**3*I-305*np.pi**3/12)),'flux_error':float(abs(4*np.pi*J-4*np.pi))}

def payload():
 from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
 _,cycles=cycle_basis();chord=next(i for i,row in enumerate(cycles) if row[0]==1 and not np.any(row[1:]));cochain=np.zeros(160,dtype=int);cochain[chord]=1;assert np.array_equal(cochain@cycles,np.r_[1,np.zeros(80,dtype=int)])
 d=exact_data();controls=[numerical_replay(h) for h in [1e-5,1e-6]];assert max(x['junction_error'] for x in controls)<1e-9
 return {'status':'PASS','result_scope':'PASS_EXACT_QUANTIZED_MAXWELL_SUPPORTED_CLOSED_HISTORY_WALL','geometry':'Two identical regular caps, ds²=dxi²+r(xi)²dphi²+b(xi)²(dx²+dy²), three angles have period2pi. Collapse onlyphi at each pole. Double at a pure-tension wall with b=5/4. TopologyS2 x T2=(S2 x S1) x S1.','bulk_solution':'bprime=c r, (bprime)²=F(b)=C/b-rho b²/(3kappa²)-q²/(2kappa² b²), r=sqrtF/c. F>0 for1<b<=5/4. Positive Maxwell action, F2=(q/b²) orthonormal_volume(S2).','native_one_handle_quotient':{'cochain':cochain.tolist(),'cycle_periods':(cochain@cycles).tolist(),'scope':'Explicit primitive H1 projection selects one native cycle and kills the other80. Identifying its circle with the S1 factor is an added quotient choice, not equivalence to the81-handle history.'},'exact_certificate':d,'numerical_controls':controls,'wall_junction':'Reflection doubling gives Kplus-Kminus=-T/(2kappa²)gamma. Both independent eigenvalues match because Fprime(bwall)=2F(bwall)/bwall. The worldvolume mean curvature and mean four-form flux vanish, satisfying the charged wall force equation.','charged_membrane_sector':'Take Lorentzian fplus=1/5,fminus=-1/5, membrane charge2/5 and common bare vacuum2/25. Totalrho=bare+f²/2=1/10 on both sides. Electric four-form uses the usual Wick-continued Euclidean saddle; Maxwell2-form is real magnetic. Totalmembrane4-flux cancels by reflection. Nontrivial Dirac integer1 is onS2 magneticflux, not a claimed nonzero4-form sector.','Lorentzian_continuation':'Continue one torus coordinate to time: ds²=dxi²+r²dphi²+b²(dy²-dt²). F2 remains real magnetic on the spatialS2, and the wall is timelike. Canonical Maxwell and positive wall tension satisfy the Lorentzian null-energy condition; no Ellis phantom throat is introduced.','sequestering_extension':'One may set the distinct global sequestering inputs Qseq=mu4*(305pi³/12), Qhat=-(515pi³/24) for M0²=1. These are chosen sectors, not a prediction of the cosmological constant.','fluctuation_boundary':'The naive one-junction radial potential has negative curvature-512/625, but a displacement with frozen bulk violates the second spatial junction. It is not an established negative eigenmode. Full coupled bulk/wall perturbations and determinant remain open.','scope':['Explicit exact closed one-handle history wall with quantized magnetic flux and positive sources, not the full81-handle history or an emergent topology.','Tension, bare vacuum, gauge coupling, Newton scale, torus periods and history choice are inputs; no observedCC or mass prediction.','Static saddle existence is established; full fluctuation spectrum, nucleation rate and state selection are not.'],'prior_owners':['analysis/w33_pass11317_history_einstein_obstruction.py','analysis/w33_pass11319_smooth_neck_stress.py','analysis/w33_pass11302_sequestered_membrane_junction.py'],'primary_sources':['https://arxiv.org/abs/1604.04000','https://arxiv.org/html/1505.01492v2']}
if __name__=='__main__':
 d=payload();(ROOT/'data/w33_pass11324_quantized_warped_history_wall.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['exact_certificate'])
