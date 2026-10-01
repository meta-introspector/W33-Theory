#!/usr/bin/env python3
"""Curved Dirac heat-response benchmark separating gravity from finite mass-volume terms.

Prior: BT1130 product heat bookkeeping; BT1033 geometric spectral-action route.
Here the external torus is input. This does not reconstruct spacetime from W33.
"""
from pathlib import Path
import json,argparse
import numpy as np
import sympy as s
from scipy.special import iv
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_20261001_isovolume_dirac_gravity_response.json'

def shells(K):
    x=np.arange(-K,K+1,dtype=int)
    a,b,c=np.meshgrid(x,x,x,indexing='ij');counts=np.bincount((a*a+b*b+c*c).ravel())
    q=np.flatnonzero(counts)
    return q.astype(float),counts[q].astype(float)

def heat_response(t,K=42):
    q,w=shells(K);n=np.arange(-K,K+1,dtype=float)[:,None]
    E=n*n+q[None,:];Ep=(n+1)**2+q[None,:];d=2*n+1
    ee=np.exp(-t*E);dd=ee*(-np.expm1(-t*d))/d
    scalar=E+Ep-.5
    # Exact second-order perturbation structure; finite unbounded sums are cut off.
    term=-4*t*(3*E+.125)*ee+t*(scalar*scalar+q[None,:])*dd
    return float(np.sum(term*w[None,:]))

def block_control(t=.2,K=14,q=4.,eps=.001):
    n=np.arange(-K,K+1,dtype=float);nm=(n[:,None]+n[None,:])/2
    diff=np.abs(n[:,None]-n[None,:]).astype(int)
    def trace(e):
        a=iv(0,4*e)**.25 * iv(diff,-e)
        L=nm*a;R=np.sqrt(q)*a
        D=np.block([[L,R],[R,-L]])
        lam=np.linalg.eigvalsh(D)
        return float(2*np.sum(np.exp(-t*lam*lam))) # four spin components
    h0=trace(0)
    finite=(trace(eps)+trace(-eps)-2*h0)/(2*eps**2)
    half=(trace(eps/2)+trace(-eps/2)-2*h0)/(2*(eps/2)**2)
    numerical=(4*half-finite)/3
    E=n*n+q;diag=2.5*E
    for shift in (-1,1):
        mask=(n+shift>=-K)&(n+shift<=K)
        diag+=mask*((n+shift/2)**2+q)/4
    exact_structure=-4*t*np.sum(diag*np.exp(-t*E))
    m=n[:-1];Ea=m*m+q;Eb=(m+1)**2+q;d=2*m+1
    dd=np.exp(-t*Ea)*(-np.expm1(-t*d))/d
    exact_structure+=t*np.sum(((Ea+Eb-.5)**2+q)*dd)
    err=abs(numerical-exact_structure)
    assert err<1e-5
    return {'finite_Galerkin_modes':2*K+1,'transverse_squared_momentum':q,
      'heat_time':t,'finite_epsilon':eps,'Richardson_response':float(numerical),
      'perturbation_response':float(exact_structure),'absolute_error':float(err)}

def gaussian_coefficients():
    u,t=s.symbols('u t');a=u*(1-u);c=s.Rational(1,2)-2*a
    Q=24+t*(s.Rational(15,2)-24*a)+t*t*c*c
    expr=s.expand(Q*s.series(s.exp(-t*a),t,0,7).removeO())
    out=[s.integrate(expr.coeff(t,k),(u,0,1)) for k in range(7)]
    assert out==[24,-s.Rational(1,2),0,s.Rational(1,140),-s.Rational(1,1260),s.Rational(1,18480),-s.Rational(1,360360)]
    return list(map(str,out))

def payload():
    rows=[]
    spectra={'W33_vertex_Laplacian':{0:1,10:24,16:15},
             'separate_legacy_480_mode_fixture':{0:82,4:320,10:48,16:30}}
    for t in [.2,.1,.05,.025]:
        v=heat_response(t,42);v2=heat_response(t,46)
        ratio=-t*v/np.pi**2
        rows.append({'t':t,'epsilon_squared_heat_response':v,
          'normalized_curvature_coefficient':float(ratio),
          'cutoff_42_46_difference':abs(v2-v),
          'finite_factors':{name:float(sum(mult*np.exp(-t*m2) for m2,mult in spec.items())) for name,spec in spectra.items()}})
    assert max(x['cutoff_42_46_difference'] for x in rows)<1e-7
    assert abs(rows[-1]['normalized_curvature_coefficient']-1)<5e-6
    return {
      'schema':'w33.20261001.isovolume-dirac-gravity-response.v1',
      'status':'PASS_CURVED_DIRAC_RESPONSE_WITH_SEPARATE_FINITE_MASS_FACTOR',
      'geometry':{'input':'external four-torus of periods 2*pi, periodic spin structure',
        'metric':'g_epsilon=exp(2 sigma) g_flat',
        'sigma':'epsilon*cos(x1) - (1/4)*log(I0(4*epsilon))',
        'volume':'exactly (2*pi)^4 for every epsilon',
        'integrated_scalar_curvature':'3*(2*pi)^4*epsilon^2 + O(epsilon^4)',
        'Dirac_on_flat_L2':'Dtilde=exp(-sigma/2) D0 exp(-sigma/2)',
        'Fourier_matrix':'D_nm=a_(n-m)[gamma1*(n+m)/2 + gamma_perp dot p_perp], a=exp(-sigma)'},
      'analytic':{'coefficient':'[epsilon^2] Tr exp(-t Dtilde^2) = -pi^2/t + pi^2*t/140 - pi^2*t^2/1260 + O(t^3), plus exponentially small Poisson images',
        'Gaussian_integral_coefficients':gaussian_coefficients(),
        'a0_response':0,'a4_response':0,
        'interpretation':'negative R coefficient is the standard Euclidean Einstein-Hilbert convention with positive inverse Newton coupling',
        'no_mass_volume_conflation':'C2(g)=N*A2(g)-F2*A0(g). The -F2*A0 term multiplies volume, not scalar curvature. Exact isovolume deformations remove it from the curvature response.'},
      'numerical_curvature_replay':rows,'independent_finite_epsilon_control':block_control(),
      'finite_product':{'formula':'Theta_product(t,g)=Theta_M(t,g)*sum_j multiplicity_j*exp(-t*m_j^2)',
        'spectra_are_distinct_declared_fixtures':spectra,
        'positivity':'each positive-multiplicity Hermitian mass spectrum gives a strictly positive finite heat factor; it cannot reverse the conventional EH sign',
        'physical_scale':'dimensionless mass-squared fixtures; conversion to measured units is an additional input'},
      'boundary':[
        'A numerical curved-continuum response benchmark with an analytic perturbative coefficient, not interval-certified lattice-sum error bounds.',
        'The smooth four-manifold and metric family are supplied, not emergent from W33.',
        'No discrete curved refinement convergence theorem, Lorentzian dynamics or observed Newton constant is proved.',
        'Finite spectra are model inputs, not measured particle masses; the 40 and 480 carriers are not conflated.',
        'The Euclidean conformal-factor instability and quantum gravitational measure are not resolved.'
      ],
      'prior_owners':['analysis/BT1130_ricci_flat_seed_paradox_resolution.md','analysis/BT1033_spectral_action_term_by_term_geometric.md'],
      'external_sources':['https://arxiv.org/abs/hep-th/9606001','https://arxiv.org/abs/1409.4983'],
      'checks':{'exact_isovolume_family':True,'analytic_curvature_coefficient':True,'independent_block_control':True,'two_cutoff_control':True,'distinct_finite_carriers':True}
    }

def compare_certificate(a,b):
    """Keep analytic fields exact; allow BLAS variation in numerical controls."""
    if isinstance(a,float) or isinstance(b,float):
        return abs(float(a)-float(b)) <= 1e-7+1e-10*abs(float(b))
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(compare_certificate(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list):
        return len(a)==len(b) and all(compare_certificate(x,y) for x,y in zip(a,b))
    return a==b

def main():
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+'\n'
    if args.check:assert compare_certificate(json.loads(txt),json.loads(OUT.read_text()))
    else:OUT.write_text(txt)
    print(p['status'])
if __name__=='__main__':main()
