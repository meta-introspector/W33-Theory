#!/usr/bin/env python3
"""Nearest-neighbor curved Wilson tower with the actual signed W33 mass fiber.

The torus, metric and dimension remain inputs. Fixed-t refinement is a numerical
control, not a certified continuum limit or emergent spacetime theorem.
Prior: Wilson controls Pass4057/4169, BT1033 geometric route, dated Dirac benchmark.
"""
from pathlib import Path
import sys,json,argparse,itertools
import numpy as np
from scipy.special import iv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_exact_cartan_mirror_completion as C
from w33_20261001_isovolume_dirac_gravity_response import heat_response,compare_certificate
OUT=ROOT/'data/w33_20261001_wilson_gravity_refinement.json'

def transverse_shells(N,wilson):
    assert N%2==0 and N>=6
    h=2*np.pi/N;p=2*np.pi*np.fft.fftfreq(N);si=np.sin(p)/h;wi=wilson*(1-np.cos(p))/h
    top=N//4 if wilson==0 and N%4==0 else N//2
    triples=np.array(list(itertools.combinations_with_replacement(range(top+1),3)))
    if top==N//4:
        signs=np.prod(np.where((triples==0)|(triples==top),2,4),axis=1)
    else:
        signs=np.prod(np.where((triples==0)|(triples==N//2),1,2),axis=1)
    perms=np.where(triples[:,0]==triples[:,2],1,np.where((triples[:,0]==triples[:,1])|(triples[:,1]==triples[:,2]),3,6))
    multiplicities=signs*perms;assert sum(multiplicities)==N**3
    return si,wi,np.sum(si[triples]**2,axis=1),np.sum(wi[triples],axis=1),multiplicities

def block_perturbation(N,t,qs,wp,weights,wilson=1):
    h=2*np.pi/N;p=2*np.pi*np.fft.fftfreq(N);si=np.sin(p)/h;wi=wilson*(1-np.cos(p))/h
    out=0.
    for n in range(N):
        a=si[n];W=wi[n]+wp;E=a*a+qs+W*W;H2=19*E/8
        for m in ((n-1)%N,(n+1)%N):
            Wm=wi[m]+wp;Em=si[m]**2+qs+Wm*Wm;dot=a*si[m]+qs+W*Wm
            H2=H2+Em/16+dot/4
        out+=float(np.sum(weights*(-4*t*H2*np.exp(-t*E))))
        m=(n+1)%N;Wm=wi[m]+wp;Em=si[m]**2+qs+Wm*Wm;dot=a*si[m]+qs+W*Wm
        H1norm=((E+Em+2*dot)**2+4*(E*Em-dot*dot))/4
        d=Em-E;dd=np.zeros_like(E);mask=abs(d)>1e-9
        dd[mask]=(np.exp(-t*E[mask])-np.exp(-t*Em[mask]))/d[mask]
        dd[~mask]=t*np.exp(-t*(E[~mask]+Em[~mask])/2)
        out+=float(np.sum(weights*t*H1norm*dd))
    return out

def response(N,t,wilson=1):
    _,_,q,w,multiplicities=transverse_shells(N,wilson)
    return block_perturbation(N,t,q,w,multiplicities,wilson)

def gammas():
    I=np.eye(2);X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1])
    gs=[np.kron(X,I),np.kron(Y,I),np.kron(Z,X),np.kron(Z,Y),np.kron(Z,Z)]
    for i,a in enumerate(gs):
        for j,b in enumerate(gs):assert np.max(abs(a@b+b@a-2*int(i==j)*np.eye(4)))==0
    return gs

def independent_control(N=12,t=.2,eps=.001,indices=(0,1,2)):
    gs=gammas();h=2*np.pi/N;p=2*np.pi*np.fft.fftfreq(N);si=np.sin(p)/h;wi=(1-np.cos(p))/h
    base=np.array([sum(si[k]*gs[j+1] for j,k in enumerate(indices))+si[n]*gs[0]+(wi[n]+sum(wi[k] for k in indices))*gs[4] for n in range(N)])
    D0=np.zeros((4*N,4*N),complex)
    for n in range(N):D0[4*n:4*n+4,4*n:4*n+4]=base[n]
    x=2*np.pi*np.arange(N)/N
    def heat(e):
        sigma=e*np.cos(x)-np.log(np.mean(np.exp(4*e*np.cos(x))))/4
        assert abs(np.mean(np.exp(4*sigma))-1)<1e-12
        b=np.fft.fft(np.exp(-sigma/2))/N
        B=np.array([[b[(n-m)%N] for m in range(N)] for n in range(N)])
        BB=np.kron(B,np.eye(4));D=BB@D0@BB
        assert np.max(abs(D-D.conj().T))<1e-12
        return float(np.sum(np.exp(-t*np.linalg.eigvalsh(D)**2)))
    h0=heat(0)
    finite=lambda e:(heat(e)+heat(-e)-2*h0)/(2*e*e)
    extrap=(4*finite(eps/2)-finite(eps))/3
    q=np.array([sum(si[k]**2 for k in indices)]);w=np.array([sum(wi[k] for k in indices)])
    exact=block_perturbation(N,t,q,w,np.ones(1))
    assert abs(extrap-exact)<2e-5
    return {'N':N,'transverse_indices':list(indices),'finite_epsilon':eps,
            'Richardson_response':extrap,'perturbation_response':exact,'error':abs(extrap-exact)}

def finite_fiber(t):
    _,op=C.operators();Q=C.plane()/np.sqrt(3);z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3)
    A=op(Q@q)+.001*np.outer((Q@q).conj(),(Q@q).conj())
    DF=np.block([[np.zeros_like(A),A],[A.conj().T,np.zeros_like(A)]])
    grading=np.diag(np.r_[np.ones(81),-np.ones(81)])
    assert np.max(abs(grading@DF+DF@grading))==0
    m2=np.linalg.svd(A,compute_uv=False)**2
    theta=float(2*np.sum(np.exp(-t*m2)))
    return {'finite_Hilbert_dimension':162,'heat_factor':theta,
      'grading':'Gamma_F=diag(I81,-I81), {Gamma_F,D_F}=0 exactly by block construction',
      'product_Dirac':'D_total=D_Wilson(g) tensor Gamma_F + I tensor D_F',
      'product_square':'D_total^2=D_Wilson(g)^2 tensor I + I tensor D_F^2',
      'point_of_internal_grading':'The Wilson regulator breaks external chiral anticommutation; using the internal grading preserves the product heat identity without assuming it.',
      'scope':'The actual signed81 mass map plus the prior one-null-mode lift; not the unrelated40/480 fixtures. The extra two gradient lifts are a separate EFT option.'}

def payload(levels=(48,96,192,384),t=.2):
    continuum=heat_response(t,46);rows=[]
    for N in levels:
        val=response(N,t);rows.append({'N':N,'spacing':2*np.pi/N,'response':val,
                                      'curvature_ratio':-t*val/np.pi**2,'continuum_error':abs(val-continuum)})
    assert all(rows[i+1]['continuum_error']<rows[i]['continuum_error'] for i in range(len(rows)-1))
    rich=(4*rows[-1]['response']-rows[-2]['response'])/3
    error=abs(rich-continuum);assert error<.05
    naive_rows=[{'N':N,'response':response(N,t,wilson=0)} for N in (96,192,384)]
    for row in naive_rows:row['ratio_to_single_continuum_species']=row['response']/continuum
    naive_extrap=(4*naive_rows[-1]['response']-naive_rows[-2]['response'])/3
    naive_ratio=naive_extrap/continuum
    assert abs(naive_ratio-16)<.2
    fiber=finite_fiber(t)
    return {'schema':'w33.20261001.wilson-gravity-refinement.v1','status':'PASS_DISCRETE_CURVED_REFINEMENT_CONTROL_WITH_ACTUAL_W33_MASS_FIBER',
      'geometry':{'input':'four-torus periods2pi; conformal metric varying along x1',
        'finite_isovolume':'sigma_j=epsilon cos(2pi*j/N)-(1/4)log(mean exp(4epsilon cos)); h^4 sum exp(4sigma)=(2pi)^4',
        'Wilson_operator':'D0(p)=sum_mu gamma_mu sin(hp_mu)/h + gamma5 sum_mu[1-cos(hp_mu)]/h',
        'curved_operator':'D(g)=b D0 b, b=exp(-sigma/2); Hermitian nearest-neighbor stencil',
        'refinement':'N -> 2N at fixed periods, t and mass map; no unscaled graph Laplacian substitution'},
      'perturbation':{'H2_diagonal_scalar':'19E_n/8 + sum_neighbors[E_m/16 + v_n dot v_m/4]',
        'Tr_H1nm_H1mn':'[(E_n+E_m+2 v_n dot v_m)^2+4(E_n E_m-(v_n dot v_m)^2)]/4',
        'divided_difference':'(exp(-tE_n)-exp(-tE_m))/(E_m-E_n), equal-energy limit t exp(-tE_n)',
        'transverse_compression':'exact signed-permutation multiplicities sum to N^3'},
      'fixed_heat_time':t,'continuum_response':continuum,'refinement_rows':rows,
      'Richardson_response':rich,'Richardson_absolute_error':error,
      'independent_finite_matrix_control':independent_control(),
      'doubling_control':{'Wilson0_refinement_rows':naive_rows,'Richardson_response':naive_extrap,
        'Richardson_ratio_to_single_continuum_species':naive_ratio,'expected_species':16,
        'coarse_grid_warning':'N96 gives a wrong-sign curvature coefficient; the single-grid doubling tolerance was rejected. Only the finer sequence and its extrapolation are used as evidence.'},
      'W33_mass_fiber':fiber,'product_curvature_response':fiber['heat_factor']*rich,
      'boundary':['Finite numerical refinement evidence, not interval tail bounds or a convergence theorem.','Four dimensions, torus topology, metric and spacing are inputs; spacetime emergence and Lorentzian gravity remain open.','Wilson heavy modes cause substantial finite-spacing artifacts; extrapolation is a control, not a proof.','No measured Newton constant, cosmological constant or particle assignments are inferred.'],
      'prior_owners':['analysis/w33_20261001_isovolume_dirac_gravity_response.py','analysis/w33_20261001_dflat_cartan_mass_loop.py','analysis/w33_20261001_exact_cartan_mirror_completion.py','analysis/BT4169_BT4176_discrete_c2_hawking_backreaction_gray_levi_casimir_axion.md','analysis/BT1033_spectral_action_term_by_term_geometric.md'],
      'external_sources':['https://arxiv.org/abs/hep-th/9606001','https://arxiv.org/abs/1409.4983']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();out=payload()
    if args.check:assert compare_certificate(out,json.loads(OUT.read_text()))
    else:OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])
