"""Constructed continuations of11438--11442; no claim of completed TOE.

Native tensors, overlap action, mediator masses and Steane code retain prior owners.
"""
import json, hashlib
from pathlib import Path
from itertools import combinations, product
import numpy as np
from scipy.linalg import expm
import w33_pass11438_11442_coercive_curved_recovery as N
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11443_11447_preserved_gauge_lorentz_recovery.json'


def gauge_fixed_space():
    import w33_pass11438_finite_native_model as F
    p=N.prior()['native_alignment'];phi=N.M.P.read(p['phi']).ravel()
    act=F.H@phi;R=np.r_[(-act.imag).T,act.real.T]
    z=np.linalg.svd(R,full_matrices=True)[2][-8:].T
    K=np.einsum('ah,aij->hij',z,F.H)
    w,E=np.linalg.eigh(np.einsum('hij,hjk->ik',K,K))
    return K,E[:,w<1e-9]


def preserved_vacuum(x=None,H=None):
    import w33_pass11438_finite_native_model as F
    from scipy.linalg import block_diag
    from scipy.optimize import minimize
    K,B=gauge_fixed_space();E=np.block([[B.real,-B.imag],[B.imag,B.real]]);T=block_diag(E,E)
    with F.native_context():
        if x is None:
            p=N.prior()['native_alignment'];start=T.T@F.pack(N.M.P.read(p['phi']).ravel(),N.M.P.read(p['psi']).ravel())
            def obj(y):
                v,g=F.fun(T@y);return v,T.T@g
            sol=minimize(obj,start,jac=True,method='L-BFGS-B',options=dict(maxiter=350,gtol=1e-8,ftol=1e-14,maxcor=30,maxls=30));x=T@sol.x
        v,g=F.fun(x);a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
        orbit=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
        rank=int(np.linalg.matrix_rank(orbit,tol=1e-8))
        assert rank==78 and v<-2.94 and np.linalg.norm(g)<5e-6
        # Identify closure and absence of family components, not a physical colour embedding.
        basis=K.reshape(8,-1).T
        closure=max(np.linalg.norm((1j*(u@v-v@u)).ravel()-basis@np.linalg.lstsq(basis,(1j*(u@v-v@u)).ravel(),rcond=None)[0]) for u in K for v in K)
        row=dict(status='PASS',coordinates=x.tolist(),value=v,gradient_norm=float(np.linalg.norm(g)),fixed_complex_dimension=B.shape[1],common_gauge_rank=rank,unbroken_generator_count=86-rank,fixed_residual=float(max(np.linalg.norm(K@a),np.linalg.norm(K@b))),lie_closure_residual=float(closure),CP_flux=N.bounded_cross(a.reshape(27,3),b.reshape(27,3))[1])
        if H is None:
            step=2e-5;eye=np.eye(324)
            H=np.column_stack([(F.fun(x+step*d)[1]-F.fun(x-step*d)[1])/(2*step) for d in eye]);H=(H+H.T)/2
        if H is not None:
            Q=np.linalg.svd(orbit,full_matrices=True)[0][:,rank:];ev=np.linalg.eigvalsh(Q.T@H@Q)
            row.update(normal_eigenvalues=ev.tolist(),negative_normal_count=int(sum(ev < -1e-3)),positive_normal_count=int(sum(ev>1e-3)))
            w,V=np.linalg.eigh(Q.T@H@Q);errors=[]
            for j in range(8):
                d=Q@V[:,j];re=(F.fun(x+1e-5*d)[1]-F.fun(x-1e-5*d)[1])/2e-5;errors.append(float(np.linalg.norm(re-H@d)))
            row['smaller_step_direction_errors']=errors
        row['scope']='A lower-energy stationary candidate in the common fixed space of eight compact generators. The subgroup is not yet identified with physical colour; transverse stability and full global optimality must be assessed separately.'
        return row


def relaxed_soft_modes():
    import w33_pass11438_finite_native_model as F
    p=json.loads(N.OUT.read_text())['coercive_alignment']['finite_field_witness'];x=np.array(p['coordinates']);H=np.array(p['full324_Hessian'])
    a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:]
    orbit=np.concatenate([(1j*(F.H@a)).real,(1j*(F.H@a)).imag,(1j*(F.H@b)).real,(1j*(F.H@b)).imag],axis=1).T
    Q=np.linalg.svd(orbit,full_matrices=True)[0][:,86:];w,V=np.linalg.eigh(Q.T@H@Q);hard=Q@V[:,4:];inv=(hard/w[4:])@hard.T
    rows=[]
    with F.native_context():
        val,g=F.fun(x)
        for j in range(4):
            d=Q@V[:,j]
            for h in [.02,.04]:
                delta=-inv@((F.fun(x+h*d)[1]+F.fun(x-h*d)[1])/2-g)
                raw=(F.fun(x+h*d)[0]+F.fun(x-h*d)[0])/2-val
                relaxed=(F.fun(x+h*d+delta)[0]+F.fun(x-h*d+delta)[0])/2-val-g@delta
                rows.append(dict(mode=j,h=h,raw_over_h4=raw/h**4,relaxed_over_h4=relaxed/h**4))
    return rows


def orbit_measure():
    # Exact determinant cocycle on gauge orbits; it deliberately does not supply
    # the missing local current on the quotient or transitions between flux sectors.
    charges=[1,-4,2,-3,6,0];mult=[6,3,3,2,1,1];rng=np.random.default_rng(11444);alpha=rng.normal(size=81)
    linear=sum(q*n for q,n in zip(charges,mult));cubic=sum(q**3*n for q,n in zip(charges,mult))
    cocycle=np.prod([np.exp(2j*q*alpha.sum())**n for q,n in zip(charges,mult)])
    bad=[1,1,-2];assert sum(bad)==0 and sum(q**3 for q in bad)!=0
    links=N.uniform_flux_links(34);u,v=links
    plaquette=u*np.roll(v,-1,axis=0)*np.roll(u.conj(),-1,axis=1)*v.conj()
    rows=[]
    for q in [1,2,3,4,6]:
        angle=np.angle(plaquette**q);rows.append(dict(charge=q,flux=float(angle.sum()/(2*np.pi)),max_plaquette=float(np.max(abs(1-plaquette**q)))))
    assert max(r['max_plaquette'] for r in rows)<1/30
    return dict(status='PASS',linear_anomaly=linear,cubic_anomaly=cubic,orbit_cocycle_error=float(abs(cocycle-1)),unit_flux_L=34,charged_fluxes=rows,negative_control=dict(charges=bad,linear=0,cubic=-6),scope='Explicit nonlocal gauge-orbit determinant cocycle and nonzero admissible flux links. The anomalous negative control also passes the orbit cocycle; hence this does not establish a local anomaly-safe measure, global quotient section or sector gluing. Luescher theorem remains external.')


def scalar_matrix(X,mass=0.):
    """Lorentzian affine finite element: half integral ((df)^2-m^2 f^2)."""
    A=np.c_[np.ones(5),X];grad=np.linalg.inv(A)[1:,:]
    vol=abs(np.linalg.det(X[1:]-X[0]))/24
    eta=np.diag([-1.,1,1,1]);massmat=vol*(np.ones((5,5))+np.eye(5))/30
    return vol*grad.T@eta@grad-mass**2*massmat,vol


def scalar_refinement(X,weights,mass):
    Y=np.vstack([X,np.asarray(weights)@X]);K=np.zeros((6,6));vol=0.
    for omitted in range(5):
        ids=[j for j in range(5) if j!=omitted]+[5];a,v=scalar_matrix(Y[ids],mass);K[np.ix_(ids,ids)]+=a;vol+=v
    assert abs(K[5,5])>1e-10
    return K[:5,:5]-np.outer(K[:5,5],K[5,:5])/K[5,5],vol,float(K[5,5])


def lorentzian_matter():
    boundary=np.array(N.M.previous()['curved_gauge_islands']['boundary_lengths']);X=N.M.P.local_simplex_coordinates(boundary,[0,1,2,3,4])
    eta=np.diag([-1.,1,1,1]);G=X[1:]@eta@X[1:].T;inertia=np.linalg.eigvalsh(G)
    boost=np.eye(4);boost[:2,:2]=[[np.cosh(.4),np.sinh(.4)],[np.sinh(.4),np.cosh(.4)]]
    rows=[]
    for weights in [np.ones(5)/5,np.array([.1,.15,.2,.25,.3])]:
        for m in [0.,.2]:
            coarse,vol=scalar_matrix(X,m);fine,v,pivot=scalar_refinement(X,weights,m)
            second,_,_=scalar_refinement(X@boost.T,weights,m)
            rows.append(dict(weights=weights.tolist(),mass=m,volume_error=abs(v-vol),coarse_difference=float(np.linalg.norm(fine-coarse)),boost_error=float(np.linalg.norm(second-fine)),center_pivot=pivot))
    assert sum(inertia<0)==1 and sum(inertia>0)==3
    assert max(r['coarse_difference'] for r in rows if r['mass']==0)<1e-10
    assert max(r['boost_error'] for r in rows)<1e-10
    return dict(status='PASS',coordinates=X.tolist(),Gram_eigenvalues=inertia.tolist(),rows=rows,scope='Off-shell Lorentzian scalar action and exact elimination of an interior scalar on a supplied Wick-metric native simplex. Signature is supplied. Massive Schur action changes under refinement. Lorentzian gravitational boost angles, causal branches, varying geometry and physical Newton coupling remain open.')


def threshold_running():
    heavy=json.loads((ROOT/'data/w33_pass11423_11427_flux_index_curved_uv_noise.json').read_text())['heavy_extension']
    # Q:240 Dirac doublets =>480 triplets; U,D:480; F_u,F_d:6.
    nf_heavy=80*3*2+80*3*2+3*2;assert heavy['heavy_Dirac_count']==3*nf_heavy
    p=N.M.previous();A=N.M.P.load()['A'];masses=np.sort(np.concatenate([N.M.full_quark_masses(A,N.M.P.read(p['pin_vacua'][k]),.81) for k in ['pin_up','pin_down']]))
    assert len(masses)==nf_heavy+6
    b=lambda n:11-2*n/3
    def evolve(mu0,mu1,g):
        thresholds=masses[(masses>mu0)&(masses<mu1)];points=np.r_[mu0,thresholds,mu1];inv=1/g**2
        for lo,hi in zip(points[:-1],points[1:]):
            n=int(sum(masses<=np.sqrt(lo*hi)));inv+=b(n)*np.log(hi/lo)/(8*np.pi**2)
        return inv
    rows=[]
    mu0=float(masses[-1]*1.01)
    for g in [.1,.2,.5,1.]:
        logratio=8*np.pi**2/(-b(len(masses))*g*g)
        rows.append(dict(supplied_g=g,log_Landau_ratio=logratio,Landau_ratio=float(np.exp(logratio))))
    return dict(status='PASS',heavy_representations=heavy['heavy_representations'],heavy_triplet_species=nf_heavy,total_triplet_species=len(masses),masses=masses.tolist(),high_energy_b0=b(len(masses)),normalization='beta(g)=-b0*g^3/(16*pi^2); T(fund)=1/2; all stated quarks elementary Dirac SU3 fundamentals; no additional coloured scalars/vectors.',examples=rows,threshold_inverse_g2_at_100=dict(g_at_5=.2,value=evolve(5,100,.2)),scope='Threshold-dependent one-loop gauge running of the supplied heavy extension, not observed coupling prediction. Negative b0 obstructs asymptotic freedom. Pole extrapolation is illustrative and fails near strong coupling. Full joint pin/link/Majorana vacuum and UV completion remain open.')


def native_clock():
    P=N.M.P;B=P.read(N.M.previous()['native_noise']['invariant_frame']);Q=B.conj().T@P.read(P.prior()['local_gates']['logical_basis']);D=np.array(P.load()['D'],float)
    Hp=B.conj().T@(D[:40].T@D[:40])@B;Hl=B.conj().T@(D[40:].T@D[40:])@B
    U=expm(-.02j*Hp)@expm(-.04j*Hl)@expm(-.02j*Hp);ph=np.angle(np.linalg.eigvals(U));nonzero=ph[abs(ph)>1e-9];rows=[]
    for m in [256,1024,4096]:
        false=max(abs(np.sin(m*nonzero/2)/(m*np.sin(nonzero/2)))**2)
        rows.append(dict(clock_powers=m,max_false_accept=float(false),envelope=float(1/(m*m*np.min(np.sin(nonzero/2)**2)))))
    return dict(phases=ph.tolist(),joint_dark_dimension=int(sum(np.linalg.eigvalsh(Hp+Hl)<1e-9)),logical_clock_residual=float(np.linalg.norm((U-np.eye(12))@Q)),zero_phase_detection=rows,scope='Fejer zero-phase acceptance distinguishes the two-dimensional logical dark space in this isolated native12-cell. Controlled clock, phase readout and reset remain supplied hardware primitives; repeated clock cost and noisy readout are not a threshold proof.')


def route_optimization():
    A=np.array(N.M.P.load()['A']);n=len(A);dist=np.where(A,1,1000);np.fill_diagonal(dist,0)
    for k in range(n):dist=np.minimum(dist,dist[:,k,None]+dist[None,k,:])
    weights=np.array([1,1,2,1,2,2,3]);best=None
    for hub in range(n):
        for flag in np.where(A[hub])[0]:
            ids=[j for j in range(n) if j not in [hub,flag]];ids.sort(key=lambda j:(dist[hub,j],j));data=np.zeros(7,dtype=int)
            for j,node in zip(np.argsort(-weights,kind='stable'),ids[:7]):data[j]=node
            costs=6*(dist[hub,data]-1)+1;value=int(6*(weights@costs)+36)
            if best is None or value<best['baseline_routed_CNOTs']:
                best=dict(hub=hub,flag=int(flag),data=data.tolist(),distances=dist[hub,data].tolist(),costs=costs.tolist(),baseline_routed_CNOTs=value,conditional_routed_CNOTs=int(value+2*(weights@costs)))
    # Degree4: flag occupies one neighbour; three heaviest data have weight7.
    assert best['baseline_routed_CNOTs']==6*(7+7*5)+36==288
    best.update(lower_bound=288,conditional_ticks=best['conditional_routed_CNOTs']*107635,scope='Optimal supplied independent-cell placement for the previous flagged extraction count and SWAP routing model. No noisy routing threshold or fault propagation guarantee follows.')
    return best


def erasure_decoder():
    # Hamming columns; CSS syndrome is (H z, H x). For <=2 known erased
    # positions, all 4^e Pauli branches have distinct syndromes.
    H=np.array([[(j>>k)&1 for j in range(1,8)] for k in range(3)])
    tables=[]
    for e in [1,2]:
        for sites in combinations(range(7),e):
            decoder={}
            for labels in product(range(4),repeat=e):
                x=np.zeros(7,dtype=int);z=x.copy()
                for j,p in zip(sites,labels):x[j]=p&1;z[j]=(p>>1)&1
                syndrome=tuple(np.r_[H@z%2,H@x%2]);assert syndrome not in decoder
                decoder[syndrome]=labels
            tables.append(dict(sites=list(sites),Pauli_branches=len(decoder)))
    return dict(status='PASS',supports=tables,process='Known erased sites reset to0; independently twirl each erased qubit to I/2; measure six Steane stabilizers; use erasure-mask-specific syndrome to invert the unique Pauli branch. Reset+twirl equals the erasure marginal channel; distinct syndromes establish deterministic correction.',native_clock=native_clock(),route=route_optimization(),scope='Operational Pauli/CSS recovery for one or two known erasures under supplied reset/twirl/readout; gate compilation prior-owned. Extra unknown faults violate the generic2t+e<3 guarantee. Native controlled spectroscopy and reset synthesis, routed noisy channels and threshold remain open.')


def run(x=None,H=None):
    result=dict(status='PASS',passes=[11443,11444,11445,11446,11447],reservation='bad2527a6',preserved_gauge=preserved_vacuum(x,H),relaxed_soft=relaxed_soft_modes(),measure=orbit_measure(),lorentzian=lorentzian_matter(),running=threshold_running(),recovery=erasure_decoder())
    names=['data/w33_pass11438_11442_coercive_curved_recovery.json','data/w33_pass11423_11427_flux_index_curved_uv_noise.json']
    result['source_sha256']={name:hashlib.sha256(json.dumps(json.loads((ROOT/name).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for name in names}
    OUT.write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--candidate');p.add_argument('--hessian');a=p.parse_args()
    x=np.load(a.candidate)['x'] if a.candidate else None;H=np.load(a.hessian)['H'] if a.hessian else None
    r=run(x,H);print({k:v.get('status','AUDIT') for k,v in r.items() if isinstance(v,dict) and k!='source_sha256'})
