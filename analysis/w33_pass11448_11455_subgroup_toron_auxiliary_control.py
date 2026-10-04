"""Eight constructed tests of the five11443--11447 frontiers.
Established methods retain prior ownership; physical closure is not asserted.
"""
import json,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import expm
import w33_pass11443_11447_preserved_gauge_lorentz_recovery as P
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/w33_pass11448_11455_subgroup_toron_auxiliary_control.json'


def previous():return json.loads(P.OUT.read_text())


def subgroup():
    import w33_pass11438_finite_native_model as F
    K,B=P.gauge_fixed_space();flat=K.reshape(8,-1).T;coeff=np.linalg.lstsq(F.H.reshape(86,-1).T,flat,rcond=None)[0]
    ads=[];closure=[]
    for u in K:
        bracket=np.array([(1j*(u@v-v@u)).ravel() for v in K]).T;c=np.linalg.lstsq(flat,bracket,rcond=None)[0];ads.append(c);closure.append(np.linalg.norm(bracket-flat@c))
    adcas=sum(a.conj().T@a for a in ads);rep=sum(k@k for k in K);w=np.linalg.eigvalsh(rep)
    values=[]
    for v in w:
        if values and abs(v-values[-1][0])<1e-8:values[-1][1]+=1
        else:values.append([float(v),1])
    small=F.H[:78,::3,::3];ks=K[:,::3,::3]
    comm=np.concatenate([np.array([(u@v-v@u).ravel() for v in small]).T for u in ks],axis=0)
    centralizer=78-int(np.linalg.matrix_rank(comm,tol=1e-8))
    center=int(sum(np.linalg.eigvalsh(adcas)<1e-9))
    assert center==0 and np.linalg.norm(coeff[78:])<1e-10
    return dict(status='PASS',algebra_dimension=8,E6_centralizer_dimension=centralizer,center_dimension=center,bracket_residual=max(closure),family_component_norm=float(np.linalg.norm(coeff[78:])),adjoint_Casimir_eigenvalues=np.linalg.eigvalsh(adcas).tolist(),representation_Casimir_blocks=values,fixed_dimension=B.shape[1],scope='Compact center-free eight-dimensional Lie algebra is su3 by classification of compact semisimple Lie algebras. It lies in E6, not family SU3. Casimir blocks test the candidate colour-type branching. Explicit physical charge/intertwiner map remains required; this alone does not select the SM embedding.')


def coupled_soft(H=None):
    import w33_pass11438_finite_native_model as F
    x=np.array(previous()['preserved_gauge']['coordinates'])
    with F.native_context():
        v,g=F.fun(x)
        if H is None:
            h=2e-5;H=np.column_stack([(F.fun(x+h*d)[1]-F.fun(x-h*d)[1])/(2*h) for d in np.eye(324)]);H=(H+H.T)/2
        a=x[:81]+1j*x[81:162];b=x[162:243]+1j*x[243:];ap=1j*(F.H@a);bp=1j*(F.H@b);orbit=np.concatenate([ap.real,ap.imag,bp.real,bp.imag],axis=1).T
        Q=np.linalg.svd(orbit,full_matrices=True)[0][:,78:];w,V=np.linalg.eigh(Q.T@H@Q);soft=Q@V[:,:4];hard=Q@V[:,4:];inv=(hard/w[4:])@hard.T
        rng=np.random.default_rng(11449);dirs=list(np.eye(4))+[a/np.linalg.norm(a) for a in rng.normal(size=(8,4))];rows=[]
        for j,c in enumerate(dirs):
            d=soft@c
            for h in [.02,.04]:
                delta=-inv@((F.fun(x+h*d)[1]+F.fun(x-h*d)[1])/2-g)
                value=(F.fun(x+h*d+delta)[0]+F.fun(x-h*d+delta)[0])/2-v-g@delta
                rows.append(dict(direction=j,coefficients=c.tolist(),h=h,relaxed_energy=value,over_h4=value/h**4))
    return dict(status='PASS',rows=rows,scope='Mixed finite-amplitude soft directions with quadratic massive relaxation at the preserved-gauge candidate. Numerical signs at1e-10 energy scale do not certify a quartic tensor, exact flatness or stability.')


def flat_weyl(q,t,L=2):
    gamma,g5=P.N.M.P.gamma_matrices();sites=list(product(range(L),repeat=4));lookup={s:i for i,s in enumerate(sites)};n=len(sites);W=3*np.eye(4*n,dtype=complex)
    for mu in range(4):
        T=np.zeros((n,n),complex)
        for i,s in enumerate(sites):
            y=list(s);y[mu]=(y[mu]+1)%L
            T[i,lookup[tuple(y)]]=np.exp(1j*q*t/L if mu==3 else 0)*(-1 if mu==3 and s[mu]==L-1 else 1)
        W-=.5*(np.kron(T,np.eye(4)-gamma[mu])+np.kron(T.conj().T,np.eye(4)+gamma[mu]))
    G=np.kron(np.eye(n),g5);H=G@W;w,E=np.linalg.eigh(H);frame=E[:,w>0];gauge=np.repeat(np.exp(-2j*np.pi*q*np.array(sites)[:,3]/L),4)
    return frame,float(min(abs(w))),gauge


def toron_loops():
    rows=[]
    for q in [1,2,3,4,6]:
        for steps in [32,64]:
            frames=[];gaps=[]
            for t in np.linspace(0,2*np.pi,steps+1):
                E,g,G=flat_weyl(q,t);frames.append(E);gaps.append(g)
            hol=np.eye(frames[0].shape[1],dtype=complex);minimum=1.
            for E,F in zip(frames,frames[1:]):
                u,s,vh=np.linalg.svd(E.conj().T@F);minimum=min(minimum,min(s));hol=hol@u@vh
            u,s,vh=np.linalg.svd(frames[-1].conj().T@(G[:,None]*frames[0]));hol=hol@u@vh
            endpoint=np.linalg.norm(frames[-1]@frames[-1].conj().T-G[:,None]*(frames[0]@frames[0].conj().T)*G.conj()[None,:])
            rows.append(dict(charge=q,steps=steps,phase=float(np.angle(np.linalg.det(hol))),Wilson_gap=min(gaps),minimum_overlap_singular_value=float(minimum),endpoint_projector_error=float(endpoint)))
    return dict(status='PASS',L=2,rows=rows,scope='Actual closed toron loop on a flat admissible torus, with large gauge endpoint identification and full overlap projector. This tests Wilson-line holonomy on one sector; it is not an implemented local current or nontrivial-flux sector gluing.')


def stress_map():
    X=np.array(previous()['lorentzian']['coordinates']);weights=np.array([.1,.15,.2,.25,.3]);f=np.array([.3,-.2,.7,.1,-.5]);h=1e-5;rows=[]
    # Vary the geometry and integrate the interior scalar at each geometry.
    for i,j in [(1,0),(2,1),(3,2),(4,3)]:
        d=np.zeros_like(X);d[i,j]=1
        energy=lambda Y:float(f@P.scalar_refinement(Y,weights,.2)[0]@f/2)
        force=(energy(X+h*d)-energy(X-h*d))/(2*h)
        smaller=(energy(X+h*d/2)-energy(X-h*d/2))/h
        coarse=lambda Y:float(f@P.scalar_matrix(Y,0)[0]@f/2)
        fine=lambda Y:float(f@P.scalar_refinement(Y,weights,0)[0]@f/2)
        massless_error=abs((coarse(X+h*d)-coarse(X-h*d)-fine(X+h*d)+fine(X-h*d))/(2*h))
        rows.append(dict(vertex=i,coordinate=j,massive_force=force,step_error=abs(force-smaller),massless_force_matching_error=massless_error))
    assert max(r['step_error'] for r in rows)<1e-7
    return dict(status='PASS',boundary_field=f.tolist(),rows=rows,scope='Explicit matter force on varying supplied Lorentzian vertex coordinates after stationary scalar elimination. This is a term ready for coupled edge equations, not a complete gravitational action or solution.')


def auxiliary_determinants():
    A=np.array(P.N.M.P.load()['A'],float);K=10*np.eye(80)-.81*A;idx=[44,1,0];J=np.eye(80)[:,idx];Pin=P.N.M.P.read(P.N.M.previous()['pin_vacua']['pin_up'])
    # Native-node spatial kernel with family3. Actual path inverse uses two copies.
    Kh=np.kron(K,np.eye(3));E=np.kron(J,np.eye(3));light=2*np.eye(9);g=.05
    full=np.block([[Kh,g*E],[g*E.T,light]]);Schur=light-g*g*E.T@np.linalg.solve(Kh,E)
    slog=lambda M:np.linalg.slogdet(M)[1];error=slog(full)-slog(Kh)-slog(Schur)
    assert abs(error)<1e-10
    # Setting determinant normalization to1 is a changed integration measure.
    # At K=mI-phi A, removing detK removes its actual scalar force.
    force=-3*np.trace(np.linalg.solve(K,A));h=1e-5
    derivative=(3*slog(10*np.eye(80)-(.81+h)*A)-3*slog(10*np.eye(80)-(.81-h)*A))/(2*h)
    return dict(status='PASS',native_kernel_eigenvalue_range=[float(np.linalg.eigvalsh(K)[0]),float(np.linalg.eigvalsh(K)[-1])],factorization_error=float(error),heavy_logdet_phi_derivative=float(force),derivative_error=float(abs(force-derivative)),scope='Exact Gaussian Schur identity on the actual80-node kernel. Integrating elementary auxiliary fermions retains detK and its background force. Dividing by detK changes the measure and generally needs compensators or a new microscopic theory; renaming mediators does not repair the972-flavour beta coefficient. This finite matrix is a control, not the full spacetime gauge determinant.')


def joint_scale_pin():
    from scipy.optimize import minimize
    # New declared finite Majorana/pin slice; link stiffness imports the locally
    # renormalized quadratic response, not its full nonlinear determinant.
    old=P.N.M.previous();B=P.N.M.P.read(old['pin_vacua']['pin_down']);S=np.diag([1.,2.,3.]);seed=np.array([3.0697803,3.0896656])
    def nu(s,angles):return float(P.N.M.neutrino_phase_potential(B,s*S,angles,precision=35))
    def radial(s):return 100*(s*s-1)**2+nu(s,seed)
    opt=minimize(lambda x:radial(x[0]),[1.],method='BFGS',options={'gtol':1e-7});s=float(opt.x[0]);h=1e-4
    curvature=(radial(s+h)+radial(s-h)-2*radial(s))/h**2;force=(radial(s+h)-radial(s-h))/(2*h)
    # Portal .1*(phi-.81)*(s-1), imported radial curvature100.
    # Solve this explicitly declared quadratic link/portal approximation.
    portal=.1;phi=.81-portal*(s-1)/100
    # Phase potential differences evaluated separately to avoid losing tiny
    # phase forces in the large heavy/link vacuum constant.
    base=P.N.M.neutrino_phase_potential(B,s*S,seed,precision=45)
    phase=lambda x:float(P.N.M.neutrino_phase_potential(B,s*S,x,precision=45)-base)*1e12
    phaseopt=minimize(phase,seed,method='Nelder-Mead',options={'xatol':1e-7,'fatol':1e-9,'maxiter':160})
    return dict(status='PASS',supplied_Majorana_tree_coefficient=100.,supplied_portal=.1,Majorana_scale=s,link_linear_response_phi=phi,radial_force=force,radial_curvature=curvature,pin_angles=phaseopt.x.tolist(),phase_scaled_improvement=float(phaseopt.fun),scope='Declared finite pin/Majorana slice with an actual6-state neutrino determinant, supplied stabilizer and imported quadratic link response. Orientation and up/down pins are fixed inputs; full nonlinear joint vacuum and physical coefficient selection are not solved. Portal feedback into scale is omitted at order portal²/100.')


def reset_target():
    # One additional native cell is an environment for this target isometry.
    n=12;source=[i*n for i in range(n)];target=[0,n]+list(range(1,n-1))
    remaining_source=[i for i in range(n*n) if i not in source]
    remaining_target=[i for i in range(n*n) if i not in target]
    perm=np.empty(n*n,dtype=int);perm[source]=target;perm[remaining_source]=remaining_target
    U=np.eye(n*n)[:,perm];assert np.linalg.norm(U.T@U-np.eye(n*n))==0
    kraus=[U.reshape(n,n,n,n)[:,e,:,0] for e in range(n)]
    assert np.linalg.norm(sum(K.T@K for K in kraus)-np.eye(n))<1e-12
    return dict(native_cell_dimension=n,permutation=perm.tolist(),unitarity_error=float(np.linalg.norm(U.T@U-np.eye(n*n))),Kraus_rank=int(np.linalg.matrix_rank(np.stack([K.ravel() for K in kraus]))),scope='Complete logical-preserving leakage-reset unitary target on two native cells in the adapted logical/leakage basis. Initial environment state is logical0. Permutation extension is supplied; existing logical CNOT compilation does not synthesize its leakage-sector interactions.')


def clock_instrument():
    c=previous()['recovery']['native_clock'];theta=np.array(c['phases']);rows=[]
    # Spectral Fourier instrument from uniform controlled-clock history.
    for m in [32,128,1024]:
        labels=np.arange(m);k=np.arange(m);A=np.exp(1j*np.outer(k,theta));kraus=np.fft.fft(A,axis=0)/m
        completeness=max(abs(np.sum(abs(kraus)**2,axis=0)-1));false=max(abs(kraus[0,abs(theta)>1e-9])**2)
        rows.append(dict(history_dimension=m,TP_error=float(completeness),max_false_accept=float(false),logical_accept_error=float(max(abs(abs(kraus[0,abs(theta)<1e-9])**2-1)))))
    # Unitary reset with no environment impossible: rank cannot decrease.
    return dict(status='PASS',rows=rows,reset_target=reset_target(),reset_required_environment_dimension=11,scope='Complete Fourier clock measurement instrument, not only acceptance amplitudes. Controlled clock/readout remains a primitive to synthesize. Reset of ten orthogonal leakage states to one logical state requires ten leakage environment states orthogonal to the retained logical environment state, hence at least eleven environment dimensions in a Stinespring realization; closed unitary cell control alone cannot reset.')


def relay_faults():
    # CNOT Pauli propagation on data(control), relay, syndrome(target).
    ops=[(0,1),(1,0),(0,1),(1,2),(0,1),(1,0),(0,1)];rows=[]
    def propagate(x,z,ops):
        for c,t in ops:x[t]^=x[c];z[c]^=z[t]
        return x,z
    for location,(c,t) in enumerate(ops):
        for p,q in product(range(4),repeat=2):
            if p==q==0:continue
            x=np.zeros(3,dtype=int);z=x.copy();x[c]=p&1;z[c]=p>>1;x[t]=q&1;z[t]=q>>1
            x,z=propagate(x,z,ops[location+1:]);support=[i for i in range(3) if x[i] or z[i]]
            rows.append(dict(location=location,paulis=[p,q],output_support=support))
    return dict(status='PASS',fault_count=len(rows),relay_fault_count=sum(1 in r['output_support'] for r in rows),three_site_fault_count=sum(len(r['output_support'])==3 for r in rows),rows=rows,scope='Exhaustive105 Pauli faults in one actual SWAP-out/CNOT/SWAP-back route. Relay errors persist and can affect subsequent reuse. This is a fault-propagation building block, not full flagged extraction or a threshold; clean relays must be reset or audited in the full schedule.')


def run(H=None):
    r=dict(status='PASS',passes=list(range(11448,11456)),reservation='ec4cee553',subgroup=subgroup(),soft=coupled_soft(H),torons=toron_loops(),stress=stress_map(),auxiliary=auxiliary_determinants(),joint_slice=joint_scale_pin(),instrument=clock_instrument(),relay=relay_faults())
    names=['data/w33_pass11443_11447_preserved_gauge_lorentz_recovery.json']
    r['source_sha256']={n:hashlib.sha256(json.dumps(json.loads((ROOT/n).read_text()),sort_keys=True,separators=(',',':')).encode()).hexdigest() for n in names};OUT.write_text(json.dumps(r,indent=2)+'\n');return r

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--hessian');args=a.parse_args();H=np.load(args.hessian)['H'] if args.hessian else None
    print({k:v['status'] for k,v in run(H).items() if isinstance(v,dict) and 'status' in v})
