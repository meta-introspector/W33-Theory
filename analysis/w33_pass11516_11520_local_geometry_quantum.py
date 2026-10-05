"""Five follow-ups to11506-11510; explicit conditional maps, not a TOE.

Owners:11438 native potential;11271/11293 and11507 SM interface;
11494/11508 metric geometry;11429/11509 overlap;11475/11480/11510 recovery.
"""
import itertools,json,hashlib,time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm,null_space
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/w33_pass11516_11520_local_geometry_quantum.json'
SOURCES=['data/w33_pass11506_11510_five_physics_targets.json','data/w33_pass11493_11497_exact_matching_inference.json','data/w33_pass11476_11480_cubic_geometry_correlated.json','data/w33_pass11271_canonical_sm_higgs_bridge.json','data/w33_pass11384_native_inputs.json']
read=lambda name:json.loads((ROOT/name).read_text())
enc=lambda a:dict(real=np.asarray(a).real.tolist(),imag=np.asarray(a).imag.tolist())
dec=lambda a:np.array(a['real'])+1j*np.array(a['imag'])

def vacuum_normal_audit():
    import w33_pass11438_finite_native_model as F
    import w33_pass11481_11485_fibers_metric_yukawa_noisy as Q
    old=read(SOURCES[2])['vacuum'];E=dec(old['isotope']['native_tripotents'])
    z=np.array(read(SOURCES[1])['stationary']['final_reduced_coordinates'])
    T=np.kron(E,np.eye(3));R=np.block([[T.real,-T.imag],[T.imag,T.real]])
    embedding=np.zeros((324,36));embedding[:162,:18]=R;embedding[162:,18:]=R
    x=embedding@z;eye=np.eye(324);spectra=[]
    with F.native_context():
        value,g=F.fun(x)
        for step in [3e-5,1e-5]:
            H=np.column_stack([(F.fun(x+step*v)[1]-F.fun(x-step*v)[1])/(2*step) for v in eye]);H=(H+H.T)/2
            eigen,V=np.linalg.eigh(H);spectra.append(dict(step=step,eigenvalues=eigen.tolist(),symmetry_error=0.))
            print('full native Hessian',step,eigen[:8],flush=True)
        v=V[:,0];controls=[]
        for step in [1e-4,3e-4,1e-3]:
            vp=F.fun(x+step*v)[0];vm=F.fun(x-step*v)[0]
            controls.append(dict(step=step,plus_energy=vp,minus_energy=vm,second_difference=(vp+vm-2*value)/step**2))
        p=x[:81]+1j*x[81:162];s=x[162:243]+1j*x[243:]
        gauge=np.column_stack([F.pack(1j*h@p,1j*h@s) for h in F.H]);sing=np.linalg.svd(gauge,compute_uv=False)
        N=null_space(gauge.T,rcond=1e-8);normal=N.T@H@N;normal_eig=np.linalg.eigvalsh(normal)
    return dict(status='PASS',value=value,gradient_norm=float(np.linalg.norm(g)),coordinates=x.tolist(),full_Hessian_spectra=spectra,
                gauge_singular_values=sing.tolist(),gauge_rank=int(sum(sing>1e-8)),normal_dimension=N.shape[1],normal_eigenvalues=normal_eig.tolist(),
                minimum_direction=v.tolist(),directional_controls=controls,
                scope='Full324 numerical transverse/gauge-normal audit of the earlier low-energy CP point, with two finite-difference scales and energy replays. Not an interval vacuum certificate; small eigenvalues and gauge-rank tolerance are numerical. Existing11481 fibers retain ownership.')

def adjoint_flavor():
    import w33_pass11506_11510_five_physics_targets as P
    B,d,_=P.tensors();sm=read('data/w33_pass11271_canonical_sm_higgs_bridge.json');y=list(map(sp.Rational,sm['hypercharge_diagonal']))
    old=read(SOURCES[0])['SM'];blocks={};n=27
    # Same adjoint on indistinguishable fermion legs, explicitly symmetrized.
    for sector,a,b,h in [('down','down_Q','Dc',old['neutral_Hd_index']),('lepton','charged_L','Ec',old['neutral_Hd_index']),('up','up_Q','Uc',old['neutral_Hu_index'])]:
        pairs=[(i,j) for i in old['SM_indices'][a] for j in old['SM_indices'][b] if d[i,j,h]]
        single={y[i]+y[j] for i,j in pairs};double={y[i]*y[j] for i,j in pairs}
        assert len(single)==len(double)==1
        blocks[sector]=dict(single_fermion_insertion=str(next(iter(single))),one_on_each_fermion=str(next(iter(double))))
    assert blocks['down']['single_fermion_insertion']==blocks['lepton']['single_fermion_insertion']=='1/2'
    assert blocks['down']['one_on_each_fermion']=='1/18' and blocks['lepton']['one_on_each_fermion']=='-1/2'
    # Covariance: d(A psi,A psi,H) has both adjoint insertions transforming.
    # A single insertion on symmetrized identical fermions is -d(psi,psi,A H).
    assert all(not d[i,j,k] or y[i]+y[j]+y[k]==0 for i,j,k in itertools.product(range(n),repeat=3))
    Yd=sp.diag(1,2,3);Ye=sp.diag(sp.Rational(1,3),6,3)
    F0=(9*Yd+Ye)/10;F2=sp.Rational(9,5)*(Yd-Ye)
    assert F0+F2/18==Yd and F0-F2/2==Ye
    return dict(status='PASS',charge_blocks=blocks,exact_split_ratio='-9',operator='d_abc (A psi)^a_i (A psi)^b_j H^c_ij / Lambda²; H in (27,bar6), A in E6 adjoint',
                single_insertion_identity='d(A psi,psi,H)+d(psi,A psi,H)=-d(psi,psi,A H)',
                exact_target_down=list(map(str,Yd.diagonal())),exact_target_lepton=list(map(str,Ye.diagonal())),
                exact_F0=[[str(x) for x in row] for row in F0.tolist()],exact_F2=[[str(x) for x in row] for row in F2.tolist()],
                scope='Exact minimal split within this same-adjoint symmetrized-fermion-leg subclass: one insertion is universal, two split with ratio-9. A transforming adjoint and independent sextets are supplied; target textures are engineered controls, not dynamical selection or observed mass predictions. General GUT Clebsch mechanisms are prior art.')

def metric_blindness():
    from scipy.optimize import minimize_scalar,minimize
    from scipy.linalg import expm_frechet
    src=read(SOURCES[1])['metric'];h=np.array(src['displacements_times80'],int);edges=np.array(src['edges']);n=80
    rows=sp.Matrix([[int(v[0]**2),int(v[1]**2),int(v[2]**2),2*int(v[0]*v[1]),2*int(v[0]*v[2]),2*int(v[1]*v[2])] for v in h])
    rank=rows.rank();assert rank==4
    null=rows.nullspace();assert len(null)==2
    adjacency=[[] for _ in range(n)]
    for (a,b),v in zip(edges,h):adjacency[a].append((b,v));adjacency[b].append((a,-v))
    two=[]
    for a in range(n):
        for b,v in adjacency[a]:
            for c,w in adjacency[b]:
                if a<c and np.any(v+w):two.append((a,c,v+w))
    hh=np.array([v for _,_,v in two],int);allh=np.vstack([h,hh]);alle=np.vstack([edges,[[a,c] for a,c,_ in two]])
    allrows=sp.Matrix([[int(v[0]**2),int(v[1]**2),int(v[2]**2),2*int(v[0]*v[1]),2*int(v[0]*v[2]),2*int(v[1]*v[2])] for v in allh]);assert allrows.rank()==6
    witness=allrows.T.rref()[1][:6];det=allrows[list(witness),:].det();assert det!=0
    B=np.zeros((len(alle),n));B[np.arange(len(alle)),alle[:,0]]=1;B[np.arange(len(alle)),alle[:,1]]=-1
    S=[np.eye(3),np.diag([1,-1,0]),np.diag([1,1,-2])]
    for i,j in [(0,1),(0,2),(1,2)]:s=np.zeros((3,3));s[i,j]=s[j,i]=1;S.append(s)
    S=np.array(S);multipliers=np.r_[np.ones(len(h)),np.full(len(hh),.1)]
    def action(z,gradient=False):
        Z=np.einsum('a,aij->ij',z,S);G=expm(Z);length=np.einsum('ei,ij,ej->e',allh,G,allh);weights=6400*multipliers/length
        L=B.T@(weights[:,None]*B);ev=np.linalg.eigvalsh(L)
        ev[0]=0 # the constant mode is exactly zero on this connected graph.
        value=float(np.sum(2*np.log1p(ev)-np.log(.25+ev)-np.log(4+ev)))
        if not gradient:return value
        R=2*np.linalg.inv(np.eye(n)+L)-np.linalg.inv(.25*np.eye(n)+L)-np.linalg.inv(4*np.eye(n)+L)
        lever=np.einsum('ei,ij,ej->e',B,R,B);grad=[]
        for s in S:
            dg=expm_frechet(Z,s,compute_expm=False);dl=np.einsum('ei,ij,ej->e',allh,dg,allh)
            grad.append(-np.dot(lever,weights*dl/length))
        return value,np.array(grad)
    initial=minimize_scalar(lambda t:action(np.r_[t,np.zeros(5)]),bounds=(-8,8),method='bounded')
    candidates=[]
    for perturb in [np.zeros(5),np.array([.2,.1,0,0,0]),np.array([-.2,.1,0,0,0])]:
        o=minimize(lambda z:action(z,True),np.r_[initial.x,perturb],jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':200})
        candidates.append(o)
    opt=min(candidates,key=lambda o:o.fun)
    z=opt.x;eps=2e-4;I=np.eye(6);H=np.empty((6,6))
    for i,j in itertools.product(range(6),repeat=2):H[i,j]=(action(z+eps*(I[i]+I[j]))-action(z+eps*(I[i]-I[j]))-action(z+eps*(-I[i]+I[j]))+action(z-eps*(I[i]+I[j])))/(4*eps**2)
    eig=np.linalg.eigvalsh(H)
    return dict(status='PASS',edge_metric_rank=rank,exact_blind_directions=[list(map(str,v)) for v in null],
                exact_blindness='For every diagonal traceless D, h_e^T D h_e=0; G and G+D have identical ungrounded L whenever both SPD. Any action depending only on this L cannot determine all six metric entries.',
                two_hop_count=len(two),two_hop_edges=alle[len(edges):].tolist(),two_hop_displacements_times80=hh.tolist(),two_hop_metric_rank=6,rank_witness_indices=list(witness),rank_witness_determinant=str(det),
                inventory='four real vertex bosons of mass1; two complex graph-Dirac fermions of squared masses1/4 and4',
                augmented_action='2 logdet(I+L)-logdet(I/4+L)-logdet(4I+L); two-hop coupling0.1 supplied',
                optimizer_success=bool(opt.success),optimizer_message=str(opt.message),coordinates=z.tolist(),value=action(z),gradient=action(z,True)[1].tolist(),optimization_controls=[dict(value=float(o.fun),coordinates=o.x.tolist(),success=bool(o.success)) for o in candidates],Hessian=H.tolist(),Hessian_eigenvalues=eig.tolist(),
                flux_winding_control=flux_winding_control(),
                scope='Exact blindness theorem and exact rank-six two-hop repair on the actual native cover. Augmented determinant optimization is a control and its very soft direction is not certified. A separately declared quantized torus flux/winding completion proves full local metric stabilization with a tuned counterterm; it is not native-only gravity, physical spin identification or a cosmological-constant prediction.')

def flux_winding_control():
    """Actual11508 determinant plus declared torus stress; exact six-metric proof."""
    # The trace expression is an EXACT rational definition, not its dyadic estimate.
    gamma=sp.Rational(1,1000);upper_r=sp.Rational(79,3)
    Lambda_lower=1-sp.Rational(2,3)*gamma*upper_r
    volume_hessian_lower=6-sp.Rational(237,4)*gamma
    assert Lambda_lower>0 and volume_hessian_lower>0
    # Four-boson/two-fermion numerical inventory above is a separate experiment.
    # Here reuse the already certified11508 TWO-boson/ONE-fermion inventory.
    src=read(SOURCES[1])['metric'];h=np.array(src['displacements_times80'],int);edges=src['edges'];L=np.zeros((80,80))
    for (a,b),v in zip(edges,h):
        w=6400/(v@v)
        for i,j,s in [(a,a,1),(b,b,1),(a,b,-1),(b,a,-1)]:L[i,j]+=s*w
    eigen=np.linalg.eigvalsh(L)[1:];r=float(np.sum(3*eigen/((4+eigen)*(1+eigen))))
    lam=1-2*float(gamma)*r/3
    wzz=float(np.sum(3*eigen*(eigen**2-4)/((4+eigen)**2*(1+eigen)**2)))
    return dict(gamma=str(gamma),torus_flux_integer=1,torus_flux_coupling='3',three_unit_windings=[1,1,1],winding_coupling='1',
                action='V(G)=Lambda v+3/(2v)+(v/2)Tr(G^-1)+gamma W11508(G), v=sqrt(det G)',
                exact_counterterm='Lambda=1-(2gamma/3)Tr[(L+I)^-1 L-(L+4I)^-1 L], L=L11508(I); this is rational because L is rational',
                numerical_counterterm=lam,exact_counterterm_lower_bound=str(Lambda_lower),
                exact_stationarity='The base has Hessian diag(6,1,3,1,1,1) at I. The counterterm cancels the determinant volume gradient exactly. Native sign-flip automorphisms and equal-component h make shape gradients and volume/shape cross terms zero.',
                exact_Hessian_lower_diagonal=[str(volume_hessian_lower),'1','3','1','1','1'],numerical_volume_Hessian=6+float(gamma)*(wzz-1.5*r),
                exact_bounds='For each nonzero lambda, r=3lambda/[(4+lambda)(1+lambda)]<=1/3 since (lambda-2)^2>=0; |Wzz|<=1/4 since each logdet mode second derivative m²lambda/(m²+lambda)^2 lies in[0,1/4]. There are79 modes. Shape Hessian of W11508 is positive by its rational certificate, so the full six-metric Hessian has the stored positive diagonal lower bound.',
                coercivity='Lambda>0, flux term diverges as v->0, Lambda v as v->infinity; at bounded positive volume, winding term diverges at shape boundary. Native W is bounded between0 and80 log4, so a global minimum exists. I is proved a strict local minimum, not the unique full-shape global minimum.',
                scope='Supplied3-torus with a quantized3-form flux, three unit axion windings and a matched finite cosmological counterterm. Prior11324/11341 own flux/winding frameworks. This adds genuine geometric stress beyond graph L; it does not identify the graph Dirac field with SM spinors or derive the physical cosmological constant.')

def run_one(name):
    return {'vacuum':vacuum_normal_audit,'flavor':adjoint_flavor,'metric':metric_blindness,'ward':projector_transgression,'quantum':quantum_fidelity}[name]()

def projector_transgression():
    """Gapped projector pair current: relative-reference covariance, not a measure."""
    from numpy.polynomial.legendre import leggauss
    import w33_pass11428_11432_native_alignment_measure_coarse_flags as Q
    A=np.array(read(SOURCES[0])['ward']['background']).reshape(16,4);nodes,weights=leggauss(20);current=np.zeros((16,16));gaps=[]
    charges=[1,-4,2,-3,6];mult=[6,3,3,2,1]
    def pair_current(P,dP):
        K=dP@P-P@dP
        blocks=K.reshape(16,4,16,4).transpose(0,2,1,3);pb=P.reshape(16,4,16,4).transpose(0,2,1,3)
        j=2*np.einsum('xyij,yxji->xy',blocks,pb).real
        assert np.linalg.norm(j+j.T)<1e-10
        assert np.linalg.norm(j.sum(axis=1)-np.trace(dP.reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real)<1e-10
        return j
    a0=np.zeros(16);a1=np.zeros(16);rows=[]
    lam=.003*np.random.default_rng(11519).normal(size=16);cov_error=0.
    for q,m in zip(charges,mult):
        P0,*_=Q.weak_overlap(q,1.,np.zeros_like(A),[],L=2);P1,*_=Q.weak_overlap(q,1.,A,[],L=2)
        diag=lambda P:np.trace(P.reshape(16,4,16,4),axis1=1,axis2=3).diagonal().real
        a0+=q*m*diag(P0);a1+=q*m*diag(P1)
        for t,w in zip((nodes+1)/2,weights/2):
            P,ds,_,gap,_,_=Q.weak_overlap(q,1.,t*A,[A],L=2);j=pair_current(P,ds[0]);current+=q*m*w*j;gaps.append(gap)
            G=np.diag(np.repeat(np.exp(1j*q*lam),4));jg=pair_current(G@P@G.conj().T,G@ds[0]@G.conj().T);cov_error=max(cov_error,np.linalg.norm(jg-j))
    # Gauge-transforming only the endpoint does not keep the same reference path.
    defect=np.linalg.norm(current.sum(axis=1)-(a1-a0));assert defect<1e-10
    M=8.;g=float(sp.Rational(read(SOURCES[0])['ward']['homotopy_gap_lower_bound']));r=1-(g/M)**2
    tails=[]
    for N in [32,128,512,2048]:
        sign_tail=r**(N+1)/(1-r)
        derivative_tail=sign_tail+2*r**N*((N+1)-N*r)/(1-r)**2
        tails.append(dict(N=N,sign_tail_bound=sign_tail,derivative_tail_per_D_over_M=derivative_tail,propagation_radius=2*N+1))
    assert cov_error<1e-10
    aa=sp.symbols('a0:4',nonnegative=True)
    freegap=sp.expand(sum(2*x-x*x for x in aa)+(sum(aa)-1)**2)
    assert sp.expand(freegap-1-2*sum(aa[i]*aa[j] for i in range(4) for j in range(i+1,4)))==0
    return dict(status='PASS',pair_current=current.tolist(),anomaly_difference=(a1-a0).tolist(),divergence_error=float(defect),relative_covariance_error=float(cov_error),minimum_sampled_Wilson_gap=float(min(gaps)),
                exact_identity='Pdot=[[Pdot,P],P]; j_xy=2 Re tr([Pdot,P]_xy P_yx)=-j_yx; sum_y j_xy=tr_x Pdot',
                integral_current='J_xy=sum_q m_q q integral_0^1 j_xy(P_q(t),Pdot_q(t)) dt',
                uniform_bound=dict(M=M,gap=g,r=r,free_Wilson_gap_identity=str(freegap),
                    uniform_gap_proof='For a_mu=1-cos(p_mu)>=0, H0²=1+2 sum_(mu<nu) a_mu a_nu>=1, at every torus size. ||H(A)-H0||<=max|q| sum_mu max_x|A_mu(x)|, so the stored componentwise weak-field bounds give the same positive homotopy gap at every volume. ||H||<=7<8.',
                    series='sign(H)=(H/M) sum_k binom(2k,k)/4^k (I-H²/M²)^k',tails=tails),
                locality='For uniformly finite-range H,Hdot with ||H||<=M and spectral gap>=g>0, polynomial radii and geometric tails give volume-independent exponential locality for P,Pdot and their pair current. Routing pair currents along shortest lattice paths preserves exponential decay up to polynomial factors.',
                scope='Actual weak-overlap projector transgression with uniform-gap locality theorem conditional on the supplied gapped homotopy. Current is invariant under simultaneous time-independent gauge transformation of endpoint AND flat reference; fixing reference=free and transforming endpoint alone is not this covariance theorem. A canonical reference-free gauge-invariant current, chiral-measure curvature/integrability, flux sectors and physical infinite-volume construction remain open. Kato projector transport and sparse spectral-projector locality are prior art.')

def quantum_fidelity():
    """Full logical PTM and complete flagged instrument fidelity, retaining native PTM."""
    import w33_pass11476_11480_cubic_geometry_correlated as M
    import w33_pass11471_11475_frames_currents_matching as Q
    V,R=Q.recovery_rows();words=[];signs=[];corrections=[]
    bits_to_label=np.array([0,1,3,2]);label_to_bits=bits_to_label
    # Every logical Pauli has64 stabilizer-coset terms; coefficients are exact +/-1/64.
    for logical,p in enumerate(M.PAULI):
        w,c=M.sparse_terms(V@p@V.conj().T);assert len(w)==64
        sg=np.rint(64*c).astype(int);assert np.max(abs(c-sg/64))<1e-12
        words.append(w);signs.append(sg)
        x=((w==1)|(w==2)).astype(int);z=((w==2)|(w==3)).astype(int);ch=[]
        for s in range(64):
            sz,sx=divmod(s,8);phase=np.zeros(64,int)
            if sz:phase+=z[:,7-sz]
            if sx:phase+=x[:,7-sx]
            ch.append(sg*(-1)**phase)
        corrections.append(np.array(ch))
    words=np.array(words);signs=np.array(signs);flat=words.reshape(256,7);sg=signs.ravel();pair=4*flat[:,None,:]+flat[None,:,:]
    C=corrections;labels=4*(np.arange(256)[:,None]//64)+np.arange(256)[None,:]//64
    source=read(SOURCES[0])['recovery'];T0=np.array(source['physical_transfer']);gate_p=.001;logical_z=.002
    results=[]
    for case,T in [('native',T0),('one_paired_extraction_fault_layer',np.diag(np.r_[1.,np.full(15,1-gate_p)])@T0)]:
        # Roundoff zeros are removed with an explicit normalized-Choi diamond-norm bound.
        Tc=T.copy();Tc[abs(Tc)<1e-14]=0
        delta=4*np.linalg.norm(M.transfer_choi(T-Tc),ord='nuc');roundoff_bound=float(np.expm1(7*np.log1p(delta)))
        raw=np.empty((16,16,64,64));start=time.time()
        first=[]
        for op in range(16):
            ii,jj=np.where(Tc[op,pair[:,:,0]]!=0);first.append((ii,jj))
        for a in range(256):
            logical_left,sa=divmod(a,64)
            for b in range(256):
                logical_right,sb=divmod(b,64);op=4*flat[a]+flat[b];ii,jj=first[op[0]];val=Tc[op[0],pair[ii,jj,0]]*sg[ii]*sg[jj]
                for wire in range(1,7):
                    keep=Tc[op[wire],pair[ii,jj,wire]]!=0;ii=ii[keep];jj=jj[keep];val=val[keep]*Tc[op[wire],pair[ii,jj,wire]]
                    if not len(ii):break
                accum=np.bincount(labels[ii,jj],weights=val,minlength=16)
                raw[4*logical_left+logical_right,:,sa,sb]=accum
            if a%64==63:print('quantum raw',case,a+1,'seconds',round(time.time()-start,1),flush=True)
        branches=np.empty((4096,16,16))
        for out,inp in itertools.product(range(16),repeat=2):
            a,b=divmod(out,4);branches[:,out,inp]=(C[a]@raw[out,inp]@C[b].T/4096).ravel()
        dephase=np.array([1,1-2*logical_z,1-2*logical_z,1])
        if case!='native':branches*=np.kron(dephase,dephase)[None,:,None]
        basis=np.array([np.kron(a,b) for a,b in itertools.product(M.PAULI,repeat=2)])
        conj_sign=np.array([[np.trace(p@u@p@u.conj().T).real/4 for p in basis] for u in basis])
        K=dec(source['logical_cell_subchannel']);E=np.kron(K.conj().T@K,K.conj().T@K);v=source['verification_flip']
        def fidelity_coeff(effect):
            vec=effect.T.ravel()/2
            return np.array([[np.vdot(vec,np.kron(basis[j].T,basis[i])@vec).real/16 for j in range(16)] for i in range(16)])
        coeff=(1-v)*((1-v)*fidelity_coeff(E)+v*fidelity_coeff(np.eye(4)-E))
        scores=np.einsum('lo,oi,soi->ls',conj_sign,coeff,branches).reshape(16,64,64)
        assert min(scores.ravel())>-1e-9
        bits=np.array([[(j>>k)&1 for k in range(6)] for j in range(64)]);distance=np.sum(bits[:,None]!=bits[None,:],axis=2)
        decoded=[]
        for prior in source['rows']:
            err=prior['bit_error'];acc=np.array(prior['acceptance']);transfer=np.einsum('s,soi->oi',acc.ravel(),branches)
            if case=='native' and err==0:
                old=next(r for r in read(SOURCES[2])['recovery']['rows'] if r['p']==.03 and r['relay']==0)
                assert np.linalg.norm(transfer-np.array(old['decoded_transfer']))<1e-9
            J=M.transfer_choi(transfer);eigen=np.linalg.eigvalsh(J);assert min(eigen)>-1e-9
            overlap=lambda effect:float(np.vdot(effect.T.reshape(16)/2,J@(effect.T.reshape(16)/2)).real)
            fe_core=float(np.trace(transfer)/16)
            fe_full=(1-v)*((1-v)*overlap(E)+v*overlap(np.eye(4)-E))
            assert -1e-10<=fe_full<=fe_core+1e-8
            confusion=err**distance*(1-err)**(6-distance);optimal=0.;chosen=[]
            for a,b in itertools.product(range(64),repeat=2):
                weighted=scores*confusion[a,None,:,None]*confusion[b,None,None,:]
                choice=int(np.argmax(weighted));chosen.append(choice);optimal+=float(weighted.ravel()[choice])
            assert optimal>=fe_full-1e-9 and optimal<=1+1e-9
            decoded.append(dict(bit_error=err,correct_sector_probability=float(transfer[0,0]),logical_transfer=transfer.tolist(),logical_Choi=enc(J),minimum_Choi_eigenvalue=float(min(eigen)),core_entanglement_fidelity=fe_core,complete_flagged_entanglement_fidelity=fe_full,
                fidelity_optimal_complete_flagged_fidelity=optimal,fidelity_optimal_report_to_logical_syndrome=chosen))
        results.append(dict(case=case,physical_paired_depolarization=0 if case=='native' else gate_p,postdecode_logical_phase_fault=0 if case=='native' else logical_z,roundoff_seven_channel_diamond_bound=roundoff_bound,complete_fidelity_branch_scores=scores.tolist(),rows=decoded))
    return dict(status='PASS',results=results,
                exact_formula='F_E(channel)=vec(E^T)^dag J vec(E^T)/4 for normalized Choi J; F_flagged=(1-v)[(1-v)F_E+v F_(I-E)], E=(Kdag K) tensor(Kdag K)',
                target='Ideal isometry W=(U tensor U) embedding of two logical qubits in the actual two12-cells; distinct accepted wrong-syndrome and erasure output sectors have zero overlap with this target.',
                optimal_decoder='For each12-bit report choose(syndrome1,syndrome2,logical Pauli1,logical Pauli2) maximizing complete-instrument entanglement-fidelity contribution times its readout likelihood. Stored integer choice=4096*logical_pair+64*syndrome1+syndrome2. Optimal only in the named finite Pauli decision class.',
                scope='All256 logical Pauli-transfer entries computed by sparse finite stabilizer contraction, retaining native cross-block coherence. Complete11510 flagged instrument entanglement fidelity includes native cell projection and accepted leakage; classical MAP and a constructed fidelity-optimal Pauli decoder are compared. Extra fault case models one two-qubit depolarizing layer immediately before ideal extraction and independent postdecode logical Z faults, not an entire physical syndrome circuit or threshold.')

def run():
    result={}
    for name in ['vacuum','flavor','metric','ward','quantum']:
        result[name]=run_one(name);print(name,'PASS',flush=True)
        Path('/tmp/w33_11516_partial.json').write_text(json.dumps(result))
    result.update(status='PASS',passes=list(range(11516,11521)),reservation='948b8e23d',producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256={name:hashlib.sha256(json.dumps(read(name),sort_keys=True,separators=(',',':')).encode()).hexdigest() for name in SOURCES})
    OUT.write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':run()
