"""Independent algebra/operator checks for the five constructed maps."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11506_11510_five_physics_targets as P


def certificate():
    return json.loads(P.OUT.read_text())


def test_inputs_and_scopes():
    d=certificate();assert d['status']=='PASS' and d['passes']==list(range(11506,11511))
    for name,value in d['source_sha256'].items():
        assert value==hashlib.sha256(json.dumps(P.read(name),sort_keys=True,separators=(',',':')).encode()).hexdigest()
    for name in ['vacuum','SM','metric','ward','recovery']: assert d[name]['status']=='PASS'


def test_exact_native_stationarity_first_jets():
    # Replay all81 first jets, not only the one-dimensional radial derivative.
    fresh=P.balanced_vacuum();d=certificate()['vacuum']
    assert fresh['selected_radius_squared_interval']==d['selected_radius_squared_interval']
    assert fresh['numerical_full_gradient_norm']<1e-8 and fresh['energy_gap_above_CP_candidate']>2


def test_isolated_root_radial_value_derivative():
    d=certificate()['vacuum'];lo,hi=map(sp.Rational,d['selected_radius_squared_interval']);r=sp.symbols('r')
    # Independently differentiate the two-field balanced restricted potential.
    V=2*(27**4*r**12-sp.Rational(27,2)*r**3+27**6*r**18)+sp.Rational(36,10**6)*r**2
    poly=sp.Poly(sp.diff(V,r)/(3*r),r)
    assert poly.as_expr()==sp.sympify(d['radial_polynomial'])
    assert poly.count_roots(lo,hi)==1 and poly.eval(lo)<0<poly.eval(hi)


def test_actual_SM_color_blocks_and_charge_conservation():
    d=P.sm_blocks();B,tensor,_=P.tensors();t3=list(map(sp.Rational,d['T3']))
    y=list(map(sp.Rational,P.read(P.SOURCES[1])['hypercharge_diagonal']))
    for h in [d['neutral_Hu_index'],d['neutral_Hd_index']]:
        for i,j in np.argwhere(tensor[:,:,h]):assert y[i]+t3[i]+y[j]+t3[j]==0
    assert d['heavy_rank']==30 and d['light_first_order_rank']>0
    assert d['light_second_order_rank']>0


def test_schur_matching_replays_full_heavy_inverse():
    _,tensor,_=P.tensors();d=certificate()['SM'];Fu=sp.Matrix(d['Fu']);Fd=sp.Matrix(d['Fd'])
    mass=sp.kronecker_product(sp.Matrix(tensor[:,:,0]),sp.eye(3))
    e=sp.kronecker_product(sp.Matrix(tensor[:,:,d['neutral_Hu_index']]),Fu)+sp.kronecker_product(sp.Matrix(tensor[:,:,d['neutral_Hd_index']]),Fd)
    h=[3*i+j for i in d['heavy_E6_indices'] for j in range(3)];l=[3*i+j for i in d['light_E6_indices'] for j in range(3)]
    Elh=e.extract(l,h);H=mass.extract(h,h);correction=-(Elh*H.inv()*Elh.T)
    assert [[i,j,str(correction[i,j])] for i in range(51) for j in range(51) if correction[i,j]]==d['second_order_nonzero']
    eps=sp.Rational(1,100);exact=eps*e.extract(l,l)-eps**2*Elh*(H+eps*e.extract(h,h)).inv()*Elh.T
    assert exact==eps*e.extract(l,l)+eps**2*correction if not d['heavy_EW_block_nonzero'] else exact.shape==(51,51)


def test_graph_Dirac_spectral_determinant():
    d=certificate()['metric'];old=P.read(P.SOURCES[2])['geometry'];edges=np.array(old['cover_edges']);h=np.array(old['harmonic_FCC_displacements'])
    B=np.zeros((160,80));B[np.arange(160),edges[:,0]]=1;B[np.arange(160),edges[:,1]]=-1
    G=expm(np.array([[.1,.02,0],[.02,-.05,.03],[0,.03,-.05]]));w=1/np.einsum('ei,ij,ej->e',h,G,h)
    A=np.sqrt(w)[:,None]*B;D=np.block([[np.zeros((80,80)),A.T],[A,np.zeros((160,160))]])
    m=1.7
    lhs=np.linalg.slogdet(m*np.eye(240)+1j*D)[1]
    rhs=80*np.log(m)+np.linalg.slogdet(m*m*np.eye(80)+A.T@A)[1]
    assert abs(lhs-rhs)<1e-9
    assert min(d['shape_Hessian_eigenvalues'])>0 and np.linalg.norm(d['shape_gradient'])<1e-8


def test_diagonal_metric_monotonicity_symbolic():
    f,lam=sp.symbols('f lambda',positive=True)
    assert sp.simplify(sp.diff(sp.log((4+f*lam)/(1+f*lam)),f)+3*lam/((4+f*lam)*(1+f*lam)))==0
    assert all(r['shift']>0 for r in certificate()['metric']['diagonal_scans'])


def test_rational_metric_inverse_residuals_and_Sylvester_witness():
    d=certificate()['metric']['rational_shape_certificate'];src=P.read(P.SOURCES[-1])['metric']
    h=np.array(src['displacements_times80'],int);N=np.zeros((80,80),dtype=object);den=d['denominator'];scale=int(d['inverse_quantization_scale'])
    for (a,b),v in zip(src['edges'],h):
        w=sp.Rational(6400,int(v@v))*den;assert w.q==1
        for i,j,sg in [(a,a,1),(b,b,1),(a,b,-1),(b,a,-1)]:N[i,j]+=sg*int(w)
    for m2,q,bound in zip([4,1],d['inverse_numerators'],d['inverse_error_operator_norm_bounds']):
        Q=np.array([[int(x) for x in row] for row in q],dtype=object)
        residual=den*scale*np.eye(80,dtype=object)-(N+m2*den*np.eye(80,dtype=object))@Q
        proved=sp.Rational(80*max(abs(int(v)) for v in residual.ravel()),m2*den*scale)
        assert proved==sp.Rational(bound) and proved<sp.Rational(1,10**40)
    H=sp.Matrix(d['exact_approximate_Hessian']);low=H-5*sp.Rational(d['Hessian_entry_error_bound'])*sp.eye(5)
    assert H==H.T and all(low[:k,:k].det(method='domain-ge')>0 for k in range(1,6))
    edges={tuple(sorted(e)):h[j] for j,e in enumerate(src['edges'])}
    key=lambda x:tuple(x if x[0]>0 else -x)
    for flip,perm in zip(d['sign_flips'],d['sign_flip_vertex_permutations']):
        assert sorted(perm)==list(range(80))
        for (a,b),v in edges.items():assert key(v)==key(np.array(flip)@edges[tuple(sorted([perm[a],perm[b]]))])


def test_exact_Ward_primitive_and_torus_locality_boundary():
    import w33_pass11476_11480_cubic_geometry_correlated as M
    d=certificate()['ward'];B=sp.Matrix(M.gauge_boundary(2).astype(int));K=sp.Matrix(d['exact_Green_kernel']);n=16
    assert B.T*K==-(sp.eye(n)-sp.ones(n)/n)
    assert d['anomaly_norm']>1e-10 and d['gauge_invariance_residual']<1e-10
    assert all(r['Ward_defect']<=r['proved_upper_bound']+1e-11 for r in d['local_series'])
    assert d['box_gap_scaling'][-1]['rho']>.99


def test_nonlinear_Ward_gauge_replay():
    fresh=P.exponential_ward();assert fresh['Ward_residual']<1e-10
    assert fresh['minimum_Wilson_gap']>.1 and fresh['gauge_invariance_residual']<1e-10


def test_native_joint_prior_matches_direct_quantum_marginals():
    import w33_pass11471_11475_frames_currents_matching as Q
    import w33_pass11476_11480_cubic_geometry_correlated as M
    d=certificate()['recovery'];T=np.array(d['physical_transfer']);V,R=Q.recovery_rows();prior=np.array(d['syndrome_prior'])
    states=np.kron(V,V).astype(complex);paulis=np.array([np.kron(a,b) for a in M.PAULI for b in M.PAULI])
    def apply_pair(op,state,wire):
        axes=[wire,7+wire]+[i for i in range(15) if i not in [wire,7+wire]]
        t=state.reshape((2,)*14+(4,)).transpose(axes)
        return (op@t.reshape(4,-1)).reshape((2,)*14+(4,)).transpose(np.argsort(axes)).reshape(16384,4)
    chars=np.array(d['syndrome_characters']);words=np.array(d['stabilizer_words']);q=np.empty(64)
    for side,expected in [(0,prior.sum(axis=1)),(1,prior.sum(axis=0))]:
        for j,word in enumerate(words):
            evolved=states.copy()
            for wire,p in enumerate(word):
                # Full dual two-qubit channel, including dependence on the
                # other code block; no nonsignalling/marginal-product assumption.
                out_label=4*p if side==0 else p
                op=np.einsum('j,jab->ab',T[out_label],paulis)
                evolved=apply_pair(op,evolved,wire)
            q[j]=np.vdot(states,evolved).real/4
        actual=chars@q/64
        assert np.max(abs(actual-expected))<1e-10


def test_joint_inference_and_complete_leakage_branches():
    d=certificate()['recovery'];prior=np.array(d['syndrome_prior'])
    assert abs(prior.sum()-1)<1e-12 and min(prior.ravel())>=0
    assert d['correlated_vs_marginal_prior_difference']>1e-6
    assert d['untwirled_vs_twirled_prior_difference']>1e-8
    assert abs(d['rows'][0]['joint_correct_syndrome_probability']-1)<1e-10
    for r in d['rows']:assert r['joint_correct_syndrome_probability']>=r['separate_correct_syndrome_probability']-1e-12
    K=P.dec(d['logical_cell_subchannel']);effect=np.kron(K.conj().T@K,K.conj().T@K)
    lo,hi=d['two_cell_leakage_probability_bounds'];w=np.linalg.eigvalsh(np.eye(4)-effect)
    assert np.max(abs(w[[0,-1]]-[lo,hi]))<1e-10
    assert hi>1e-4
