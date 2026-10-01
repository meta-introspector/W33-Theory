#!/usr/bin/env python3
"""Exact sparse D-flat plane, universal mirror law, and normal-slice mass completion.

Prior owners: Pass11260/61 semisimple classification/coordinates, Pass11269
mirrors, and the dated global-selector and dimension-five mass packets.
The new plane replaces a numerical normalization witness, not its old gauge map.
"""
from pathlib import Path
import sys, json, argparse
from collections import Counter
import numpy as np
import sympy as s
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_pass11261_g26_cartan_coordinate_map as P
import w33_pass11255_11259_g26_common as G
from w33_20261001_isovolume_dirac_gravity_response import compare_certificate
OUT=ROOT/'data/w33_20261001_exact_cartan_mirror_completion.json'
SUPPORT=((0,52,80),(5,37,75),(26,40,54))

def plane():
    Q=np.zeros((81,3),dtype=int)
    for j,indices in enumerate(SUPPORT):Q[list(indices),j]=1
    return Q

def operators():
    old=P._load(ROOT/'analysis/w33_pass11218_e8_cubic_cartan_kernel.py','exact_sparse_old')
    records=old.load_parent().build_records()
    def op(v):
        A=np.zeros((81,81),dtype=np.result_type(v,int))
        for u,x,o,c in records:A[o,x]+=c*v[u]
        return A
    return [op(q) for q in plane().T],op

def components(adj):
    unseen=set(range(len(adj)));out=[]
    while unseen:
        seed=min(unseen);unseen.remove(seed);block=[seed]
        for i in block:
            for j in np.flatnonzero(adj[i]):
                if int(j) in unseen:unseen.remove(int(j));block.append(int(j))
        out.append(sorted(block))
    return out

def exact_blocks(Ms):
    adj=sum(abs(m) for m in Ms);blocks=components(adj);rows=[];pair_counts=Counter()
    for indices in blocks:
        ix=np.array(indices);n=len(ix);sub=[m[np.ix_(ix,ix)] for m in Ms]
        if n==3:
            active=[j for j,m in enumerate(sub) if np.any(m)]
            assert len(active)==1
            A=sub[active[0]]
            assert np.array_equal(A.T@A,3*np.eye(3,dtype=int)-np.ones((3,3),int))
            rows.append({'indices':indices,'size':3,'coordinate':active[0],
                         'proof':'cycle skew block: eigenvalues of A^T A are 0,3,3'})
            continue
        color={0:0};queue=[0]
        for i in queue:
            for j in np.flatnonzero(adj[np.ix_(ix,ix)][i]):
                if int(j) not in color:color[int(j)]=1-color[i];queue.append(int(j))
                assert color[int(j)]!=color[i]
        l=[i for i in range(n) if color[i]==0];r=[i for i in range(n) if color[i]==1]
        assert len(l)==len(r)==n//2
        Ps=[m[np.ix_(l,r)] for m in sub];active=[j for j,p in enumerate(Ps) if np.any(p)]
        for p in (Ps[j] for j in active):assert np.array_equal(p.T@p,np.eye(n//2,dtype=int))
        Us=[Ps[active[0]].T@Ps[j] for j in active[1:]]
        if n==6:
            assert len(active)==2
            U=Us[0];assert np.array_equal(np.linalg.matrix_power(U,3),-np.eye(3,dtype=int))
            assert np.trace(U)==np.trace(U@U)==0
            pair_counts[tuple(active)]+=1
            proof='relative signed permutation has each eigenvalue -omega^k once'
        else:
            assert n==18 and active==[0,1,2]
            U,V=Us
            assert np.array_equal(U@V,V@U)
            assert np.array_equal(U@U@U,np.eye(9,dtype=int))
            assert np.array_equal(V@V@V,np.eye(9,dtype=int))
            for a in range(3):
                for b in range(3):
                    assert np.trace(np.linalg.matrix_power(U,a)@np.linalg.matrix_power(V,b))==(9 if a==b==0 else 0)
            proof='regular Z3 x Z3 character: each joint character occurs once'
        rows.append({'indices':indices,'size':n,'active_coordinates':active,
                     'left_indices':ix[l].tolist(),'right_indices':ix[r].tolist(),
                     'relative_signed_permutations':[u.tolist() for u in Us],'proof':proof})
    assert sorted(len(b) for b in blocks)==[3]*3+[6]*9+[18]
    assert pair_counts=={(0,1):3,(0,2):3,(1,2):3}
    return rows

def gradient_data():
    z=s.Symbol('zeta');mod=s.Poly(s.cyclotomic_poly(9,z),z)
    red=lambda f:s.rem(s.Poly(s.expand(f),z),mod).as_expr()
    conj=lambda f:red(f.subs(z,z**8))
    raw=s.Matrix([1,z,z**8]);at=dict(zip(G.VARS,raw))
    gradients=[s.Matrix([red(s.diff(f,x).subs(at)) for x in G.VARS]) for f in G.invariants()[:2]]
    C=s.Matrix.hstack(raw.applyfunc(conj),*gradients)
    gram=(C.applyfunc(conj).T*C).applyfunc(red)
    assert gram==s.diag(3,3888,17714700)
    assert [gram[i+1,i+1]/3**(d-1) for i,d in enumerate([6,12])]==[16,100]
    assert red(C.det())!=0
    gs=[np.array([complex(s.N(x.subs(z,s.exp(2*s.pi*s.I/9)))) for x in g])/3**((d-1)/2) for g,d in zip(gradients,[6,12])]
    return gs,{'raw_gradient_vectors':[list(map(str,g)) for g in gradients],
      'raw_three_covector_Gram':list(map(str,gram.diagonal())),
      'unit_T_three_covector_Gram':[1,16,100],
      'EFT_operators':['(c5/Lambda)(Phi^dag chi)(Phi^dag psi)',
        '(c6/Lambda^9)(dI6[chi])(dI6[psi])',
        '(c12/Lambda^21)(dI12[chi])(dI12[psi])'],
      'light_masses_at_T':['|c5|r^2/Lambda','16|c6|r^10/Lambda^9','100|c12|r^22/Lambda^21'],
      'relative_to_heavy_scale_r':['(r/Lambda)^1','16(r/Lambda)^9','100(r/Lambda)^21'],
      'definition_of_I':'Explicit signed-cubic contraction invariants from w33_20261001_global_e6_cartan_covariants.py, independently proved to restrict to u6,u12 on this plane.',
      'scope':'This certificate independently proves the normal-slice Gram and separation; the global polynomial circuits and analytic differential evaluator are certified in the companion global-covariant packet. No measured masses or EFT coefficients are fixed.'}

def payload():
    Q=plane();Ms,op=operators();B=np.load(ROOT/'artifacts/e6_27rep_basis_export/E6_basis_78.npy')
    assert np.max(abs(B-np.rint(B)))==0
    B=B.real.astype(int);V=Q.reshape(27,3,3)
    assert np.array_equal(Q.T@Q,3*np.eye(3,dtype=int))
    for b in B:
        assert not np.any(np.einsum('ica,ij,jcb->ab',V,b,V))
    for j in range(3):
        for k in range(3):assert np.array_equal(V[:,:,j].T@V[:,:,k],np.eye(3,dtype=int)*int(j==k))
    for m in Ms:assert not np.any(m@Q) and np.array_equal(m.T,-m)
    blocks=exact_blocks(Ms);gs,grad=gradient_data()
    E=Q/np.sqrt(3);z=np.exp(2j*np.pi/9);q=np.array([1,z,z.conjugate()])/np.sqrt(3)
    phi=E@q;M=op(phi);covectors=[phi.conj()]+[E@g for g in gs]
    coefficients=[.001,.002/16,.003/100]
    lift=sum(c*np.outer(v,v) for c,v in zip(coefficients,covectors))
    lifted=M+lift;old=np.linalg.svd(M,compute_uv=False);new=np.linalg.svd(lifted,compute_uv=False)
    assert np.max(abs(M.conj().T@lift))<1e-12
    assert np.max(abs(old[:78]-new[:78]))<1e-12
    assert np.max(abs(np.sort(new[78:])-np.array([.001,.002,.003])))<1e-12
    assert np.linalg.matrix_rank(lifted,tol=1e-10)==81
    gen=B[0]+B[0].T;U27=expm(.17j*gen);U3=expm(.23j*np.diag([1,-1,0]));U=np.kron(U27,U3)
    transformed=op(U@phi);cov=float(np.max(abs(transformed-U.conj()@M@U.conj().T)))
    assert cov<1e-12
    return {'schema':'w33.20261001.exact-cartan-mirror-completion.v1','status':'PASS_EXACT_SPARSE_DFLAT_MIRROR_LAW_AND_NORMAL_SLICE_MASS_COMPLETION',
      'plane':{'support':list(map(list,SUPPORT)),'unnormalized_Gram':(Q.T@Q).tolist(),
        'exact_E6_cross_moments_zero':True,'exact_SU3_cross_moments_zero':True,
        'exact_pairwise_brackets_zero':True,'kinetic_isometry':'E=Q/sqrt(3), E^dag E=I',
        'Cartan_justification':'Moment-zero implies a closed complex orbit, hence semisimple in the theta representation; generic bracket rank78 gives the minimal three-dimensional centralizer. Use the prior rank-three G26 classification.',
        'old_plane_boundary':'Alternative exact canonical plane; not an algebraic certificate of the old large numerical gauge transformation.'},
      'mass_blocks':blocks,
      'universal_mass_law':{'scope':'Exact for every complex q on Phi=Qq/sqrt(3), including degeneracies',
        'squared_singular_values':'|q_i|^2 twice; |q0+omega^a q1+omega^b q2|^2/3 twice; |qi-omega^k qj|^2/3 six times; three zeros',
        'mirror_reading':'12 normalized MUB overlaps twice and 9 normalized SIC overlaps times2/3 six times',
        'trace_squared':'20||q||^2','trace_fourth':'8||q||^4','previous_scope':'Five numerical directions in the dated full-normalization packet; this packet certifies the law on the new exact plane.'},
      'covariant_completion':grad,
      'finite_completion_control':{'rank_before_after':[78,81],'new_light_masses_sorted':np.sort(new[78:]).tolist(),
        'heavy_mass_error':float(np.max(abs(old[:78]-new[:78]))),'unitary_gauge_covariance_error':cov},
      'boundary':['Global I6/I12 polynomial circuits are supplied by the companion global-covariant packet.','This nonsupersymmetric EFT completion does not fix Standard Model assignments, scales, coupling values, or radiative stability.'],
      'prior_owners':['analysis/w33_pass11260_exact_semisimple_cartan.py','analysis/w33_pass11261_g26_cartan_coordinate_map.py','analysis/w33_pass11269_g26_qutrit_dictionary.py','analysis/w33_20261001_g26_global_polynomial_vacuum.py','analysis/w33_20261001_dflat_cartan_mass_loop.py'],
      'external_sources':['https://m.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=2123&wshow=paper','https://arxiv.org/abs/hep-th/9506098']}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();p=payload()
    if args.check:assert compare_certificate(p,json.loads(OUT.read_text()))
    else:OUT.write_text(json.dumps(p,indent=2,sort_keys=True)+'\n')
    print(p['status'])
