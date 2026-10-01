#!/usr/bin/env python3
"""The new G26 T-magic orbit feeds the existing T port only with a retained frame.

Prior owners: balanced_g26_tmagic_vacuum owns the stationary point/orbit;
w33_qutrit_t_teleportation_port and Pass411 own injection and corrections.
The added result is the 72 -> 8 x 9 resource/frame contract and its lost-frame
entanglement-breaking channel. Pauli twirling itself is a standard theorem.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
import w33_20261001_balanced_g26_tmagic_vacuum as B
import w33_qutrit_t_teleportation_port as PORT
OUT=ROOT/'data/w33_20261001_tmagic_orbit_factory_contract.json'


def exact_pauli_twirl():
    """Verify the universal symbolic channel, including off-diagonal inputs."""
    w=sp.Symbol('w');mod=sp.Poly(w*w+w+1,w)
    def red(q):return sp.rem(sp.Poly(sp.expand(q),w),mod).as_expr()
    X=sp.Matrix([[0,0,1],[1,0,0],[0,1,0]])
    Z=sp.diag(1,w,w*w)
    r=sp.Matrix(3,3,sp.symbols('r0:9'))
    avg=sp.zeros(3)
    for a in range(3):
        for b in range(3):
            U=X**a*Z**b
            # w -> w^2 is complex conjugation; symbols in r are untouched.
            Ud=U.applyfunc(lambda q:q.subs(w,w*w)).T
            avg+=U*r*Ud/9
    avg=avg.applyfunc(red)
    assert avg==sp.trace(r)*sp.eye(3)/3
    return '1/9 sum_(a,b) X^a Z^b rho Z^-b X^-a = Tr(rho) I/3'


def to_choi(state):
    """SUM from state tensor |0>: sum_j state[j] |j,j>."""
    pair=np.zeros(9,complex)
    for j in range(3):pair[4*j]=state[j]
    return pair


def transfer(resource,a,b):
    beta=np.array(PORT.bell_state(a,b)).reshape(3,3)
    pair=resource.reshape(3,3)
    return np.einsum('iq,qj->ji',beta.conj(),pair)


def partial_transpose(rho):
    return rho.reshape(3,3,3,3).transpose(0,3,2,1).reshape(9,9)


def payload():
    twirl=exact_pauli_twirl()
    G,stab=B.projective_group()
    v=np.array(PORT.T_PLUS)
    X=np.array(PORT.X);Z=np.array(PORT.Z)
    paul=[np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b)
          for a in range(3) for b in range(3)]
    rays=[];frames=[]
    for i,g in enumerate(G):
        q=g@v
        if all(abs(abs(np.vdot(q,r))-1)>1e-8 for r in rays):
            rays.append(q);frames.append(i)
    assert len(rays)==72
    left=set(range(72));orbits=[]
    while left:
        j=min(left)
        orb=sorted(i for i in left if any(abs(abs(np.vdot(rays[i],p@rays[j]))-1)<1e-8 for p in paul))
        assert len(orb)==9
        orbits.append(orb);left.difference_update(orb)
    comps=json.loads((ROOT/'data/w33_pass11262_g26_mub_s4_yukawa_bridge.json').read_text())['mub_action']['components']
    axes=[]
    for q in rays:
        pur=[sum(abs(np.vdot(stab[j],q))**4 for j in c) for c in comps]
        hit=[i for i,p in enumerate(pur) if abs(p-1/3)<1e-8]
        assert len(hit)==1;axes.append(hit[0])
    orbit_axes=[];twirl_errors=[]
    for orb in orbits:
        assert len({axes[i] for i in orb})==1
        orbit_axes.append(axes[orb[0]])
        rho=sum(np.outer(rays[i],rays[i].conj()) for i in orb)/9
        twirl_errors.append(float(np.max(abs(rho-np.eye(3)/3))))
    assert [orbit_axes.count(i) for i in range(4)]==[2]*4
    # Full operator identities, independent of input ket or external entanglement.
    residual=0.0;choi_error=0.0
    for q,i in zip(rays,frames):
        canonical=G[i].conj().T@q
        resource=to_choi(canonical)
        choi_error=max(choi_error,float(np.max(abs(resource-np.array(PORT.t_choi_resource())))))
        for a in range(3):
            for b in range(3):
                K=transfer(resource,a,b);C=np.array(PORT.correction(a,b))
                residual=max(residual,float(np.max(abs(C@K-np.array(PORT.T_GATE)/3))))
    assert residual<1e-12 and choi_error<1e-12
    # Forgetting the Pauli frame gives a separable pair, not the pure T-Choi pair.
    lost=np.zeros((9,9),complex)
    for j in range(3):lost[4*j,4*j]=1/3
    assert np.min(np.linalg.eigvalsh(partial_transpose(lost)))>=0
    # Direct Bell-channel replay on all nine matrix units proves dephasing.
    channel_error=0.0
    for i in range(3):
        for j in range(3):
            E=np.zeros((3,3),complex);E[i,j]=1
            got=np.zeros((3,3),complex)
            for k in range(3):
                pair=np.zeros(9,complex);pair[4*k]=1/np.sqrt(3)
                for a in range(3):
                    for b in range(3):
                        C=np.array(PORT.correction(a,b));K=C@transfer(pair,a,b)
                        got+=K@E@K.conj().T
            want=E if i==j else np.zeros((3,3),complex)
            channel_error=max(channel_error,float(np.max(abs(got-want))))
    assert channel_error<1e-12
    return {
      'schema':'w33.20261001.tmagic-orbit-factory-contract.v1',
      'status':'PASS_G26_T_MAGIC_RESOURCE_REQUIRES_PAULI_FRAME_AND_ALL_72_FRAMES_FEED_EXISTING_T_PORT',
      'prior_owners':['analysis/w33_20261001_balanced_g26_tmagic_vacuum.py',
                      'analysis/w33_qutrit_t_teleportation_port.py',
                      'analysis/w33_pass411_qutrit_magic_injection.py'],
      'orbit':{'rays':72,'pauli_orbits':8,'pauli_orbit_sizes':[9]*8,
               'pauli_orbits_per_MUB_axis':[2]*4,'pauli_orbit_axes':orbit_axes,
               'frame_records':[{'ray':j,'projective_clifford_index':i,'MUB_axis':axes[j]}
                                for j,i in enumerate(frames)]},
      'exact_theorem':twirl,
      'frame_loss':{'each_pauli_orbit_average':'I/3',
                   'axis_only_average':'I/3',
                   'even_pauli_orbit_label_without_displacement':'I/3',
                   'choi_resource':'1/3 sum_j |jj><jj|',
                   'implemented_channel':'rho -> diag(rho)',
                   'entanglement_breaking':True,
                   'meaning':'A selected MUB axis or scalar invariant value does not specify usable magic. The relative Pauli displacement must be retained or fixed by a physical reference.'},
      'retained_frame':{'operation':'apply the stored Clifford inverse to the selected ray, then SUM(state tensor |0>)',
                        'resource':'(I tensor T)|Phi3>',
                        'gate':'T=diag(1,zeta9,zeta9^-1)',
                        'bell_branches_checked':648,
                        'probability_per_branch':'1/9',
                        'corrected_operator':'T/3 for every ray and Bell branch',
                        'operator_error':residual,'choi_preparation_error':choi_error},
      'numerical_controls':{'maximum_pauli_orbit_twirl_error':max(twirl_errors),
                            'dephasing_channel_operator_error':channel_error},
      'boundary':'This connects a conditional prepared ray to the prior ideal injection protocol. The barrier supplies an exact local stationary orbit, not a demonstrated state factory, cooling dynamics, noise threshold, or laboratory reference frame. The 72-ray census is numerical; the universal Pauli-twirl identity is symbolic.',
      'checks':{'pauli_twirl_symbolic':True,'eight_orbits_of_nine':True,
                'two_pauli_orbits_per_axis':True,'all648_branch_operator_maps':True,
                'frame_loss_dephases_input':True}
    }


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    p=payload();txt=json.dumps(p,indent=2,sort_keys=True)+'\n'
    if a.check:
        assert json.loads(OUT.read_text())==p
    else:OUT.write_text(txt)
    print(json.dumps({'status':p['status'],'rays':72,'branches':648,'frame_loss':'dephasing'}))
if __name__=='__main__':main()
