#!/usr/bin/env python3
"""Native parabolic harmonic periods -> connected rank-three Levi covers.

Prior owners: BT861 (Steinberg H1); Pass11070 (the two radicals);
Pass5677 (voltage covers). New object: the saturated harmonic period quotient.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from bt861_code_register_is_steinberg import canon
from w33_pass11389_prime_period_family import prime_period_witness

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/w33_pass11389_parabolic_spatial_cover.json"
FCC_CHANGE = sp.Matrix([[-2, 0, 1], [1, 1, -1], [2, -1, 0]])
FCC_ROOTS = sp.Matrix([[1, 0, -1], [-1, 1, -1], [0, -1, 0]])


def geometry():
    points = sorted({canon(x) for x in product(range(3), repeat=4) if any(x)})
    index = {p: i for i, p in enumerate(points)}

    def symp(x, y):
        return (x[0]*y[2]-x[2]*y[0]+x[1]*y[3]-x[3]*y[1]) % 3

    lines = [t for t in combinations(range(40), 4)
             if all(symp(points[i], points[j]) == 0 for i, j in combinations(t, 2))]
    lindex = {frozenset(l): i for i, l in enumerate(lines)}
    edges = [(p, 40+l) for l, ps in enumerate(lines) for p in ps]
    eindex = {e: i for i, e in enumerate(edges)}
    boundary = sp.zeros(80, 160)
    for i, (a, b) in enumerate(edges):
        boundary[a, i], boundary[b, i] = -1, 1

    def permutation(matrix):
        J = np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
        assert np.array_equal((matrix.T @ J @ matrix) % 3, J % 3)
        pp = [index[canon(tuple(int(a) % 3 for a in matrix @ np.array(p)))]
              for p in points]
        return tuple(pp + [40+lindex[frozenset(pp[p] for p in l)] for l in lines])

    assert len(points) == len(lines) == 40 and len(edges) == 160
    return points, lines, edges, eindex, boundary, permutation


def radical(permutation, kind):
    matrices = []
    for a, b, c in product(range(3), repeat=3):
        if kind == "line":
            m = np.eye(4, dtype=int)
            m[:2, 2:] = [[a, b], [b, c]]
        else:
            m = np.array([[1,a,c,b],[0,1,b,0],[0,0,1,0],[0,0,-a,1]])
        matrices.append(permutation(m))
    group = set(matrices)
    assert len(group) == 27
    compose = lambda a, b: tuple(a[b[i]] for i in range(80))
    assert all(compose(a,b) in group for a in group for b in group)
    center = sum(all(compose(a,b) == compose(b,a) for b in group) for a in group)
    assert center == (27 if kind == "line" else 3)
    return matrices, center


def build_cover(kind="line"):
    points, lines, edges, ei, D, permutation = geometry()
    group, center = radical(permutation, kind)
    seen, orbits = set(), []
    for i, (a, b) in enumerate(edges):
        if i not in seen:
            orbit = sorted({ei[(g[a],g[b])] for g in group})
            orbits.append(orbit)
            seen.update(orbit)
    O = sp.zeros(160, len(orbits))
    for j, orbit in enumerate(orbits):
        for i in orbit:
            O[i,j] = 1
    basis = (D*O).nullspace()
    C = O*sp.Matrix.hstack(*basis)
    for j in range(C.cols):
        C[:,j] *= sp.ilcm(*[x.q for x in C[:,j]])
    assert C.cols == 3 and D*C == sp.zeros(80,3)

    adjacency = [[] for _ in range(80)]
    for i, (a,b) in enumerate(edges):
        adjacency[a].append((b,i,1))
        adjacency[b].append((a,i,-1))
    phi, seen, queue, tree = sp.zeros(80,3), {0}, [0], []
    for a in queue:
        for b, i, sign in adjacency[a]:
            if b not in seen:
                seen.add(b)
                queue.append(b)
                tree.append(i)
                phi[b,:] = phi[a,:] + sign*C[i,:]
    assert len(tree) == 79
    periods = C-D.T*phi
    assert all(x.q == 1 for x in periods)
    lattice = hermite_normal_form(periods.T.applyfunc(sp.Integer))
    V, H = periods*lattice.inv().T, C*lattice.inv().T
    assert all(x.q == 1 for x in V)
    assert hermite_normal_form(V.T.applyfunc(sp.Integer)) == sp.eye(3)
    assert V[tree,:] == sp.zeros(79,3)
    assert D*H == sp.zeros(80,3)
    E = H*FCC_CHANGE.inv().T*FCC_ROOTS.T
    positions = phi*lattice.inv().T
    assert H-V == D.T*positions
    K = H.T*H/80
    assert E.T*E/80 == sp.eye(3)*sp.Rational(27,3200)
    assert FCC_CHANGE.det() == -1 or FCC_CHANGE.det() == 1
    assert FCC_CHANGE.T*K.inv()*FCC_CHANGE == FCC_ROOTS.T*FCC_ROOTS*sp.Rational(3200,27)

    vseen, vorbits, vmap = set(), [], {}
    for i in range(80):
        if i not in vseen:
            orb = sorted({g[i] for g in group})
            for v in orb:
                vmap[v] = len(vorbits)
            vorbits.append(orb)
            vseen.update(orb)
    quotient = [(vmap[edges[o[0]][0]], vmap[edges[o[0]][1]], len(o)) for o in orbits]
    assert len(vorbits) == 14 and len(orbits) == 16
    return dict(points=points, lines=lines, edges=edges, D=D, H=H, V=V, E=E,
                K=K, C=C, lattice=lattice, group=group, center=center,
                quotient=quotient, vertex_orbits=vorbits, positions=positions,
                permutation=permutation, ei=ei)


def bloch(edges, voltage, momentum):
    M = 4*np.eye(80, dtype=complex)
    for (a,b), vector in zip(edges, np.asarray(voltage, dtype=float)):
        phase = np.exp(1j*np.dot(vector,momentum))
        M[a,b] -= phase
        M[b,a] -= phase.conjugate()
    return M


def line_levi(cover):
    H, Q = cover["H"], cover["H"].T*cover["H"]
    reps = {}
    for x in product(range(3), repeat=4):
        A = sp.Matrix(2,2,x)
        if A.det() % 3:
            Ai = A.inv().T.applyfunc(lambda v: int(v.p*pow(int(v.q),-1,3)) % 3)
            m = np.zeros((4,4), dtype=int)
            m[:2,:2], m[2:,2:] = np.array(A,dtype=int), np.array(Ai,dtype=int)
            g = cover["permutation"](m)
            gH = sp.zeros(160,3)
            for i,(a,b) in enumerate(cover["edges"]):
                gH[cover["ei"][(g[a],g[b])],:] = H[i,:]
            R = Q.inv()*H.T*gH
            assert H*R == gH and R.T*Q*R == Q
            assert all(v.q == 1 for v in R)
            reps[tuple(R)] = R
    assert len(reps) == 24
    # Norm one proves absolute irreducibility of the actual three-space.
    assert sum(sp.trace(R)**2 for R in reps.values()) == 24
    return list(reps.values())


def quartic(cover):
    D, E, edges = cover["D"], cover["E"], cover["edges"]
    L0 = D*D.T
    Lplus = (L0+sp.ones(80)/80).inv()-sp.ones(80)/80
    values = []
    for direction in [(1,0,0),(1,1,1),(1,1,0),(1,2,3)]:
        t = E*sp.Matrix(direction)
        v2 = sp.zeros(80,1)
        for i,(a,b) in enumerate(edges):
            v2[a] += t[i]**2/2
            v2[b] += t[i]**2/2
        values.append(-sum(x**4 for x in t)/(12*80)-(v2.T*Lplus*v2)[0]/80)
    axis, mixed = values[0], (values[1]-3*values[0])/3
    assert values[2] == 2*axis+mixed
    assert values[3] == axis*(1+16+81)+mixed*(4+9+36)
    return axis, mixed


def matrix_json(M):
    return [[str(x) for x in M.row(i)] for i in range(M.rows)]


def internal_sector(cover):
    """Exact flat-band gauge proof and an actual H27 qutrit intertwiner."""
    orbits, H, edges, D = (cover[k] for k in ("vertex_orbits","H","edges","D"))
    PN = sp.zeros(80)
    orbit_of = {}
    for j,orbit in enumerate(orbits):
        assert all(all(x.q==1 for x in cover["positions"][a,:]-cover["positions"][orbit[0],:])
                   for a in orbit)
        for a in orbit:
            orbit_of[a] = j
            for b in orbit:
                PN[a,b] = sp.Rational(1,len(orbit))
    Q = sp.eye(80)-PN
    assert PN*PN == PN and sp.trace(PN)==14 and sp.trace(Q)==66
    # Nontrivial radical representations vanish on singleton vertex orbits.
    active = {j for j,o in enumerate(orbits) if len(o)>1}
    adjacency = {j:[] for j in active}
    orbit_edges = {}
    for i,(a,b) in enumerate(edges):
        ja,jb = orbit_of[a],orbit_of[b]
        orbit_edges.setdefault((ja,jb),[]).append(i)
    for (a,b),indices in orbit_edges.items():
        if a in active and b in active:
            vector = H[indices[0],:]
            assert all(H[i,:]==vector for i in indices)
            adjacency[a].append((b,vector))
            adjacency[b].append((a,-vector))
    assert len(active)==9 and sum(map(len,adjacency.values()))==16
    start = min(active)
    potentials = {start:sp.zeros(1,3)}
    queue = [start]
    for a in queue:
        for b,vector in adjacency[a]:
            candidate = potentials[a]+vector
            if b not in potentials:
                potentials[b]=candidate
                queue.append(b)
            else:
                assert potentials[b]==candidate
    assert len(potentials)==9
    psi = sp.zeros(80,3)
    for a in range(80):
        if orbit_of[a] in potentials:
            psi[a,:]=potentials[orbit_of[a]]
    for (ja,jb),indices in orbit_edges.items():
        if ja in active and jb in active:
            for i in indices:
                a,b=edges[i]
                assert H[i,:]==psi[b,:]-psi[a,:]
        else:
            # The omitted phase-bearing edge blocks kill every internal state.
            B=sp.zeros(80)
            for i in indices:
                a,b=edges[i]
                B[a,b]=1
            assert B*Q==sp.zeros(80) and B.T*Q==sp.zeros(80)
    L0=D*D.T
    W=L0-4*sp.eye(80)
    assert W*(W*W-6*sp.eye(80))*Q==sp.zeros(80)
    assert sp.trace(W*Q)==0 and sp.trace(W*W*Q)==240
    grading=sp.diag(*([1]*40+[-1]*40))
    assert grading*W+W*grading==sp.zeros(80)
    zero=Q-W*W*Q/6
    assert zero*zero==zero and sp.trace(zero)==26

    group=cover["group"]
    character=[sum(i==g[i] for i in range(80)) for g in group]
    linear=[]
    for u,v in product(range(3),repeat=2):
        coefficients=[0,0,0]
        for i,(a,b,c) in enumerate(product(range(3),repeat=3)):
            coefficients[-(u*a+v*b)%3]+=character[i]
        assert coefficients[1]==coefficients[2]
        multiplicity=sp.Rational(coefficients[0]-coefficients[1],27)
        assert multiplicity==(14 if (u,v)==(0,0) else 3)
        linear.append(dict(character=[u,v],multiplicity=int(multiplicity)))
    assert character[0]==80 and character[1]==character[2]==17
    # The two central-nontrivial irreps have dimension3 and multiplicity7.
    assert sp.Rational(character[0]-character[1],9)==7
    def P(g):
        a=np.zeros((80,80),dtype=complex)
        a[np.array(g),np.arange(80)]=1
        return a
    X,Y,Z=P(group[9]),P(group[3]),P(group[1])
    Zi=sp.Matrix(Z.real.astype(int))
    assert sp.trace(Zi*zero)==-1
    assert sp.trace(grading*zero)==2 and sp.trace(grading*Zi*zero)==11
    # Central-nontrivial zero rank=(26-(-1))/3=9, index=(2-11)/3=-3.
    omega=np.exp(2j*np.pi/3)
    central=(np.eye(80)+omega.conjugate()*Z+omega*Z@Z)/3
    joint=central@(np.eye(80)+Y+Y@Y)/3
    ew,ev=np.linalg.eigh(joint)
    S=ev[:,ew>.5]
    assert S.shape==(80,7)
    T=np.hstack([S,X@S,X@X@S])
    assert np.linalg.norm(T.conj().T@T-np.eye(21))<1e-12
    gamma7=S.conj().T@np.asarray(grading,dtype=float)@S
    gw,gv=np.linalg.eigh(gamma7)
    plus,minus=gv[:,gw>.5],gv[:,gw<-.5]
    assert plus.shape==(7,3) and minus.shape==(7,4)
    a7=4*np.eye(7)-S.conj().T@np.asarray(L0,dtype=float)@S
    B=plus.conj().T@a7@minus
    left,singular,right=np.linalg.svd(B,full_matrices=True)
    assert np.max(abs(singular-np.array([np.sqrt(6),np.sqrt(6),0])))<1e-12
    lifted=B+0.1*np.outer(left[:,-1],right[2,:])
    assert np.max(abs(np.linalg.svd(lifted,compute_uv=False)-np.array([np.sqrt(6),np.sqrt(6),.1])))<1e-12
    samples=[]
    for momentum in ([0,0,0],[.11,.23,-.17],[1.1,-2.3,.87]):
        M=bloch(edges,H,momentum)
        assert all(np.linalg.norm(M[np.ix_(g,g)]-M)<1e-12 for g in group)
        reduced=T.conj().T@M@T
        A=reduced[:7,:7]
        residual=float(np.linalg.norm(reduced-np.kron(np.eye(3),A)))
        assert residual<2e-12
        # Exact forest gauge above proves allmomentum flatness, not this sample.
        gauge=np.diag(np.exp(-1j*np.asarray(psi,dtype=float)@momentum))
        q=np.asarray(Q,dtype=float)
        gauge_residual=float(np.linalg.norm((M-gauge@np.asarray(L0,dtype=float)@gauge.conj().T)@q))
        assert gauge_residual<2e-12
        expected=np.array([4-np.sqrt(6)]*2+[4]*3+[4+np.sqrt(6)]*2)
        assert np.max(abs(np.linalg.eigvalsh(A)-expected))<2e-12
        samples.append(dict(momentum=momentum,qutrit_factor_residual=residual,internal_gauge_residual=gauge_residual))
    return dict(vertex_orbits=orbits,forest_phase_potential=matrix_json(psi),
      invariant_sector_dimension=14,internal_sector_dimension=66,
      exact_internal_spectrum={"4-sqrt(6)":20,"4":26,"4+sqrt(6)":20},
      all_momentum_flat_band_proof=True,linear_characters=linear,
      central_traces=[80,17,17],central_sector_ranks=[38,21,21],
      qutrit_irrep_multiplicity_each=7,qutrit_reduced_spectrum={"4-sqrt(6)":2,"4":3,"4+sqrt(6)":2},
      qutrit_seed_real=S.real.tolist(),qutrit_seed_imag=S.imag.tolist(),
      group_permutations=group,operator_samples=samples,
      completion=dict(operator="L_completed(p)=L_native(p)+(27/12800)*ell_FCC(p)*Q_internal",
        fcc_positive_root_directions=[[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],
        laplacian="ell_FCC=sum_positive_roots(2-2*cos(p.dot(root)))=4|p|^2+O(p^4)",
        status="Explicit additional local-in-space coupling, supplied by a common-acoustic-speed requirement; not implied by native edges or symmetry alone.",
        common_quadratic_speed="27/3200"),
      gapless_alternative=dict(operator="L_gapless=P_N L_native P_N+Q(4I-L_native)^2Q+(27/12800)ell_FCC Q",
        exact_zero_projector="Q-(L0-4I)^2Q/6",internal_zero_rank=26,
        zero_projector_central_traces=[26,-1,-1],zero_projector_grading_traces=[2,11,11],
        central_qutrit_zero_rank_each=9,central_qutrit_point_zero_rank_each=3,
        central_qutrit_line_zero_rank_each=6,central_qutrit_index_each=-3,
        multiplicity_space_index=-1,protected_qutrit_copies_each=1,accidental_extra_zero_copies_each=2,
        rectangular_qutrit_block_shape=[3,4],native_singular_values=[float(x) for x in singular],
        rank_lift_sample=dict(epsilon="1/10",perturbed_singular_values=[float(x) for x in np.linalg.svd(lifted,compute_uv=False)],
          squared_gaps=[0,"1/100","1/100",6,6,6,6],scope="Supplied rankone index-preserving deformation; small epsilon is not predicted."),
        scope="Supplied squared-adjacency wave action. Index protection only within H27-equivariant bipartite first-order operators; generic wave mass counterterms need not preserve it. No physical fermion chirality or radiative mass protection inferred."))


def produce():
    line, point = build_cover("line"), build_cover("point")
    reps = line_levi(line)
    internal=internal_sector(point)
    family=[prime_period_witness(q) for q in (2,3,5)]
    axis, mixed = quartic(line)
    samples = []
    for direction in [(1,0,0),(1,1,1),(1,2,3)]:
        k = np.array(direction)*0.01
        actual = float(np.linalg.eigvalsh(bloch(line["edges"],line["E"],k))[0])
        quadratic = 27/3200*float(k@k)
        fourth = float(axis)*sum(k**4)+float(mixed)*sum(k[i]**2*k[j]**2 for i,j in combinations(range(3),2))
        assert abs(actual-quadratic-fourth) < 2e-11
        samples.append(dict(momentum=list(k),eigenvalue=actual,quadratic=quadratic,quartic=fourth))
    pK = point["K"]
    assert pK == line["K"]
    positions = Counter(tuple(sp.frac(x) for x in line["positions"].row(i)) for i in range(80))
    assert len(positions) == 14
    response = sp.Matrix([[x*x,y*y,z*z,x*y,x*z,y*z]
                          for x,y,z in line["E"].tolist()])
    assert response.rank() == 4
    missing = response.nullspace()
    assert len(missing) == 2
    # Another line: symplectic exchange of the two Lagrangian planes.
    exchange = np.array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    g = line["permutation"](exchange)
    moved = sp.zeros(160,3)
    for i,(a,b) in enumerate(line["edges"]):
        moved[line["ei"][(g[a],g[b])],:] = line["H"][i,:]
    assert line["D"]*moved == sp.zeros(80,3) and moved.T*moved == line["H"].T*line["H"]
    gap = float(np.linalg.eigvalsh(bloch(line["edges"],line["V"],[0,0,0]))[1])
    assert abs(gap-(4-np.sqrt(6))) < 1e-12
    return dict(status="PASS", reservation="4093ef96e", pass_number=11389,
      prior_owners=["bt861_code_register_is_steinberg.py","w33_pass11070_dual_parabolic_flag_diamond.py","PASS5675_5682_external_prior_art.md"],
      scope="Selected parabolic, equal-weight Levi Laplacian, connected periodic cover. Wave dynamics and line selection are supplied; no Einstein equations or physical units derived.",
      native=dict(vertices=80,edges=160,cycle_rank=81),
      line=dict(radical_order=27,center_order=line["center"],fixed_cycle_rank=3,
                vertex_orbit_sizes=[len(o) for o in line["vertex_orbits"]],quotient_edges=line["quotient"],
                harmonic_raw_gram=matrix_json(line["C"].T*line["C"]),raw_period_lattice=matrix_json(line["lattice"]),
                edges=line["edges"],integer_voltage=matrix_json(line["V"]),harmonic_voltage=matrix_json(line["H"]),
                cell_positions_in_primitive_periods=matrix_json(line["positions"]),
                distinct_positions_modulo_deck=14,position_multiplicity_counts=dict(Counter(str(n) for n in positions.values())),
                fcc_harmonic_displacements=matrix_json(line["E"]),diffusion_tensor=matrix_json(line["K"]),
                primitive_fcc_change=matrix_json(FCC_CHANGE),fcc_roots=matrix_json(FCC_ROOTS),
                cartesian_quadratic_coefficient="27/3200",cartesian_quartic_axis=str(axis),cartesian_quartic_mixed=str(mixed),
                optical_gap="4-sqrt(6)",levi_order=24,levi_character_norm=1,
                levi_trace_counts=dict(Counter(str(sp.trace(r)) for r in reps)),levi_cohomology_matrices=[matrix_json(r) for r in reps],
                second_line_harmonic_check=True,bloch_samples=samples),
      metric_response=dict(conductance_tangent_rank=4,
          invisible_symmetric_tensors=[matrix_json(t) for t in missing],
          scope="Envelope derivative at the equal-weight cell: deltaK=sum_e delta_w_e E_e E_e^T/80. Two trace-free diagonal metric directions are absent at first order; no full gravitational metric dynamics follows."),
      point=dict(radical_order=27,center_order=point["center"],fixed_cycle_rank=3,
                 same_diffusion_tensor=True,edges=point["edges"],harmonic_voltage=matrix_json(point["H"]),
                 cell_positions_in_primitive_periods=matrix_json(point["positions"]),
                 fcc_harmonic_displacements=matrix_json(point["E"]),
                 integer_voltage=matrix_json(point["V"])),
      internal_h27=internal,small_prime_period_family=family,
      all_field_rank_law="For a selected line in W(3,q), the unipotent orbit quotient has q+1 length4 paths and hence q harmonic cycle directions. Stronger integral A_q periods checked here only at primes2,3,5.",
      global_symmetry_boundary="A nontrivial free-abelian regular Levi cover whose kernel is invariant under full PSp(4,3) has deck rank81, by prior irreducible Steinberg H1. The rank3 period kernel preserves a selected parabolic, not full G.")


if __name__ == "__main__":
    result = produce()
    OUT.write_text(json.dumps(result,indent=2)+"\n")
    print(result["status"], "rank3 connected cover; FCC periods; isotropic acoustic coefficient27/3200")
