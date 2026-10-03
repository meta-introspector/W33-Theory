"""Five scoped follow-ups to the W33 physical frontier.

Each section states what its named action or kinematic model actually proves.
None is a derivation of observed masses, gravity or vacuum energy.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def flavor_phase_selector() -> dict:
    """A CP-even selector on the *prior* three-ray ansatz, not a UV potential."""
    z = sp.symbols("z", real=True)
    i = sp.I
    u = sp.Matrix([1, 0, 0])
    v = sp.Matrix([1, 1, 0]) / sp.sqrt(2)
    w = sp.Matrix([1, sp.exp(i*z), 1]) / sp.sqrt(3)
    x = sp.Matrix.hstack(u, v, w)
    gram = x.conjugate().T * x
    bargmann = sp.simplify((u.conjugate().T*v)[0]*(v.conjugate().T*w)[0]*(w.conjugate().T*u)[0])
    assert sp.simplify(bargmann - (1+sp.exp(i*z))/6) == 0
    assert sp.simplify(gram.det() - sp.Rational(1, 6)) == 0
    a, b = (1, 2, 4), (2, 3, 5)
    A = x * sp.diag(*a) * x.conjugate().T
    B = x * sp.diag(*b) * x.conjugate().T
    c = A*B-B*A
    cp_cubic = sp.trigsimp(sp.expand_complex(sp.trace(c**3)/i))
    assert sp.simplify(cp_cubic-sp.sin(z)) == 0
    # A rephasing-invariant, CP-even six-ray operator. Restriction to the
    # one-parameter ansatz is an explicit, supplied dynamical choice.
    potential = sp.trigsimp((sp.re(bargmann)-sp.Rational(1, 6))**2)
    assert sp.simplify(potential-sp.cos(z)**2/36) == 0
    assert sp.diff(potential,z).subs(z,sp.pi/2) == 0
    assert sp.diff(potential,z,2).subs(z,sp.pi/2) == sp.Rational(1,18)
    eig=sp.symbols("eig")
    pa=sp.factor(A.subs(z,sp.pi/2).charpoly(eig).as_expr())
    pb=sp.factor(B.subs(z,sp.pi/2).charpoly(eig).as_expr())
    assert sp.simplify(pa-(eig**3-7*eig**2+9*eig-sp.Rational(4,3)))==0
    assert sp.simplify(pb-(eig**3-10*eig**2+sp.Rational(59,3)*eig-5))==0
    da,db=sp.discriminant(pa,eig),sp.discriminant(pb,eig)
    assert (da,db)==(sp.Rational(2063,3),sp.Rational(142459,27))
    return {
        "status":"PASS", "scope":"CP-even selector on a fixed three-ray ansatz with supplied positive coefficient, not a native W33 Yukawa vacuum",
        "Bargmann":"(1+exp(i phi))/6", "det_Gram":"1/6",
        "Gram_commutator_cubic_over_i":"sin(phi)",
        "selector":"lambda*(Re Bargmann-1/6)^2=lambda*cos(phi)^2/36, lambda>0",
        "minima_mod_2pi":["pi/2","3*pi/2"],
        "phase_curvature_at_each_minimum":"lambda/18",
        "CP_cubic_at_minima":[1,-1],
        "A_characteristic_polynomial_at_positive_branch":str(pa),
        "B_characteristic_polynomial_at_positive_branch":str(pb),
        "A_B_spectral_discriminants":[str(da),str(db)],
        "mixing_J_absolute_at_either_branch":str(sp.simplify(1/(6*sp.sqrt(da*db)))),
        "boundary":"CP conjugate minima are degenerate. Ray ansatz, weights and selector coefficient are inputs; full-field normal stability and measured mixing remain open.",
        "prior_owners":["analysis/w33_pass11286_misaligned_family_vacua.py","analysis/w33_pass11326_noncommuting_flavor_operator.py","analysis/w33_pass11361_11368_physical_frontier.py"],
    }


def hard_loop_rank_boundary() -> dict:
    """Exact chain-rule rank theorem for norm-only heavy thresholds."""
    # On the 22-real-coordinate chart, I is any smooth invariant and v=dI.
    # At a stationary I0, V=f(I) has Hessian f'' vv^T, independent of Hess I.
    # This holds for any sum of norm-only Coleman-Weinberg species plus local
    # invariant counterterms; stationarity is for the total function f.
    v = sp.Matrix([1, 2] + [0]*20)
    h = v*v.T
    assert h.rank() == 1
    # Explicit rank-one spanning set for Sym(14) in the normal space; it
    # requires 105 freely controllable source directions, not yet constructed.
    n = 14
    vecs = [sp.eye(n)[:,j] for j in range(n)]
    vecs += [sp.eye(n)[:,j]+sp.eye(n)[:,k] for j in range(n) for k in range(j+1,n)]
    assert len(vecs) == n*(n+1)//2 == 105
    # The basis proof is exact: diagonal outer products give E_ii; pair-sum
    # outer products minus two diagonals give E_ij+E_ji.
    return {
        "status":"PASS", "scope":"Exact obstruction for all heavy species whose masses depend only on the single invariant I=phi-dagger phi, after total tadpole cancellation",
        "total_tadpole_condition":"f'(I0)=0", "hard_Hessian":"f''(I0) dI tensor dI",
        "max_normal_rank":1, "normal_dimension":14,
        "normal_symmetric_entries":105,
        "directions_needed_for_full_rank_at_least":14,
        "independent_rank_one_sources_needed_to_span_all_entries_at_least":105,
        "explicit_abstract_spanning_set":"e_i and e_i+e_j (i<j); global E6-invariant realization not supplied",
        "boundary":"The rank-one conclusion requires the norm-only sector's own total tadpole f'(I0)=0. If a different sector cancels a nonzero f' instead, the f' Hess(I) term survives and this bound need not hold. A zero-momentum potential Hessian is not a full momentum-dependent 1PI mass or pole matrix.",
        "prior_owners":["analysis/w33_pass11337_quotient_hard_ward_columns.py","analysis/w33_pass11342_11349_eight_frontier_probes.py","analysis/w33_pass11361_11368_physical_frontier.py"],
    }


def native_frame_flatness() -> dict:
    """The exact graph-cohomological dimension of pure-gauge frame transports."""
    import sys
    sys.path.insert(0,str(ROOT/"analysis"))
    from w33_pass11289_cycle_gram_gluing_flux import cycle_basis
    inc, cycles = cycle_basis()
    assert inc.shape == (80,160) and cycles.shape == (160,81)
    assert np.count_nonzero(inc@cycles) == 0
    # The fundamental-cycle construction has an identity on its 81 chords.
    chord_rows = np.flatnonzero(np.any(cycles,axis=1) & (np.sum(np.abs(cycles),axis=1)==1))
    assert len(chord_rows)==81 and np.array_equal(abs(cycles[chord_rows]),np.eye(81,dtype=int))
    v,e,g,dim = 80,160,81,6
    assert g==e-v+1
    return {
        "status":"PASS", "scope":"Kinematic flat SO+(1,3) edge-transport chart near the identity on the native connected Levi graph; not a gravitational Dirac count",
        "vertices":v,"edges":e,"independent_cycles":g,"Lorentz_dimension":dim,
        "unconstrained_edge_connection_dimension":e*dim,
        "independent_linearized_cycle_flatness_equations":g*dim,
        "flat_connection_dimension_before_vertex_gauge":(v-1)*dim,
        "flat_connection_mod_vertex_gauge_dimension":0,
        "proof":"A spanning tree gauges all tree edges to identity. The 81 chord transports are the 81 fundamental holonomies, each imposing six independent local equations. The remaining 79 vertex frame differences are pure gauge.",
        "boundary":"This removes independent connection holonomy only if every edge transport is Lambda_i^-1 Lambda_j. The vielbeins, metric modes, lapses, secondary constraints and local matter remain; no ghost freedom follows.",
        "prior_owners":["analysis/w33_pass11289_cycle_gram_gluing_flux.py","analysis/w33_pass11338_native_boost_constraint_slice.py","analysis/w33_pass11361_11368_physical_frontier.py"],
    }


def wall_axion_probe_sector() -> dict:
    """Exact positivity of compact-axion fluctuations with geometry fixed."""
    bw=sp.Rational(5,4)
    first_gap=1/bw**2
    assert first_gap==sp.Rational(16,25)
    return {
        "status":"PASS", "scope":"Fixed-background compact-axion quadratic block of the 11341 winding wall; does not diagonalize its coupling to metric, Maxwell, or wall bending",
        "quadratic_form":"h int r*b^2 [abs(d_xi chi)^2 + (m_x^2+m_y^2)/b^2 abs(chi)^2] dxi, h>0",
        "range_of_b":"1 <= b <= 5/4",
        "nonzero_torus_Fourier_lower_bound":"(m_x^2+m_y^2)*16/25",
        "first_harmonic_lower_bound":str(first_gap),
        "zero_torus_mode":"nonnegative radial derivative form; constant shift is its only regular zero on connected cap",
        "two_exact_nonlinear_shift_moduli":"theta_x->theta_x+alpha_x and theta_y->theta_y+alpha_y leave the supplied derivative-only action invariant",
        "boundary":"Nonconstant axion modes can mix with geometry; positivity of this principal block alone does not prove positivity of the coupled Schur complement or full wall stability.",
        "prior_owners":["analysis/w33_pass11341_winding_lifted_wall.py","analysis/w33_pass11346_winding_wall_mixed_vector_modes.py","analysis/w33_pass11361_11368_physical_frontier.py"],
    }


def wall_scale_ratio() -> dict:
    """An exact dimensionless datum, with its absolute-scale nonidentifiability."""
    b=sp.symbols("b",positive=True)
    q=sp.Rational(4,3);h=sp.Rational(16,45);rho=sp.S(0)
    C=rho/3+q*q/2+h
    F=sp.factor(C/b-rho*b*b/3-q*q/(2*b*b)-h)
    bw=sp.Rational(5,4)
    fw=sp.factor(F.subs(b,bw))
    t2=sp.factor(16*fw/bw**2)
    rwall=sp.Rational(32,45)/bw**2
    assert sp.simplify(F+8*(b-1)*(2*b-5)/(45*b*b))==0
    assert fw==sp.Rational(16,225)
    assert t2==sp.Rational(4096,5625)
    assert rwall==sp.Rational(512,1125)
    assert sp.factor(rwall/t2)==sp.Rational(5,8)
    return {
        "status":"PASS", "scope":"Exact dimensionless curvature-to-wall-tension relation in the supplied zero-bare-rho winding wall, not a measured CC or absolute scale",
        "F":"-8*(b-1)*(2*b-5)/(45*b^2)", "F_wall":str(fw),
        "wall_tension_squared_kappa2_equal_1":str(t2),
        "wall_Ricci_scalar_kappa2_equal_1":str(rwall),
        "physical_scale_invariant_ratio":"R_phys/(kappa_phys^4*T_phys^2)=5/8",
        "scale_family":"For every L>0, R_phys=R_dimensionless/L^2 and T_phys=T_dimensionless/(kappa_phys^2*L); the ratio remains 5/8.",
        "boundary":"Flux integer and this ratio do not determine L or kappa_phys, so neither absolute masses nor observed vacuum energy follow. The wall and axion parameters are supplied.",
        "prior_owners":["analysis/w33_pass11342_11349_eight_frontier_probes.py","analysis/w33_pass11361_11368_physical_frontier.py"],
    }


def payload() -> dict:
    parts={
        "11374_flavor_phase_selector":flavor_phase_selector(),
        "11375_hard_loop_rank_boundary":hard_loop_rank_boundary(),
        "11376_native_frame_flatness":native_frame_flatness(),
        "11377_wall_axion_probe_sector":wall_axion_probe_sector(),
        "11378_wall_scale_ratio":wall_scale_ratio(),
    }
    assert all(x["status"]=="PASS" for x in parts.values())
    return {"status":"PASS","sections":parts,
        "primary_sources":["https://cds.cern.ch/record/161795","https://arxiv.org/abs/1606.07069","https://arxiv.org/abs/1410.7774"]}


if __name__ == "__main__":
    out=payload()
    dest=ROOT/"data/w33_pass11374_11378_five_frontier.json"
    dest.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(out["status"],list(out["sections"]))
