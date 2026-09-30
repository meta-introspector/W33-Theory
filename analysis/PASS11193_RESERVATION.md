# Pass 11193 — every two-qutrit Clifford from one fixed W33 perfect gate, optimally in depth at most two

Producer: `analysis/w33_pass11193_optimal_perfect_gate_compiler.py`  
Certificate: `data/w33_pass11193_optimal_perfect_gate_compiler.json`  
Regression: `tests/test_w33_pass11193_optimal_perfect_gate_compiler.py`

Pass 11177 identified the local, perfect, and partial two-qutrit Clifford
classes with the three relations of the rank-three action on the 45 W33
factorisations. This pass turns that classification into an exact compiler.

## One fixed entangler is enough

Let (H) be the stabiliser of the standard tensor factorisation. It has order
(1152) and exact local ABI

```text
H = (SL(2,3) x SL(2,3)) : C2,
```

where the final (C_2) swaps the qutrits. Let (p) be the symplectic
tetracode representative obtained from the Pass 11156 time-reversing
similitude by local time reversal:

```text
p = [[1,0,1,0],
     [0,2,0,2],
     [1,0,2,0],
     [0,2,0,1]].
```

It has order four and a four-qutrit `AME(4,3)` Choi state. Exhaustive double
coset construction gives

```text
Sp(4,3) = H  disjoint_union  H p H  disjoint_union  H p h_* p H,
```

for the single fixed local word

```text
h_* = [[0,0,0,1],
       [0,0,2,0],
       [0,1,0,0],
       [2,1,0,0]].
```

All 51,840 matrices are reconstructed objectwise. The three disjoint pieces
have sizes

| optimal number of (p) gates | class | matrices |
|---:|---|---:|
| 0 | local or swap, (H) | 1,152 |
| 1 | perfect, (HpH) | 13,824 |
| 2 | partial, (Hph_*pH) | 36,864 |

The maximum is therefore two and the mean is (76/45). Optimality is not a
search claim: zero uses exactly (H), at most one uses exactly
(H\mathbin{\cup}HpH), and the third double coset is disjoint from both.
Because (-I\in H), the projective counts are exactly half.

## The rank-three graph is the compiler algebra

If (A) is adjacency in the `SRG(45,12,3,3)` perfect-relation graph, the
objectwise matrix identity is

```text
A^2 = 12 I + 3 A + 3 (J-I-A).
```

Its spectrum is `12^1, 3^20, (-3)^24`. Thus two perfect moves provide 12
returns to the same split and exactly three paths to every perfect or partial
target. A uniformly random middle local factor in `p H p` lands in the three
classes with exact probabilities

```text
local 1/12, perfect 1/4, partial 2/3.
```

This is the Hecke multiplication table made executable. The graph diameter
two is precisely the compiler depth bound.

## An exact walk on the E6 cubic's two Yukawa sectors

Pass 11185 identifies the same 45 factorisations with the 45 monomials of
the E6 cubic. Relative to one complete frame they split into five
`1.10.10` terms and forty `16.16.10` terms. The perfect-gate walk is exactly
lumpable on this partition:

```text
                         next 1.10.10   next 16.16.10
current 1.10.10               4              8
current 16.16.10              1             11
```

After division by the valency 12, the two-state transition matrix is

```text
P = [[1/3,  2/3],
     [1/12, 11/12]].
```

It has stationary distribution `(1/9,8/9)` and transient eigenvalue `1/4`.
Thus a walk begun on a singlet-sector cubic term returns to that sector after
`k` perfect steps with exact probability

```text
1/9 + (8/9)(1/4)^k.
```

This is a concrete, finite dynamics on the coupling *labels*, induced by
uniformly random local dressing of the fixed perfect interaction. It does not
assign Yukawa amplitudes and therefore does not predict a mass or mixing
angle. Its value for the TOE programme is sharper: it isolates the missing
datum. Geometry fixes the transition algebra and relaxation factor; physics
must still supply a state, weighting or action that turns label motion into
coupling values.

## Explicit SUM word and universal-computation boundary

The standard qutrit SUM symplectic matrix lies in the partial class, so its
minimal perfect-gate depth is two. The certificate freezes explicit local
matrices (L,h_*,R\in H) and verifies

```text
SUM = L p h_* p R.
```

Every local matrix is separately decoded as two one-qutrit `SL(2,3)`
operations and an optional SWAP.

The theorem closes the complete two-qutrit **Clifford** compiler using one
fixed perfect interaction. It does not make a Clifford gate universal. When
the already-certified Pass 10944 signed E6 cubic tick is adjoined, its exact
zero-controlled qutrit X plus Fourier inherits the odd-prime approximate
universality theorem of Roy, van de Wetering and Yeh (arXiv:2307.10095). No
Hamiltonian strength, fault rate, or physical implementation is inferred.

Ian Tan's arXiv:2601.19677 supplies external context for the LU uniqueness and
local symmetry of `AME(4,3)`. The repo-owned increments here are the W33
double-coset normal form, its exhaustive word compiler, the optimal depth
census, the Hecke fusion law, and the explicit SUM word.

The fixed gate's tetracode provenance is prior and cited rather than reclaimed:
`PASS10946_CLOCK_CODE_CONE_OBJECTWISE.md`,
`PASS10954_REGULAR_C8_CLOCK_COMPLETION.md`,
`PASS10955_D4_HALFSPIN_CLOCK_BRIDGE.md`, and
`PASS10970_PROJECTIVE_CLOCK_CODE_LATTICE_TOWER.md`. Pass 11091 already
separates the tetracode orientation character, while
`BT927_e8_lift_artifact_reconciliation.md` and
`PASS9173_9196_RESERVATION.md` concern different E8 glue-lift questions. The
raw `5+40` sector counts used by the new walk are explicitly Pass 11185's
result; Pass 11193 adds only the lumped dynamics and compiler interpretation.
Pass 11035's cocycle descent and Pass 11164's three-qutrit F9 tetracode gate
were also checked; neither contains this two-qutrit compiler or sector walk.
