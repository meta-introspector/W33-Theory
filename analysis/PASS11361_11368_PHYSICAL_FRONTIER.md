# Passes 11361–11368: five physical targets and three exact cross-checks

Producer: analysis/w33_pass11361_11368_physical_frontier.py. Certificate: data/w33_pass11361_11368_physical_frontier.json. These are scoped calculations, not a solution of the TOE. The packet reuses the prior qutrit rays of Pass 9949–9956, the actual W33 80-site/160-edge cycle basis, the live E6 quotient chart, and the backreacted axion-winding wall. The current paper keeps masses, mixing, dynamical space-time and the cosmological constant **OPEN**.

## 11361 — a three-ray Bargmann phase becomes a full-rank flavor CP invariant

Let the normalized, linearly independent rays be columns of $X$, set $G=X^\dagger X$, and let $P_i=|u_i\rangle\langle u_i|$. For positive real weights define $A=\sum_i a_iP_i$ and $B=\sum_i b_iP_i$. They are positive full-rank Hermitian Grams. Explicit Yukawa matrices are $Y_u=\operatorname{diag}(\sqrt{a_i})X^\dagger$ and $Y_d=\operatorname{diag}(\sqrt{b_i})X^\dagger$. With $c_{ij}=a_ib_j-a_jb_i$,

$$
\frac{\operatorname{Tr}[A,B]^3}{i}
=-6\det(G)\,\operatorname{Im}\operatorname{Tr}(P_1P_2P_3)\,c_{12}c_{13}c_{23}.
$$

Proof: $[A,B]=XKX^\dagger$ with $K_{ij}=c_{ij}G_{ij}$, so $\det[A,B]=\det(G)\det K$. The zero-diagonal anti-Hermitian $K$ has $\det K=-2i c_{12}c_{13}c_{23}\operatorname{Im}(G_{12}G_{23}G_{31})$. Both $K$ and $[A,B]$ are traceless $3\times3$ matrices, and $\operatorname{Tr}[A,B]^3=3\det[A,B]$. Hence the cubic is nonzero **if and only if** the rays span $\mathbb C^3$, their Bargmann loop has nonzero imaginary part, and all three pairwise weight minors are nonzero.

The **previously committed** optical rays $u=(1,0,0)$, $v=(1,1,0)/\sqrt2$, $w=(1,i,1)/\sqrt3$ give $\operatorname{Tr}(P_1P_2P_3)=(1+i)/6$ and $\det G=1/6$. Illustrative supplied weights $a=(1,2,4)$, $b=(2,3,5)$ give $c=(-1,-3,-2)$ and $\operatorname{Tr}[A,B]^3/i=1$ exactly. Both Grams have three distinct positive eigenvalues; the mixing invariant is about $8.75\times10^{-5}$. Conjugation reverses its sign. This supplies the algebraic map missing in 11347, while W33 dynamics does not yet select the weights, orientation, masses, or observed CKM data. The physical commutator interpretation is standard [Jarlskog 1985](https://cds.cern.ch/record/161795); the ray phase was already owned by Pass 9949–9956. Claude-track Passes 11355–11360 study a distinct degree-six Clifford-twirled invariant of tick time-reversal. The shared cubic degree does not provide a map from those ticks to the flavor Grams constructed here.

## 11362 — one real 1PI normal entry independent of Ward columns

On the actual quotient, the family orbit has rank eight in 22 real coordinates. In one **live** chart the invariant $I=\phi^\dagger\phi$ has gradient $v$, with $\|v\|^2=2$ and $T^Tv=0$ for every family tangent. Add a supplied real heavy scalar with $m_S^2=M^2+g(I-I_0)$. Its one-loop potential is $V_1=m_S^4[\log(m_S^2/\mu^2)-3/2]/(64\pi^2)$. At the named scheme point $\mu^2=M^2/e$, the tadpole vanishes and the exact normal Hessian is $H_1=g^2vv^T/(32\pi^2)$: rank one, nonzero eigenvalue $g^2/(16\pi^2)$, and $H_1T=0$. The computed Ward residual for $g=1$ is $2.8\times10^{-19}$.

This realizes one of the 105 free normal entries of 11344 as an invariant loop, rather than merely adding a matrix by hand. Its value depends on supplied matter and matching. The native E6 heavy thresholds and the remaining 104 directions are still uncomputed. A stored tangent matrix from a different singular-vector chart cannot be paired with a newly computed gradient without basis transport.

## 11363–11364 — hidden equal-lapse null and exact generic metric-cycle rank

Let $B$ be oriented incidence, $U=|B|$, $C$ an integral cycle matrix with $BC=0$, and $L=U^TN$. The named pair-metric model has
$H_N=-GA^{-1}G^T$, where $G=U\operatorname{diag}(j/\sqrt{1+j^2})C$ and $A=C^T\operatorname{diag}(L/(1+j^2)^{3/2})C>0$.
Choose a vertex potential $\eta$, set $r=B^T\eta$, require $|r_e|<L_e$, and choose $j_e=(r_e/L_e)/\sqrt{1-(r_e/L_e)^2}$. Then $C^TLj/\sqrt{1+j^2}=C^TB^T\eta=0$: this is an **actual stationary point**.

At equal lapses, $\eta$ is a third exact left null in addition to the all-ones and point/line checkerboard vectors. On each edge $(U^T\eta)_e(B^T\eta)_e=(B^T\eta^2)_e$, so $G^T\eta=0$. This explains the old rank-77 equal-lapse numerical control. At a deterministic nonconstant positive rational lapse/potential point, exact elimination over $\mathbb F_{1000003}$ gives rank 78. Two exact independent nulls bound the rational rank above by 78, so it is exactly 78. A nonzero rational minor makes rank 78 generic in this algebraic stationary family. This strengthens the prior numerical result for the named pair-**metric** cycle action; it is not a Dirac count or a refutation of the distinct vertex-frame multivielbein model. The distinction is consonant with [multigravity-cycle analyses](https://arxiv.org/abs/1410.7774).

## 11365–11366 — wound-wall continuation, positive family, and scale boundary

At the exact zero-bare-$\rho$ wall $h=16/45$, $q=4/3$, the mixed Maxwell/metric odd zero of 11346 remains exact. A restricted odd $m=0$ Schur/Ritz calculation gives approximately $0,10.1528027824,29.0427629324,56.3968737462$; bases 12 and 18 agree to $1.1\times10^{-12}$. These are **trial-space upper bounds**. They do not certify the full wall's negative-mode count. Breathing, wall displacement, axion and tensor sectors, plus nonlinear zero-mode integrability, remain open.

The supplied quantized family has $1\le q\le4/3$, $h=4(q^2-q)/5$, $c=q/5$, and $\rho=q(4-3q)/10$. Thus $q=4/3$ is its unique positive-$q$ zero of the **bare** cosmological parameter. The exact cap function is $F=q(b-1)H/(30b^2)$, with $H$ affine in $q$, $H(1)=15-b-b^2-b^3>0$, and $H(4/3)=20-8b>0$ for $1<b\le5/4$. The whole interval has a positive cap. Yet at $\rho=0$, the Ricci scalar is $R=32/(45b^2)>0$ due to supplied axion stress. No physical length is fixed: $R_{\rm phys}\propto L^{-2}$, and energy density scales as $L^{-4}$. Zero bare $\rho$ therefore means neither flat spacetime nor a measured cosmological constant. Scale invariance by itself is also insufficient to protect vacuum energy in this [primary analysis](https://arxiv.org/abs/2004.01868).

## 11367 — transverse-boost flatness on a native octagon

Every column of the native integral cycle basis is an eight-edge loop. Transverse boost generators $K_x=\sigma_x/2$ and $K_y=\sigma_y/2$ obey $[K_x,K_y]=i\sigma_z/2$. Giving four successive independent edge transports $e^{tK_x},e^{tK_y},e^{-tK_x},e^{-tK_y}$ and identity on the other edges creates rotation holonomy $I+t^2i\sigma_z/2+O(t^3)$. Vertex-frame edge transports $\Lambda_i^{-1}\Lambda_j$ telescope to identity around every loop. Independent transverse edge boosts therefore require compensating rotations before they can represent vertex frames. This is a nonlinear integrability obligation beyond 11345's rank-237 linearization, not full secondary-constraint closure or ghost freedom.

## 11368 — correct a public physical-claim conflict

The legacy site called $\log_{10}(\Lambda_{CC}/M_{Pl}^4)=-(vq+\mu-\lambda)=-127$ an exact resolution. With its own values $v=40$, $q=3$, $\mu=4$, $\lambda=2$, the integer is **122**, not 127. The same catalogue attached GeV units to dimensionless 246 and 125 matches without a mass map. The main site and 404 fallback now identify these as exploratory arithmetic, correct explicit −127 entries, and state that the current paper keeps physical masses and vacuum energy open. Other older catalogue rows remain historical proposals, not validated physics.

All seven certificate sections (covering eight reserved pass numbers) report PASS. The exact results are the CP factorization, graph-rank/null certificates, wall-family factorization and curvature identity, and transverse holonomy. The one-loop contribution is for specified extra heavy matter. Wall Ritz levels are restricted numerical controls. The five full physical targets remain open: complete wall stability, native hard 1PI matching, dynamical CKM selection, generic gravity constraints, and absolute mass/vacuum-energy determination.

Validation: all seven producer sections passed on the final source/site state, as did five direct independent regression functions and Python source compilation. The isolated pytest launcher timed out before collection and is not counted as a test pass. Tectonic rebuilt the combined PDF (895447 bytes), whose extracted text contains the new ledger and theorem. Corpus intake and publication are recorded separately in the session note.
