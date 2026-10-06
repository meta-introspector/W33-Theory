# Pass 11540 reservation — history VO(3,3) affine-polar bridge

> **Post-computation correction (Pass 11548):** this reservation records the
> intended orthogonal-group attack, but its anticipated identification of the
> repository's 648-element PSp subgroup with affine SO was not objectwise
> correct.  The final corrected result is: repo PSp linear =
> `ker(sigma)=W(D3)`, while `SO(3,3)=ker(det)` is a distinct order-24
> subgroup; their intersection is `A4=Omega(3,3)`.  See
> `PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.md`.

Reserved 2026-10-05 by the GitHub research continuation requested in chat.

Original scope:
- identify the 27-event null-history graph as the parabolic affine orthogonal polar graph VO(3,3);
- derive its adjacency spectrum by the F3^3 Fourier characters and match spectral multiplicities to quadratic norm shells;
- identify the natural affine orthogonal orders 27*48=1296 and 27*24=648;
- audit the relation of the existing oriented-history kernel 3^3:A4 (order 324) to Omega(3,3) and the finite spinor-norm quotient;
- preserve the point/line firewall and make no continuum/gravity claim.

Primary existing owners cross-checked:
- papers/forty_points/sec04_time.tex
- analysis/w33_20260924_history_bigcell_q43_compactification.py
- analysis/w33_20260924_history_invariant_cycle_orientation.py
- analysis/w33_20260924_bell_shell_quadratic_history_intertwiner.py
- data/PART_W33_PASS9741_9748_ORIENTATION_CHARACTER_WELD.json
