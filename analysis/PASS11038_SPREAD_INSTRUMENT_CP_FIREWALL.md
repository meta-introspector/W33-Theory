# Pass 11038 — W33 spreads are tomography frames, not ten raw channel contexts

Producer: `analysis/w33_pass11038_spread_instrument_cp_firewall.py`
Certificate: `data/w33_pass11038_spread_instrument_cp_firewall.json`

The exact W33 spread census is 36. Exactly 9 spreads contain the Bell line. Spread pairs meet in either 1 or 4 lines, with histogram `1:360, 4:270`.

The uploaded temporal note suggested using a spread as a ten-setting instrument family. The geometry is right, but complete positivity adds a sharp firewall.

For the two-sided Weyl basis map
`Phi_(u,v)(A)=P_u A P_v^dagger`,
the Choi operator is rank one, `|P_u>><<P_v|`. It is positive only on the diagonal `u=v`, i.e. the four projective points of the Bell line.

Therefore every spread contains exactly 4 CP basis rays and 36 non-CP basis rays.

The 9 Bell-containing spreads have one four-CP context plus nine zero-CP contexts. The other 27 spreads have four one-CP contexts plus six zero-CP contexts.

This does not make the off-Bell directions unphysical as tomography coordinates. It means a laboratory instrument must realize them through a CP completion, dilation, or positive linear reconstruction rather than treating the raw rank-one superoperators as channels.
