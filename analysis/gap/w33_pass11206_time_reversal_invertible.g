# Pass 11206: which Clifford ticks are inverted by a time reversal?  tau = (x,z) -> (x,-z) on an adapted basis is
# anti-symplectic; a tick g (projective, in PSp(2n,3)) is T-symmetric iff some element of the outer coset tau*PSp
# conjugates g to g^-1, i.e. iff g^tau is PSp-conjugate to g^-1.  Also records whether g is real inside PSp.
Run := function(n)
  local G, F, form, V, cands, om, planes, basis, i, found, u, v, B, D, t, E, ptsV, act, Ep, Gp, tp, cls, c, g, tinv,
        real, totT, totR, totBoth;
  G := Sp(2*n, 3); F := GF(3);
  form := InvariantBilinearForm(G).matrix;
  V := F^(2*n);
  cands := Filtered(Elements(V), x -> not IsZero(x));
  om := function(a, b) return a*form*b; end;
  planes := []; basis := [];
  for i in [1..n] do
    found := false;
    for u in cands do
      if found then break; fi;
      if ForAll(basis, b -> IsZero(om(u, b))) then
        for v in cands do
          if ForAll(basis, b -> IsZero(om(v, b))) and not IsZero(om(u, v)) then
            if om(u, v) <> One(F) then v := -v; fi;
            Append(basis, [u, v]); found := true; break;
          fi;
        od;
      fi;
    od;
  od;
  B := basis;
  D := DiagonalMat(List([1..2*n], k -> (-1)^(k+1) * One(F)));
  t := B^-1 * D * B;
  if t * form * TransposedMat(t) <> -form then Error("tau not anti-symplectic"); fi;
  E := Group(Concatenation(GeneratorsOfGroup(G), [t]));
  ptsV := Set(List(cands, NormedRowVector));
  act := ActionHomomorphism(E, ptsV, OnLines);
  Ep := Image(act);
  Gp := Image(act, G);
  tp := Image(act, t);
  cls := ConjugacyClasses(Gp);
  Print("n=", n, " |PSp|=", Size(Gp), " |ext|=", Size(Ep), " classes=", Length(cls), "\n");
  totT := 0; totR := 0; totBoth := 0;
  for c in cls do
    g := Representative(c);
    tinv := IsConjugate(Gp, g^tp, g^-1);
    real := IsConjugate(Gp, g, g^-1);
    if tinv then totT := totT + Size(c); fi;
    if real then totR := totR + Size(c); fi;
    Print("CLASS n=", n, " order=", Order(g), " size=", Size(c), " T_invertible=", tinv, " real=", real, "\n");
  od;
  Print("n=", n, " T-invertible elements=", totT, " fraction=", totT / Size(Gp), " real elements=", totR,
        " fraction=", totR / Size(Gp), "\n");
end;
Run(2);
Run(3);
QUIT;
