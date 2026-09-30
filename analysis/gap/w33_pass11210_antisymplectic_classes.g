# Pass 11192: conjugacy classes of the ANTI-SYMPLECTIC coset of PGSp(2n,3) = PSp(2n,3).2 (n = 2, 3; for n = 2 this is
# W(E6)), i.e. the projective classes of time reversals (antiunitary Clifford operations mod Paulis).  Adapted basis as
# in Pass 11188; tau = diag(1,-1,...) is complex conjugation.  Per outer class: order, size, M = B g B^-1 (row action).
Run := function(n)
  local F, S, form, V, cands, om, basis, i, found, u, v, B, Binv, tau, G, ptsV, act, Gp, Hp, cls, c, g, m, M;
  F := GF(3); S := Sp(2*n, 3);
  form := InvariantBilinearForm(S).matrix;
  V := F^(2*n);
  cands := Filtered(Elements(V), x -> not IsZero(x));
  om := function(a, b) return a*form*b; end;
  basis := [];
  for i in [1..n] do
    found := false;
    for u in cands do
      if found then break; fi;
      if ForAll(basis, b -> IsZero(om(u, b))) then
        for v in cands do
          if ForAll(basis, b -> IsZero(om(v, b))) and om(u, v) = One(F) then
            Append(basis, [u, v]); found := true; break;
          fi;
        od;
      fi;
    od;
  od;
  B := basis; Binv := B^-1;
  tau := DiagonalMat(List([1..2*n], k -> (-1)^(k+1)) * One(F));
  tau := Binv * tau * B;                      # complex conjugation, in GAP's coordinates
  G := Group(Concatenation(GeneratorsOfGroup(S), [tau]));
  ptsV := Set(List(cands, NormedRowVector));
  act := ActionHomomorphism(G, ptsV, OnLines);
  Gp := Image(act); Hp := Image(act, S);
  cls := Filtered(ConjugacyClasses(Gp), c -> not Representative(c) in Hp);
  Print("n=", n, " |PGSp|=", Size(Gp), " |PSp|=", Size(Hp), " outer classes=", Length(cls), "\n");
  for c in cls do
    g := Representative(c);
    m := PreImagesRepresentative(act, g);
    M := B * m * Binv;
    Print("OUTER n=", n, " order=", Order(g), " size=", Size(c), " M=", List(M, r -> List(r, IntFFE)), "\n");
  od;
end;
Run(2);
Run(3);
QUIT;
