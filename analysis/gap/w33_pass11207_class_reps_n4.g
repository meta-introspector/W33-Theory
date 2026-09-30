# Pass 11188 (n = 4): one matrix representative per conjugacy class of PSp(8,3), in the adapted symplectic basis
# (u1,v1,...,u4,v4), om(u_i,v_i) = 1.  Per class: projective order, class size, reality, M = B g B^-1 (row action).
Run := function(n)
  local G, F, form, V, cands, om, basis, i, found, u, v, ptsV, act, Gp, cls, c, g, B, Binv, m, M, real;
  G := Sp(2*n, 3); F := GF(3);
  form := InvariantBilinearForm(G).matrix;
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
  Print("BASISCHECK n=", n, " ", List(B*form*TransposedMat(B), r -> List(r, IntFFE)), "\n");
  ptsV := Set(List(cands, NormedRowVector));
  act := ActionHomomorphism(G, ptsV, OnLines);
  Gp := Image(act);
  cls := ConjugacyClasses(Gp);
  Print("n=", n, " |PSp|=", Size(Gp), " classes=", Length(cls), "\n");
  for c in cls do
    g := Representative(c);
    real := g^-1 in c;
    m := PreImagesRepresentative(act, g);
    M := B * m * Binv;
    Print("CLASS n=", n, " order=", Order(g), " size=", Size(c), " fixed=-1 real=", real,
          " M=", List(M, r -> List(r, IntFFE)), "\n");
  od;
end;
Run(4);
QUIT;
