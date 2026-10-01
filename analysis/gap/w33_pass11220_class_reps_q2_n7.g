# Pass 11220 (n = 7): one matrix per conjugacy class of Sp(2n,q) (q = 2, 5), in an adapted symplectic basis (u1,v1,...),
# om(u_i,v_i) = 1, for the arrow law A = n - c beyond qutrits.  (For q = 2, Sp = PSp; for odd q we list classes of
# PSp acting on points and print a matrix preimage.)
Run := function(n, q)
  local G, F, form, V, cands, om, basis, i, found, u, v, B, Binv, ptsV, act, Gp, cls, c, g, m, M;
  G := Sp(2*n, q); F := GF(q);
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
  ptsV := Set(List(cands, NormedRowVector));
  act := ActionHomomorphism(G, ptsV, OnLines);
  Gp := Image(act);
  cls := ConjugacyClasses(Gp);
  Print("q=", q, " n=", n, " |PSp|=", Size(Gp), " classes=", Length(cls), "\n");
  for c in cls do
    g := Representative(c);
    m := PreImagesRepresentative(act, g);
    M := B * m * Binv;
    Print("CLASS q=", q, " n=", n, " order=", Order(g), " size=", Size(c), " M=", List(M, r -> List(r, IntFFE)), "\n");
  od;
end;
Run(7, 2);
QUIT;
