# Pass 11188: one matrix representative for every conjugacy class of PSp(2n,3) (n = 2, 3), written in an adapted
# symplectic basis (u1,v1,...,un,vn) with om(u_i,v_i) = 1, so that Python can use it in the per-qutrit (x,z) ordering
# of analysis/w33_pass11180_mereology.py.  Printed per class: order, class size, number of fixed tensor
# factorisations (the Pass 11180 cross-check) and the matrix M = B g B^-1 (GAP row-vector action; Python transposes).
Run := function(n)
  local G, F, form, V, cands, om, planes, basis, i, found, u, v, ptsV, act, Gp, fac, P, e, orb, cls, c, g, fx, B, Binv,
        M, m;
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
            Add(planes, [u, v]); Append(basis, [u, v]); found := true; break;
          fi;
        od;
      fi;
    od;
  od;
  B := basis;
  Binv := B^-1;
  Print("BASISCHECK ", List(B, x -> List(B, y -> IntFFE(om(x, y)))), "\n");
  ptsV := Set(List(cands, NormedRowVector));
  act := ActionHomomorphism(G, ptsV, OnLines);
  Gp := Image(act);
  fac := [];
  for P in planes do
    for e in Elements(VectorSpace(F, P)) do
      if not IsZero(e) then AddSet(fac, Position(ptsV, NormedRowVector(e))); fi;
    od;
  od;
  orb := Orbit(Gp, fac, OnSets);
  cls := ConjugacyClasses(Gp);
  Print("n=", n, " |PSp|=", Size(Gp), " factorisations=", Length(orb), " classes=", Length(cls), "\n");
  for c in cls do
    g := Representative(c);
    fx := Number(orb, f -> OnSets(f, g) = f);
    m := PreImagesRepresentative(act, g);
    M := B * m * Binv;
    Print("CLASS n=", n, " order=", Order(g), " size=", Size(c), " fixed=", fx, " M=",
          List(M, r -> List(r, IntFFE)), "\n");
  od;
end;
Run(2);
Run(3);
QUIT;
