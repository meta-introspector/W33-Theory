# Passes 11180/11181: for every conjugacy class of PSp(2n,3) (n = 2, 3) -- class size, element order, and the number of
# tensor factorisations the class representative FIXES (= the subsystem splits in which the tick is local).
# A class with no fixed factorisation is intrinsically entangling.  Uses GAP's InvariantBilinearForm.
Run := function(n)
  local G, F, form, V, cands, om, planes, basis, i, found, u, v, ptsV, act, Gp, fac, P, e, orb, cls, c, g, fx, tot, free;
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
            Add(planes, [u, v]); Append(basis, [u, v]); found := true; break;
          fi;
        od;
      fi;
    od;
  od;
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
  tot := 0; free := 0;
  for c in cls do
    g := Representative(c);
    fx := Number(orb, f -> OnSets(f, g) = f);
    if fx = 0 then free := free + Size(c); fi;
    Print("CLASS order=", Order(g), " size=", Size(c), " fixed=", fx, "\n");
  od;
  Print("n=", n, " fixed-point-free (intrinsically entangling) elements=", free, " fraction=", free / Size(Gp), "\n");
end;
Run(2);
Run(3);
QUIT;
