# Pass 11188: one matrix representative for every conjugacy class of PSp(2n,3) (n = 2, 3), written in the adapted
# symplectic basis (u1,v1,...,un,vn) with om(u_i,v_i) = 1 used by analysis/w33_pass11180_mereology.py, so Python can
# evaluate the intrinsic arrow exactly on each class.  Per class: projective order, class size in PSp, number of fixed
# tensor factorisations (must reproduce data/w33_pass11180_gap_classes.txt), whether the class is real (g ~ g^-1), and
# the matrix M = B g B^-1 (GAP row action; Python uses S = M^T on column vectors).
Run := function(n)
  local G, F, form, V, cands, om, planes, basis, i, found, u, v, ptsV, act, Gp, fac, P, e, orb, cls, c, g, fx,
        B, Binv, m, M, real;
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
          if ForAll(basis, b -> IsZero(om(v, b))) and om(u, v) = One(F) then
            Add(planes, [u, v]); Append(basis, [u, v]); found := true; break;
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
    real := g^-1 in c;
    m := PreImagesRepresentative(act, g);
    M := B * m * Binv;
    Print("CLASS n=", n, " order=", Order(g), " size=", Size(c), " fixed=", fx, " real=", real,
          " M=", List(M, r -> List(r, IntFFE)), "\n");
  od;
end;
Run(2);
Run(3);
QUIT;
