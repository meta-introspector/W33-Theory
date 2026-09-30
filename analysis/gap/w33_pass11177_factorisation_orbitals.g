# Pass 11177: orbitals of Sp(2n,3) on two- and three-qutrit tensor factorisations.
# A factorisation = the set of projective points of its n mutually orthogonal nondegenerate planes (8 points for n = 2,
# 12 for n = 3).  Sp(2n,3) is first turned into a permutation group on the projective points (OnLines), then acts on
# factorisations by OnSets.  Uses GAP's own invariant form (InvariantBilinearForm), never a hand-written one.
Run := function(n)
  local G, F, form, V, planes, basis, cands, u, v, i, found, P, e, ptsV, pts, act, Gp, fac, orb, hom, img, stab, sub;
  G := Sp(2*n, 3);
  F := GF(3);
  form := InvariantBilinearForm(G).matrix;
  V := F^(2*n);
  planes := []; basis := [];
  cands := Filtered(Elements(V), x -> not IsZero(x));
  for i in [1..n] do
    found := false;
    for u in cands do
      if found then break; fi;
      if ForAll(basis, b -> IsZero(u*form*b)) then
        for v in cands do
          if ForAll(basis, b -> IsZero(v*form*b)) and not IsZero(u*form*v) then
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
  hom := ActionHomomorphism(Gp, orb, OnSets);
  img := Image(hom);
  stab := Stabilizer(img, 1);
  sub := List(Orbits(stab, [1..Length(orb)]), Length);
  Print("n=", n, " |Sp|=", Size(G), " points=", Length(ptsV), " factorisations=", Length(orb),
        " rank=", Length(sub), " |PSp image|=", Size(img), "\n");
  Print("subdegrees=", SortedList(sub), "\n");
  Print("point stabiliser order=", Size(stab), "\n");
  if n = 2 then Print("point stabiliser structure=", StructureDescription(stab), "\n"); fi;
end;
Run(2);
Run(3);
QUIT;
