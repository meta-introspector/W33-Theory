# Pass 11178: the 20 orbitals of Sp(6,3) on three-qutrit tensor factorisations, each with its Choi-entanglement
# signature.  For the base factorisation F0 = (P1,P2,P3) and a representative F' = (Q1,Q2,Q3) of each suborbit, the
# block (i,j) is the projection of Q_j onto P_i along the other two P's, written in symplectic bases (omega = 1) of
# both planes: code 0 = zero, 1 = rank one, 2 = invertible with det +1, 3 = invertible with det -1.
# Uses GAP's own invariant form (InvariantBilinearForm).
n := 3;
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
          v := v / om(u, v);
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
sorb := Orbits(stab, [1..Length(orb)]);
B0 := Concatenation(planes);
B0inv := B0^-1;
PlanesOf := function(s)
  local vs, pls, a, b, sp, key;
  vs := List(s, i -> ptsV[i]);
  pls := [];
  for a in vs do for b in vs do
    if a <> b then
      sp := Set(List(Filtered(Elements(VectorSpace(F, [a, b])), x -> not IsZero(x)), NormedRowVector));
      if Length(sp) = 4 and IsSubset(Set(vs), sp) and not sp in pls then Add(pls, sp); fi;
    fi;
  od; od;
  return pls;
end;
SymBasis := function(pl)
  local a, b;
  a := pl[1];
  for b in pl do
    if not IsZero(om(a, b)) then return [a, b / om(a, b)]; fi;
  od;
end;
Print("factorisations=", Length(orb), " rank=", Length(sorb), "\n");
for o in sorb do
  r := o[1];
  Q := List(PlanesOf(orb[r]), SymBasis);
  code := [];
  for i in [1..n] do
    row := [];
    for j in [1..n] do
      M := List(Q[j], w -> (w*B0inv){[2*i-1, 2*i]});
      if IsZero(M) then Add(row, 0);
      elif RankMat(M) = 1 then Add(row, 1);
      elif DeterminantMat(M) = One(F) then Add(row, 2);
      else Add(row, 3); fi;
    od;
    Add(code, row);
  od;
  Print("SUB ", Length(o), " ", code, "\n");
od;
QUIT;
