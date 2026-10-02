F := CyclotomicField(9); z := E(9);
Xs := [[0,0,1],[1,0,0],[0,1,0]]*One(F);
T := DiagonalMat([1, z, z^8]);
Sg := DiagonalMat([1, 1, z^3]);
G := Group(Xs, Sg, T);
Sc := Filtered(Elements(Centre(G)), g -> IsDiagonalMat(g) and g[1][1]=g[2][2] and g[2][2]=g[3][3]);
P := G / Subgroup(G, Sc);
Print("order ", Size(P), " id ", IdGroup(P), "\n");
for k in [7..10] do
  H := SmallGroup(81,k);
  Print(k, " order9 ", Number(Elements(H), g -> Order(g)=9), " exp ", Exponent(H), " class ", NilpotencyClassOfGroup(H), "\n");
od;
# U81 law: build from Codex mul via pc presentation is harder; count invariants instead
QUIT;
