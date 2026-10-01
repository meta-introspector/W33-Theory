# Pass 11225: every 3-dimensional irreducible representation of every subgroup (up to conjugacy) of W(E6) = Aut(W33)
# and of Sp(4,3) (the linear symmetry of two qutrits).  Prints the image group (all matrices, cyclotomic entries) once
# per (subgroup class, 3-dim irrep), with its structure description, for the residual-symmetry mixing scan.
SetUserPreference("AtlasRep", "AtlasRepAccessRemoteFiles", false);
Emit := function(label, G)
  local cc, k, H, chi, rho, gens, img, els, M, j, sd;
  cc := ConjugacyClassesSubgroups(G);
  Print("GROUP ", label, " order=", Size(G), " subgroup_classes=", Length(cc), "\n");
  for k in [1..Length(cc)] do
    H := Representative(cc[k]);
    for j in [1..Length(Irr(H))] do
      chi := Irr(H)[j];
      if chi[1] = 3 then
        rho := IrreducibleRepresentationsDixon(H, chi);
        gens := List(GeneratorsOfGroup(H), h -> h^rho);
        img := Group(gens);
        els := Elements(img);
        sd := StructureDescription(H);
        Print("REP ", label, " sub=", k, " order=", Size(H), " struct=", ReplacedString(sd, " ", ""),
              " irr=", j, " image_order=", Length(els), " image_struct=",
              ReplacedString(StructureDescription(img), " ", ""), "\n");
        for M in els do
          Print("M ", JoinStringsWithSeparator(List(Flat(M), String), ";"), "\n");
        od;
      fi;
    od;
  od;
end;
W := WeylGroup(RootSystem(SimpleLieAlgebra("E", 6, Rationals)));
Emit("WE6", Image(IsomorphismPermGroup(W)));
S := Sp(4, 3);
Emit("Sp43", Image(IsomorphismPermGroup(S)));
QUIT;
