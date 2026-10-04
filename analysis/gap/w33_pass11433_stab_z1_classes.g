Read("w33_pass11433_gens4.g");  # generators of Sp(8,3) in row convention (H, S, SUM)
G := Group(gens4);
z1 := List([1..8], i -> 0*Z(3)); z1[2] := Z(3)^0;
H := Stabilizer(G, z1, OnRight);
Print("H size ", Size(H), "\n");
cc := ConjugacyClasses(H);
Print("classes ", Length(cc), "\n");
out := OutputTextFile("w33_pass11433_stab_z1_classes_n4.txt", false);
SetPrintFormattingStatus(out, false);
for c in cc do
  AppendTo(out, List(Representative(c), r -> List(r, IntFFE)), ";", Size(c), "\n");
od;
CloseStream(out);
Print("done\n");
QUIT;
