Print("Sp43 Sylow3 ", IdGroup(SylowSubgroup(Sp(4,3),3)), "\n");
Print("PSp43 Sylow3 ", IdGroup(SylowSubgroup(PSp(4,3),3)), "\n");
Print("C3wrC3 ", IdGroup(WreathProduct(CyclicGroup(IsPermGroup,3),CyclicGroup(IsPermGroup,3))), "\n");
# Codex U81 law on H27 x F3: (a,b,c,d)(A,B,C,D) = (a+A, b+B, c+C-A*b, d+D+s*(A^2*b-2*A*c))
for s in [1,2] do
  els := Cartesian([0..2],[0..2],[0..2],[0..2]);
  mul := function(x,y) return [ (x[1]+y[1]) mod 3, (x[2]+y[2]) mod 3, (x[3]+y[3]-y[1]*x[2]) mod 3,
        (x[4]+y[4]+s*(y[1]^2*x[2]-2*y[1]*x[3])) mod 3 ]; end;
  perms := List(els, g -> PermList(List(els, x -> Position(els, mul(x,g)))));
  H := Group(perms);
  Print("U81 s=", s, " order ", Size(H), " id ", IdGroup(H), "\n");
od;
QUIT;
