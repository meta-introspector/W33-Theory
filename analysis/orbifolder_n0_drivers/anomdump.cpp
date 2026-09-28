#include <stdio.h>
#include <cstdlib>
#include "cprompt.h"
#include "cspectrum.h"
#include "canalysemodel.h"
#include "cgaugeinvariance.h"
using namespace std;
unsigned SELFDUALLATTICE;
int main(int argc, char *argv[])
{
  ifstream in(argv[1]);
  string prog = "";
  COrbifoldGroup G;
  if (!G.LoadOrbifoldGroup(in, prog)) { cout << "LOADFAIL" << endl; return 2; }
  COrbifold O(G);
  CPrint P(Tstandard, &cout);
  CGaugeIndices GI;
  SConfig c = O.StandardConfig;
  bool ok = O.CheckAnomaly(c, GI, P, true);
  cout << "\nANOMALY_FREE_OR_GS " << ok << " NSUSY " << G.GetNumberOfSupersymmetry()
       << " anomU1 " << c.SymmetryGroup.IsFirstU1Anomalous << endl;
  return 0;
}
