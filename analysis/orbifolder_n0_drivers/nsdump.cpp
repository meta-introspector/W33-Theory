#include <stdio.h>
#include <cstdlib>
#include "cprompt.h"
#include "cspectrum.h"
#include "canalysemodel.h"
using namespace std;
unsigned SELFDUALLATTICE;
int main(int argc, char *argv[])
{
  ifstream in(argv[1]);
  string prog = "";
  COrbifoldGroup G;
  if (!G.LoadOrbifoldGroup(in, prog)) { cout << "LOADFAIL" << endl; return 2; }
  cout << "NSUSY " << G.GetNumberOfSupersymmetry() << endl;
  COrbifold O(G);
  const SConfig &c = O.StandardConfig;
  cout << "GAUGE " << c.SymmetryGroup.GaugeGroup.algebra << " NU1 " << c.SymmetryGroup.GaugeGroup.u1directions.size() << endl;
  map<string,int> cnt;
  for (unsigned i = 0; i < c.Fields.size(); ++i)
  {
    const CField &F = c.Fields[i];
    ostringstream os;
    os << "m=" << (int)F.Multiplet << " dim=";
    for (unsigned d = 0; d < F.Dimensions.size(); ++d) os << (d ? "," : "") << F.Dimensions[d].Dimension << F.Dimensions[d].AdditionalLabel;
    cnt[os.str()]++;
  }
  for (auto &kv : cnt) cout << "F " << kv.first << " x" << kv.second << endl;
  cout << "NFIELDS " << c.Fields.size() << endl;
  return 0;
}
