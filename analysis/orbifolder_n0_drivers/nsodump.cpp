#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include "corbifold.h"
#include "corbifoldgroup.h"
#include "cprint.h"
#include "cgaugeinvariance.h"
using namespace std;
unsigned SELFDUALLATTICE;
static void dumpcfg(const char *tag, const SConfig &c)
{
  map<string,int> cnt;
  for (unsigned i = 0; i < c.Fields.size(); ++i)
  {
    const CField &F = c.Fields[i];
    ostringstream os;
    os << "m=" << (int)F.Multiplet << " k=" << F.SGElement.Get_k() << " m=" << F.SGElement.Get_m() << " n=" << F.SGElement.Get_n() << " dim=";
    for (unsigned d = 0; d < F.Dimensions.size(); ++d) os << (d ? "," : "") << F.Dimensions[d].Dimension << F.Dimensions[d].AdditionalLabel;
    cnt[os.str()]++;
  }
  for (map<string,int>::iterator it = cnt.begin(); it != cnt.end(); ++it) cout << tag << " " << it->first << " x" << it->second << endl;
  cout << tag << "_NFIELDS " << c.Fields.size() << endl;
}
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
  dumpcfg("F", c);
  dumpcfg("T", O.TachyonicStandardConfig);
  CPrint P(Tstandard, &cout);
  CGaugeIndices GI;
  SConfig cc = c;
  bool ok = O.CheckAnomaly(cc, GI, P, false);
  cout << "ANOMALY_OK " << ok << endl;
  return 0;
}
