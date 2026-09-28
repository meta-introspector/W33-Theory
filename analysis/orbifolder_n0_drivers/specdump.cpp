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
  ostringstream devnull;
  CPrint Quiet(Tstandard, &devnull);
  string prog = "";
  COrbifoldGroup G;
  if (!G.LoadOrbifoldGroup(in, prog)) return 2;
  COrbifold O(G);
  vector<SConfig> cfg;
  bool bSM = true, bPS = false, bSU5 = false;
  CAnalyseModel A;
  A.AnalyseModel(O, O.StandardConfig, bSM, bPS, bSU5, cfg, Quiet, 3, false);
  if (!bSM || cfg.empty()) { cout << "noSM" << endl; return 3; }
  SConfig c = cfg[0];
  unsigned L = c.use_Labels;
  cout.precision(17);
  cout << "NU1 " << c.SymmetryGroup.GaugeGroup.u1directions.size() << endl;
  cout << "ANOM " << (c.SymmetryGroup.IsFirstU1Anomalous ? 1 : 0)
       << " FI " << c.SymmetryGroup.D0_FI_term << endl;
  for (unsigned i = 0; i < c.Fields.size(); ++i)
  {
    const CField &F = c.Fields[i];
    cout << "S " << F.Labels[L] << "_" << F.Numbers[L] << " k=" << F.SGElement.Get_k() << " susy=" << (int)F.Multiplet << " dim=";
    for (unsigned d = 0; d < F.Dimensions.size(); ++d)
      cout << (d ? "," : "") << F.Dimensions[d].Dimension << F.Dimensions[d].AdditionalLabel;
    cout << " q=";
    for (unsigned a = 0; a < F.U1Charges.size(); ++a) cout << (a ? "," : "") << F.U1Charges[a];
    cout << endl;
  }
  return 0;
}
