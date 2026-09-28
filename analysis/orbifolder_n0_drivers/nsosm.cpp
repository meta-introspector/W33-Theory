// non-SUSY orbifolder: dump the SM configuration of every model in a file (factors, U(1) directions in 16D, labelled
// left fermions and scalars), full precision.
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "corbifold.h"
#include "corbifoldgroup.h"
#include "canalysemodel.h"
#include "cprint.h"
using namespace std;
unsigned SELFDUALLATTICE;
int main(int argc, char *argv[])
{
  ifstream in(argv[1]);
  string prog = "";
  ostringstream devnull;
  CPrint Quiet(Tstandard, &devnull);
  COrbifoldGroup G;
  int n = 0;
  while (G.LoadOrbifoldGroup(in, prog))
  {
    ++n;
    COrbifold O(G);
    vector<SConfig> cfg;
    bool bSM = true, bPS = false, bSU5 = false;
    CAnalyseModel A;
    A.AnalyseModel(O, O.StandardConfig, bSM, bPS, bSU5, cfg, Quiet);
    cout << "MODEL " << n << " " << G.Label << " SM " << (bSM && !cfg.empty()) << " TACH " << O.TachyonicStandardConfig.Fields.size() << endl;
    if (!bSM || cfg.empty()) continue;
    const SConfig &c = cfg[0];
    const CGaugeGroup &GG = c.SymmetryGroup.GaugeGroup;
    for (unsigned f = 0; f < GG.factor.size(); ++f)
    {
      cout << "FACTOR " << f << " " << GG.factor[f].algebra;
      for (unsigned r = 0; r < GG.factor[f].simpleroots.size(); ++r)
      {
        cout << " |";
        for (unsigned x = 0; x < 16; ++x) { char b[64]; snprintf(b, 64, " %.17g", (double)GG.factor[f].simpleroots[r][x]); cout << b; }
      }
      cout << endl;
    }
    for (unsigned u = 0; u < GG.u1directions.size(); ++u)
    {
      cout << "U1 " << u << " |";
      for (unsigned x = 0; x < 16; ++x) { char b[64]; snprintf(b, 64, " %.17g", (double)GG.u1directions[u][x]); cout << b; }
      cout << endl;
    }
    unsigned L = c.use_Labels;
    for (unsigned i = 0; i < c.Fields.size(); ++i)
    {
      const CField &F = c.Fields[i];
      if (F.Multiplet != LeftFermi && F.Multiplet != Scalar) continue;
      cout << (F.Multiplet == LeftFermi ? "L " : "SC ") << F.Labels[L] << " dim=";
      for (unsigned d = 0; d < F.Dimensions.size(); ++d) cout << (d ? "," : "") << F.Dimensions[d].Dimension << F.Dimensions[d].AdditionalLabel;
      cout << " q=";
      for (unsigned a = 0; a < F.U1Charges.size(); ++a) { char b[64]; snprintf(b, 64, "%s%.17g", a ? "," : "", (double)F.U1Charges[a]); cout << b; }
      cout << endl;
    }
  }
  return 0;
}
