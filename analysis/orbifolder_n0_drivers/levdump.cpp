// dump every field (sector k,l, multiplet, dims, U(1) charges) of the StandardConfig -- used with ORB_MASS_LEVEL
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
  int n = 0;
  while (G.LoadOrbifoldGroup(in, prog))
  {
    ++n;
    COrbifold O(G);
    const SConfig &c = O.StandardConfig;
    cout << "MODEL " << n << " " << G.Label << " NSUSY " << G.GetNumberOfSupersymmetry() << " NU1 " << c.SymmetryGroup.GaugeGroup.u1directions.size() << endl;
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
      cout << "U1STD " << u << " |";
      for (unsigned x = 0; x < 16; ++x) { char b[64]; snprintf(b, 64, " %.17g", (double)GG.u1directions[u][x]); cout << b; }
      cout << endl;
    }
    for (unsigned i = 0; i < c.Fields.size(); ++i)
    {
      const CField &F = c.Fields[i];
      cout << "S k=" << F.SGElement.Get_k() << " l=" << F.SGElement.Get_l() << " m=" << (int)F.Multiplet << " dim=";
      for (unsigned d = 0; d < F.Dimensions.size(); ++d) cout << (d ? "," : "") << F.Dimensions[d].Dimension << F.Dimensions[d].AdditionalLabel;
      cout << " q=";
      for (unsigned a = 0; a < F.U1Charges.size(); ++a) { char b[64]; snprintf(b, 64, "%s%.17g", a ? "," : "", (double)F.U1Charges[a]); cout << b; }
      cout << " osc=" << F.OsciContribution[1] << "," << F.OsciContribution[2] << "," << F.OsciContribution[3] << endl;
    }
    cout << "END " << n << " " << c.Fields.size() << endl;
  }
  return 0;
}
