// Pass 11095: Wilson-line scan with the non-SUSY orbifolder (SO(16)xSO(16) string, Z2W x Z3) keeping V0 and the W(3,3)
// A8 Kac pair V1 fixed.  usage: nsoscan base.txt tries seed prefix
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "corbifold.h"
#include "corbifoldgroup.h"
#include "corbifoldcore.h"
#include "crandommodel.h"
#include "canalysemodel.h"
#include "cinequivalentspectra.h"
#include "cspectrum.h"
#include "cprint.h"
using namespace std;
unsigned SELFDUALLATTICE;
int main(int argc, char *argv[])
{
  if (argc < 5) return 1;
  ifstream in(argv[1]);
  unsigned tries = atoi(argv[2]), seed = atoi(argv[3]);
  string prefix = argv[4];
  ostringstream devnull;
  CPrint Print(Tstandard, &devnull);
  string prog = "";
  COrbifoldGroup G;
  if (!G.LoadOrbifoldGroup(in, prog)) { cout << "LOADFAIL" << endl; return 2; }
  G.LoadedU1Generators.clear();
  vector<bool> use(9, false); use[0] = use[1] = use[2] = true;
  CRandomModel R(G.GetLattice());
  vector<CVector> Unbroken;
  R.Initiate(G, use, Unbroken);
  srand(seed);
  vector<CSector> Core;
  COrbifoldCore OC(G, Core);
  CAnalyseModel A;
  CInequivalentModels SMs, TF;
  ofstream smout((prefix + "_SM.txt").c_str()), tfsm((prefix + "_SM_tachyonfree.txt").c_str());
  unsigned ok = 0, tf = 0, sm = 0, smtf = 0, smineq = 0, prob = 0;
  for (unsigned t = 1; t <= tries; ++t)
  {
    COrbifoldGroup N = G;
    if (!N.CreateRandom(R, false)) { ++prob; continue; }
    if (N.GetModularInvariance_CheckStatus() != CheckedAndGood) { ++prob; continue; }
    COrbifold O(N, Core);
    if (O.GetCheckStatus() != CheckedAndGood) { ++prob; continue; }
    ++ok;
    bool tachfree = (O.TachyonicStandardConfig.Fields.size() == 0);
    if (tachfree) ++tf;
    bool bSM = true, bPS = false, bSU5 = false;
    vector<SConfig> cfg;
    A.AnalyseModel(O, O.StandardConfig, bSM, bPS, bSU5, cfg, Print);
    if (bSM)
    {
      ++sm;
      if (tachfree) ++smtf;
      vector<SUSYMultiplet> MT(1, LeftFermi); MT.push_back(Scalar); CSpectrum S(O.StandardConfig, MT);
      if (SMs.IsSpectrumUnknown(S, true))
      {
        ++smineq;
        ostringstream os; os << "A8SM_" << seed << "_" << t << (tachfree ? "_TF" : "_T");
        N.Label = os.str();
        N.PrintToFile(smout); smout.flush();
        if (tachfree) { N.PrintToFile(tfsm); tfsm.flush(); }
      }
    }
    if (t % 500 == 0 || t == tries)
      cout << "tries " << t << " ok " << ok << " tachyon-free " << tf << " SM " << sm << " (ineq " << smineq << ") SM&tachyon-free " << smtf << " problems " << prob << endl;
  }
  return 0;
}
