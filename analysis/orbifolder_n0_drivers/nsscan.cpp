// Pass 11093: tachyon-free filter over the Z6-I W(3,3) scan.
// usage: nsscan base.txt tries seed out_prefix            (random Wilson lines, as the original scan)
//        nsscan --eval models.txt out_prefix              (evaluate every model in a file, no randomisation)
// Tachyon test (Pass 11092 lemma): the non-SUSY twin is tachyonic iff the SUSY parent has a massless
// oscillator-excited left mover in the theta sector (k=1, l=0).  Tachyon-free twins are built (twist
// (0,1/6,1/6,2/3), same shifts and Wilson lines) and put through the Standard-Model test.
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "cprompt.h"
#include "cspectrum.h"
#include "cinequivalentspectra.h"
#include "crandommodel.h"
#include "corbifoldcore.h"
#include "canalysemodel.h"
using namespace std;
unsigned SELFDUALLATTICE;

static int theta_excited(const COrbifold &O, int &n_fp, int &n_states)
{
  int excited = 0; n_fp = 0; n_states = 0;
  const vector<CSector> &S = O.GetSectors();
  for (unsigned i = 0; i < S.size(); ++i)
  {
    if (S[i].Get_k() != 1 || S[i].Get_l() != 0) continue;
    for (unsigned j = 0; j < S[i].GetNumberOfFixedBranes(); ++j)
    {
      const CFixedBrane &FB = S[i].GetFixedBrane(j);
      ++n_fp;
      for (unsigned k = 0; k < FB.GetNumberOfInvariantStates(); ++k)
      {
        const CState &st = FB.GetInvariantState(k);
        const CMasslessHalfState &LM = FB.GetMasslessLeftMover(st.GetLeftMover().GetIndex());
        ++n_states;
        if (LM.Excitation.NumberOperator > 1e-6) ++excited;
      }
    }
  }
  return excited;
}

static string twin_file(const COrbifoldGroup &N, const string &tmp)
{
  { ofstream o(tmp.c_str()); N.PrintToFile(o); }
  ifstream i(tmp.c_str()); stringstream ss; ss << i.rdbuf(); string s = ss.str();
  size_t p = s.find("Geometry_Z6-I_G2xG2xSU3.txt");
  if (p == string::npos) return "";
  s.replace(p, string("Geometry_Z6-I_G2xG2xSU3.txt").size(), "Geometry_Z6-INS_G2xG2xSU3.txt");
  string tf = tmp + ".ns";
  ofstream o2(tf.c_str()); o2 << s; o2.close();
  return tf;
}

int main(int argc, char *argv[])
{
  if (argc < 3) { cerr << "usage" << endl; return 1; }
  ostringstream devnull;
  CPrint Print(Tstandard, &devnull);
  string prog = "";
  bool eval = (string(argv[1]) == "--eval");
  string prefix = eval ? argv[3] : argv[4];
  string tmp = prefix + ".tmp";
  ofstream tf_out((prefix + "_tachyonfree.txt").c_str()), sm_out((prefix + "_twinSM.txt").c_str()), log((prefix + ".log").c_str());
  CAnalyseModel Analyse;
  CInequivalentModels TwinSMs, TF, ParentSMs;
  ofstream psm_out((prefix + "_parentSM.txt").c_str());

  vector<COrbifoldGroup> todo;
  COrbifoldGroup G;
  unsigned tries = 0, seed = 0;
  if (eval)
  {
    ifstream in(argv[2]);
    while (true) { COrbifoldGroup g; if (!g.LoadOrbifoldGroup(in, prog)) break; todo.push_back(g); }
    tries = todo.size();
    cout << "loaded " << tries << " models" << endl;
    if (tries == 0) return 3;
  }
  else
  {
    ifstream in(argv[1]);
    if (!G.LoadOrbifoldGroup(in, prog)) { cerr << "load failed" << endl; return 2; }
    G.LoadedU1Generators.clear();
    tries = atoi(argv[2]); seed = atoi(argv[3]);
  }
  vector<bool> use(8, false); use[0] = use[1] = true;
  CRandomModel R(eval ? todo[0].GetLattice() : G.GetLattice());
  vector<CVector> Unbroken;
  if (!eval) { R.Initiate(G, use, Unbroken); srand(seed); }

  vector<CSector> Core, CoreT;
  COrbifoldCore OC(eval ? todo[0] : G, Core);
  bool haveCoreT = false;
  unsigned ok = 0, tach = 0, tfree = 0, psm = 0, tsm = 0, tsm_ineq = 0, tf_ineq = 0;
  for (unsigned t = 1; t <= tries; ++t)
  {
    COrbifoldGroup N = eval ? todo[t - 1] : G;
    if (!eval)
    {
      if (!N.CreateRandom(R, false)) continue;
      if (N.GetModularInvariance_CheckStatus() != CheckedAndGood) continue;
    }
    COrbifold O(N, Core);
    if (O.GetCheckStatus() != CheckedAndGood) continue;
    ++ok;
    int nfp = 0, nst = 0;
    int ex = theta_excited(O, nfp, nst);
    bool parentSM = true, bPS = false, bSU5 = false;
    vector<SConfig> cfg;
    Analyse.AnalyseModel(O, O.StandardConfig, parentSM, bPS, bSU5, cfg, Print);
    if (parentSM)
    {
      ++psm;
      // which SM-labelled fields of the parent sit in oscillator-excited theta-sector states?
      const SConfig &c = cfg[0];
      unsigned L = c.use_Labels;
      ostringstream ex_lab, all_lab;
      int nex = 0, nth = 0;
      for (unsigned f = 0; f < c.Fields.size(); ++f)
      {
        const CField &F = c.Fields[f];
        if ((F.SGElement.Get_k() != 1 && F.SGElement.Get_k() != 5) || F.SGElement.Get_l() != 0) continue;
        if (F.Multiplet != LeftChiral) continue;
        ++nth;
        bool exc = !F.OsciContribution.IsZero();
        all_lab << "k" << F.SGElement.Get_k() << ":" << F.Labels[L] << "_" << F.Numbers[L] << (exc ? "*" : "") << ",";
        if (exc) { ++nex; ex_lab << "k" << F.SGElement.Get_k() << ":" << F.Labels[L] << "_" << F.Numbers[L] << ","; }
      }
      log << "P " << t << " " << N.Label << " theta_fields " << nth << " excited " << nex << " : " << ex_lab.str() << " | all " << all_lab.str() << endl;
      CSpectrum SP(O.StandardConfig, LeftChiral);
      if (ParentSMs.IsSpectrumUnknown(SP, true)) { ostringstream os; os << "PSM_" << seed << "_" << t; COrbifoldGroup NN = N; if (!eval) NN.Label = os.str(); NN.PrintToFile(psm_out); psm_out.flush(); }
    }
    string lab = N.Label;
    if (ex > 0) { ++tach; log << "T " << t << " " << lab << " excited " << ex << " of " << nst << " fp " << nfp << " parentSM " << parentSM << endl; continue; }
    ++tfree;
    CSpectrum S(O.StandardConfig, LeftChiral);
    bool newtf = TF.IsSpectrumUnknown(S, true);
    if (newtf) { ++tf_ineq; ostringstream os; os << "TF_" << seed << "_" << t; if (!eval) N.Label = os.str(); N.PrintToFile(tf_out); tf_out.flush(); }
    string tfile = twin_file(N, tmp);
    ifstream tin(tfile.c_str());
    COrbifoldGroup NT;
    if (!NT.LoadOrbifoldGroup(tin, prog)) { log << "F " << t << " twin load failed" << endl; continue; }
    if (!haveCoreT) { COrbifoldCore OCT(NT, CoreT); haveCoreT = true; }
    COrbifold OT(NT, CoreT);
    if (OT.GetCheckStatus() != CheckedAndGood) { log << "F " << t << " twin orbifold failed" << endl; continue; }
    bool twinSM = true; bPS = false; bSU5 = false;
    vector<SConfig> cfg2;
    Analyse.AnalyseModel(OT, OT.StandardConfig, twinSM, bPS, bSU5, cfg2, Print);
    log << "N " << t << " " << N.Label << " fp " << nfp << " states " << nst << " parentSM " << parentSM << " twinSM " << twinSM
        << " nsusy " << NT.GetNumberOfSupersymmetry() << " newTF " << newtf << endl;
    if (twinSM)
    {
      ++tsm;
      CSpectrum ST(OT.StandardConfig, LeftChiral);
      if (TwinSMs.IsSpectrumUnknown(ST, true)) { ++tsm_ineq; NT.PrintToFile(sm_out); sm_out.flush(); }
    }
    if (t % 200 == 0 || t == tries)
      cout << "tries " << t << " ok " << ok << " tachyonic " << tach << " tachyon-free " << tfree << " (ineq " << tf_ineq
           << ") parentSM " << psm << " twinSM " << tsm << " (ineq " << tsm_ineq << ")" << endl;
  }
  cout << "FINAL tries " << tries << " ok " << ok << " tachyonic " << tach << " tachyon-free " << tfree << " (ineq " << tf_ineq
       << ") parentSM " << psm << " twinSM " << tsm << " (ineq " << tsm_ineq << ")" << endl;
  return 0;
}
