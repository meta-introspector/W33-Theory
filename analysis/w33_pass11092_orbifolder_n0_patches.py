#!/usr/bin/env python3
"""Pass 11092: the five patches that let orbifolder 1.2.1 compute N=0 (non-supersymmetric) heterotic orbifold spectra.

orbifolder 1.2.1 (Nilles, Ramos-Sanchez, Vaudrevange, Wingerter, arXiv:1110.5229) loads a twist whose sign combinations
leave no supercharge (N=0), passes modular invariance, and then segfaults or aborts in the spectrum code.  Each fault
below was localised by a debugger or a debug print, and each fix was tested against the supersymmetric control (the
N=1 spectra of all 87 Z6-I W(3,3) models are byte-identical with and without the patches).

  1. CState::FindSUSYMultiplets         -- no N=0 branch.  Fix: every right-mover is its own "multiplet", classified
                                           by helicity q_sh[0] (-1/2 left-handed Weyl fermion, +1/2 right-handed,
                                           -+1 vector, -+2 graviton, 0 real scalar pair -> labelled Hyper "H").
  2. CFixedBrane::FindSUSYMultiplets    -- RecursiveCounting over MaxDigits[0] of an empty vector.  Fix: N=0 has one,
                                           empty, supercharge combination.
  3. CSector::SortByEigenvalue          -- a right-moving oscillator of number operator N contributes -N (mod 1) to the
                                           right-mover eigenvalue, but invariance under the constructing element
                                           (level matching) needs +N.  In SUSY models no massless right-mover ever
                                           carries an oscillator, so the sign never mattered.  In the Z6-I twins the
                                           massless scalar q_sh = (0,1/6,1/6,-1/3) with N_R = 1/6 exposed it: the
                                           mismatch was exactly 1 + 2 N_R = 1/3 (mod 1), in every model.
  4. CState::CreateRepresentations      -- a gauge-neutral untwisted fermion was relabelled as a modulus and its unset
                                           q_sh read (segfault).  Fix: no moduli relabelling at N=0; q_sh = the state's
                                           single right-moving weight.
  5. COrbifold::Create (anomalous U(1)) -- the anomalous U(1) generator t = (1/12) sum p_sh (over left-handed fermions)
                                           was built only for N=1, so the anomaly check read an unrotated basis and
                                           reported 21/87 twins "not universal" (all of them are universal).  Fix:
                                           build it for N=0 too.

Usage:  py -3 analysis/w33_pass11092_orbifolder_n0_patches.py <orbifolder-1.2.1/src/orbifolder> <out_dir>
writes patched copies of cstate.cpp, cfixedpoint.cpp, csector.cpp, corbifold.cpp into out_dir.
"""
from __future__ import annotations

import sys
from pathlib import Path

PATCHES = {
    "cstate.cpp": [
        (r'''  else
  {
    cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): Case N = " << NumberOfSUSY << " SUSY not defined. Return false." << endl;''',
         r'''  else
  if (NumberOfSUSY == 0)
  {
    // PATCH 1 (Pass 11092): no supercharges -- every right-moving state is its own "multiplet",
    // classified by its spacetime helicity q[0].
    for (j = 0; j < s1; ++j)
    {
      const vector<CVector> &SetOfqCharges = this->qCharges[j];
      const double h = SetOfqCharges[0][0];
      if (fabs(h + 0.5) < prec)      Multiplets.push_back(LeftChiral);
      else if (fabs(h - 0.5) < prec) Multiplets.push_back(RightChiral);
      else if (fabs(h + 1.0) < prec) Multiplets.push_back(Vector);
      else if (fabs(h - 1.0) < prec) Multiplets.push_back(VectorCC);
      else if (fabs(h + 2.0) < prec) Multiplets.push_back(Gravity);
      else if (fabs(h - 2.0) < prec) Multiplets.push_back(GravityCC);
      else if (fabs(h) < prec)       Multiplets.push_back(Hyper);   // N=0: a real scalar pair, labelled "H"
      else                           Multiplets.push_back(NOT_DEF_SUSY);
    }
  }
  else
  {
    cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): Case N = " << NumberOfSUSY << " SUSY not defined. Return false." << endl;'''),
        (r'''      if (NewField.SGElement.IsZero() && IsSinglet(NewField.Dimensions) && NewField.U1Charges.IsZero())''',
         r'''      // PATCH 4 (Pass 11092): with N=0 there are no moduli multiplets -- keep fermions and scalars as they are
      if ((Vacuum.InvariantSupercharges.size() != 0) && NewField.SGElement.IsZero() && IsSinglet(NewField.Dimensions) && NewField.U1Charges.IsZero())'''),
        (r'''      else
      {
        for (k = 0; k < s2; ++k)
        {
          if (fabs(SetOfqCharges[k][0]) < prec)
            NewField.q_sh = SetOfqCharges[k];
        }
      }''',
         r'''      else
      if ((Vacuum.InvariantSupercharges.size() == 0) && (s2 != 0))
      {
        // PATCH 4 (Pass 11092): N=0 -- every state is its own multiplet with a single q_sh
        NewField.q_sh = SetOfqCharges[0];
      }
      else
      {
        for (k = 0; k < s2; ++k)
        {
          if (fabs(SetOfqCharges[k][0]) < prec)
            NewField.q_sh = SetOfqCharges[k];
        }
      }'''),
    ],
    "cfixedpoint.cpp": [
        (r'''  else
  {
    vector<unsigned> MaxDigits(NumberOfSUSY, 3);''',
         r'''  else
  if (NumberOfSUSY == 0)
  {
    // PATCH 2 (Pass 11092): no supercharges -- only the empty (zero) combination
    AllNumbers.push_back(vector<unsigned>());
  }
  else
  {
    vector<unsigned> MaxDigits(NumberOfSUSY, 3);'''),
    ],
    "csector.cpp": [
        (r'''        Eigenvalue = ZMxZN_Twists[k] * Weight;

        // begin: add the transformation of the oscillators
        if (with_excitation)
          Eigenvalue += OsciEigenvalues[k];''',
         r'''        Eigenvalue = ZMxZN_Twists[k] * Weight;

        // begin: add the transformation of the oscillators
        // PATCH 3 (Pass 11092): a right-moving oscillator of number operator N contributes -N (mod 1) to
        // GetTransformationProperty, but level matching (invariance under the constructing element) needs
        // +N. Massless right-movers never carry oscillators in SUSY models, so the sign never mattered there.
        if (with_excitation)
          Eigenvalue -= OsciEigenvalues[k];'''),
    ],
    "corbifold.cpp": [
        (r'''  if (CreateAnomalousU1Generator && (this->OrbifoldGroup.GetNumberOfSupersymmetry() == 1))''',
         r'''  // PATCH 5 (Pass 11092): N=0 too -- the sum runs over left-handed Weyl fermions (right-mover helicity -1/2)
  if (CreateAnomalousU1Generator && ((this->OrbifoldGroup.GetNumberOfSupersymmetry() == 1) || (this->OrbifoldGroup.GetNumberOfSupersymmetry() == 0)))'''),
    ],
}


def apply(src: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    done = {}
    for fn, pairs in PATCHES.items():
        t = (src / fn).read_text()
        for old, new in pairs:
            assert t.count(old) == 1, (fn, old[:60])
            t = t.replace(old, new)
        (out / fn).write_text(t)
        done[fn] = len(pairs)
    return done


if __name__ == "__main__":
    print(apply(Path(sys.argv[1]), Path(sys.argv[2])))
