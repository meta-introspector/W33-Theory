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
                                           PRIOR ART: the non-SUSY orbifolder (arXiv:2504.20137) carries the same fix
                                           ("wrong trafo sign for R-moving oscillator excitations in original
                                           Orbifolder"); this patch is a rediscovery and the fix is theirs.
  4. CState::CreateRepresentations      -- a gauge-neutral untwisted fermion was relabelled as a modulus and its unset
                                           q_sh read (segfault).  Fix: no moduli relabelling at N=0; q_sh = the state's
                                           single right-moving weight.
  6. CTwistVector::UpdateData (Pass 11095) -- the order of a twist was its geometric order (n v integer); on fermions it is
                                           2n when sum n v_i is odd.  Without it the local modular-invariance check drops
                                           Witten-twisted sectors (found by cross-checking against the non-SUSY
                                           orbifolder, arXiv:2504.20137).  Irrelevant for the Z6-I twins (6 v' has even sum).
  7. CTwistVector::UpdateData (Pass 11094) -- optional mass level mu from ORB_MASS_LEVEL: every non-identity sector solved
                                           at M^2/8 = mu, giving tachyonic spectra through the full projection; 7b/7c make a
                                           sector with no right-movers at that level empty instead of fatal.  Inert when
                                           the variable is unset.
  5. COrbifold::Create (anomalous U(1)) -- the anomalous U(1) generator t = (1/12) sum p_sh (over left-handed fermions)
                                           was built only for N=1, so the anomaly check read an unrotated basis and
                                           reported 21/87 twins "not universal" (all of them are universal).  Fix:
                                           build it for N=0 too.

Usage:  py -3 analysis/w33_pass11092_orbifolder_n0_patches.py <orbifolder-1.2.1/src/orbifolder> <out_dir>
writes patched copies of cstate.cpp, cfixedpoint.cpp, csector.cpp, ctwist.cpp, corbifold.cpp into out_dir.
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
        (r'''#include <iostream>''', r'''#include <iostream>
#include <cstdlib>'''),
        (r'''  const size_t s1 = this->RM_Excitations.size();
  if (s1 == 0)
  {
    cout << "\n  Warning in bool CSector::CreateMasslessRightMover(...) : Set of excitations is empty. Return false.";
    return false;
  }''',
         r'''  const size_t s1 = this->RM_Excitations.size();
  if (s1 == 0)
  {
    // PATCH 7b (Pass 11094): at a mass level mu < 0 a sector may have no right-movers at all -- that is an empty sector,
    // not a failure
    if (getenv("ORB_MASS_LEVEL") != NULL)
      return true;
    cout << "\n  Warning in bool CSector::CreateMasslessRightMover(...) : Set of excitations is empty. Return false.";
    return false;
  }'''),
        (r'''  if (((s1 == 0) || (s2 == 0)) && (this->Twist.OrderOfTwist() != 1))''',
         r'''  // PATCH 6b (Pass 11095): with PATCH 6 a pure (-1)^F sector (integer twist of odd sum, e.g. the Witten twist
  // (0,1,1,1)) has order 2, but it rotates no plane, so it has no fractional oscillators.  Test geometry, not order.
  if (((s1 == 0) || (s2 == 0)) && !(is_integer(this->Twist[1]) && is_integer(this->Twist[2]) && is_integer(this->Twist[3])))'''),
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
    "ctwist.cpp": [
        (r'''#include <cstdlib>''', r'''#include <cstdlib>
#include <cstdio>
#include <cmath>'''),
        (r'''  this->a_L = -1.0 + (0.5 * tmp);
  this->a_R = -0.5 + (0.5 * tmp);''',
         r'''  this->a_L = -1.0 + (0.5 * tmp);
  this->a_R = -0.5 + (0.5 * tmp);

  // PATCH 7 (Pass 11094): optional mass level.  With the environment variable ORB_MASS_LEVEL = "num/den" (mu < 0) every
  // non-identity sector is solved at M^2/8 = mu instead of 0 (both zero-point energies shifted by -mu), so the whole
  // projection machinery (centralisers, gamma phases, level matching) yields the TACHYONIC spectrum at that level.
  // Without the variable nothing changes.
  const char *mass_level = getenv("ORB_MASS_LEVEL");
  if (mass_level != NULL)
  {
    double num = 0.0, den = 1.0;
    if (sscanf(mass_level, "%lf/%lf", &num, &den) < 1) den = 1.0;
    bool identity = true;
    for (i = 0; i < 4; ++i)
      if (fabs((*this)[i]) > prec) identity = false;
    if (!identity)
    {
      this->a_L -= num / den;
      this->a_R -= num / den;
    }
  }'''),
        (r'''    if (Order_found)
    {
      this->Order = n;
      return true;
    }''',
         r'''    if (Order_found)
    {
      // PATCH 6 (Pass 11095): the order on SPINORS.  theta^n is a rotation by 2 pi (n v); if sum_i n v_i is odd it is
      // (-1)^F, so the twist has order 2n on fermions.  SUSY twists (sum v = 0) are unaffected.
      int s = 0;
      for (i = 1; i < 4; ++i)
        s += (int)round_double_to_int((*this)[i] * n);
      this->Order = ((s % 2) != 0) ? 2 * n : n;
      return true;
    }'''),
    ],
    "corbifold.cpp": [
        (r'''#include "chugeint.h"''', r'''#include "chugeint.h"
#include <cstdlib>'''),
        (r'''    if (Sector.GetRightMovers().size() == 0)
    {
      cout << "\n  Warning in bool COrbifold::Create() : Right-movers of the " << i << "-th sector have not been created. Return false." << endl;
      return false;
    }''',
         r'''    if (Sector.GetRightMovers().size() == 0)
    {
      // PATCH 7c (Pass 11094): at a mass level mu < 0 an empty sector is skipped, not fatal
      if (getenv("ORB_MASS_LEVEL") != NULL)
        continue;
      cout << "\n  Warning in bool COrbifold::Create() : Right-movers of the " << i << "-th sector have not been created. Return false." << endl;
      return false;
    }'''),
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
