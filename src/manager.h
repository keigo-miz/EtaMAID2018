#ifndef MANAGER_H_
#define MANAGER_H_

#include "resonance.h"

#include "TGraph.h"
#include "TFile.h"

class Manager {
 public:
  Manager();

  TComplex F1(double W, double costh);
  TComplex F2(double W, double costh);
  TComplex F3(double W, double costh);
  TComplex F4(double W, double costh);

  Resonance resonance(int i) const { return resonances_[i]; }

  // Followings are methods for debugging.

  // Generates a root file which contains CGLN amplitudes.
  //
  // eg.
  //   manager.F_web(2000);
  // in main.cc
  void F_web(int W);

  // Generates a root file which contains 21*4 TGraph's.
  // Energy dependences of multipoles for each resonance.
  //
  // eg.
  //   manager.MakeMultipoleRootFile("rt/out.root");
  // in main.cc
  void MakeMultipoleRootFile(const char *filename);

  static const int kNumResonances = 21;

 private:
  int ids_[kNumResonances];
  Resonance resonances_[kNumResonances];
};

#endif  // MANAGER_H_
