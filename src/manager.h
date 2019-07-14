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

  std::pair<TGraph, TGraph> F_web(int i, double W);

  // Generates a root file which contains 21*4 TGraph's.
  // Energy dependences of multipoles for each resonance.
  void MakeMultipoleRootFile(const char *filename);

  Resonance resonance(int i) const { return resonances_[i]; }

  static const int kNumResonances = 21;

 private:
  int ids_[kNumResonances];
  Resonance resonances_[kNumResonances];
};

#endif  // MANAGER_H_
