#ifndef MANAGER_H_
#define MANAGER_H_

#include "resonance.h"
#include "born.h"
#include "regge.h"

#include "TGraph.h"
#include "TFile.h"

class Manager {
 public:
  Manager();

  TComplex F1(double W, double costh);
  TComplex F2(double W, double costh);
  TComplex F3(double W, double costh);
  TComplex F4(double W, double costh);
  // temporarily ..
  TComplex F1_res_born(double W, double costh);
  TComplex F2_res_born(double W, double costh);
  TComplex F3_res_born(double W, double costh);
  TComplex F4_res_born(double W, double costh);

  Resonance resonance(int i) const { return resonances_[i]; }

  // Following methods are for debugging.

  // Generates a root file containing multipoles.
  //
  // eg.
  //   manager.ResonanceMultipole();
  // in main.cc
  void ResonanceMultipole();

  // Generates a root file containing Resonance CGLN amplitudes.
  //
  // eg.
  //   manager.ResonanceCGLN(2000);
  // in main.cc
  void ResonanceCGLN(int W);

  // Generates a root file containing Born CGLN amplitudes.
  //
  // eg.
  //   manager.BornCGLN(2000);
  // in main.cc
  void BornCGLN(int W);

  // Generates a root file containing Regge CGLN amplitudes.
  //
  // eg.
  //   manager.ReggeCGLN(2000);
  // in main.cc
  void ReggeCGLN(int W);

  // Generates a root file containing Res+Born CGLN amplitudes.
  //
  // eg.
  //   manager.ResBornCGLN(2000);
  // in main.cc
  void ResBornCGLN(int W);

  // Generates a root file containing Res+Born+Regge CGLN amplitudes.
  //
  // eg.
  //   manager.FullCGLN(2000);
  // in main.cc
  void FullCGLN(int W);

  static const int kNumResonances = 21;

 private:
  int ids_[kNumResonances];
  Resonance resonances_[kNumResonances];
  Born born_;
  Regge regge_;
};

#endif  // MANAGER_H_
