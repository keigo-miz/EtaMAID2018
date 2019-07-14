#ifndef MANAGER_H_
#define MANAGER_H_

#include "resonance.h"

class Manager {
 public:
  Manager();

  Resonance resonance(int i) const { return resonances_[i]; }

  static const int kNumResonances = 21;

 private:
  Resonance resonances_[kNumResonances];
};

#endif  // MANAGER_H_
