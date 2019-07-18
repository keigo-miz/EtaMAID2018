#include "manager.h"

int main() {
  Manager manager;
  manager.ResonanceMultipole();
  manager.ResonanceCGLN(1500);
  manager.ResonanceCGLN(2000);
  manager.ResonanceCGLN(2500);
  manager.ResonanceCGLN(3000);
  manager.BornCGLN(1500);
  manager.BornCGLN(2000);
  manager.BornCGLN(2500);
  manager.BornCGLN(3000);
  manager.ReggeCGLN(1500);
  manager.ReggeCGLN(2000);
  manager.ReggeCGLN(2500);
  manager.ReggeCGLN(3000);
  manager.FullCGLN(1500);
  manager.FullCGLN(2000);
  manager.FullCGLN(2500);
  manager.FullCGLN(3000);
  return 0;
}
