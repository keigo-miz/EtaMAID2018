#include <cstdio>
#include "manager.h"
#include <TGraph.h>
#include <TFile.h>

int main() {
  Manager manager;

  manager.MakeMultipoleRootFile("rt/out.root");

  return 0;
}
