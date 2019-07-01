#include <cstdio>
#include "hmaid.h"
#include <TGraph.h>
#include <TFile.h>

int main() {
  EtaMaid hmaid;

  const int num = 10000;
  double xx[num], yy0[num], yy1[num];
  for (int i = 0; i < num; i++) {
    xx[i] = 1486.2 + 0.1 * i;
    TComplex c = hmaid.Multipole(xx[i], false);
    yy0[i] = c.Re();
    yy1[i] = c.Im();
  }

  TGraph *tg0 = new TGraph(num, xx, yy0);
  TGraph *tg1 = new TGraph(num, xx, yy1);
  tg0->SetName("tg0");
  tg1->SetName("tg1");

  TFile *file0 = new TFile("rt/out.root", "recreate");
  tg0->Write();
  tg1->Write();
  file0->Close();

  return 0;
}
