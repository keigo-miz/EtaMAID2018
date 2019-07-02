#include <cstdio>
#include "hmaid.h"
#include <TGraph.h>
#include <TFile.h>

int main() {
  const int kNumRes = 21;
  int IDs[] = {1440, 1520, 1535, 1650, 1675, 1680, 1700, 1710, 1720, 1860, 1875,
               1880, 1895, 1900, 1990, 2000, 2060, 2100, 2120, 2190, 2250};

  EtaMaid hmaid[kNumRes];
  for (int i = 0; i < kNumRes; i++) hmaid[i].SetResonanceParameters(IDs[i]);

  const int kNumPts = 10000;
  double xx[kNumPts], yy[kNumRes][4][kNumPts];  // Re_E, Im_E, Re_M, Im_M
  for (int i = 0; i < kNumPts; i++) {
    xx[i] = 1486.2 + 0.1 * i;
    for (int j = 0; j < kNumRes; j++) {
      TComplex EE = hmaid[j].Multipole(xx[i], true);
      yy[j][0][i] = EE.Re();
      yy[j][1][i] = EE.Im();
      TComplex MM = hmaid[j].Multipole(xx[i], false);
      yy[j][2][i] = MM.Re();
      yy[j][3][i] = MM.Im();
    }
  }

  TGraph *tg[kNumRes][4];  // Re_E, Im_E, Re_M, Im_M
  for (int i = 0; i < kNumRes; i++) {
    for (int j = 0; j < 4; j++) {
      tg[i][j] = new TGraph(kNumPts, xx, yy[i][j]);
      tg[i][j]->SetName(Form("tg%02d_%d", i, j));
    }
  }

  TFile *file0 = new TFile("rt/out.root", "recreate");
  for (int i = 0; i < kNumRes; i++) {
    for (int j = 0; j < 4; j++) tg[i][j]->Write();
  }
  file0->Close();

  return 0;
}
