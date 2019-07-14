#include "manager.h"

Manager::Manager() {
  int IDs[] = {1440, 1520, 1535, 1650, 1675, 1680, 1700, 1710, 1720, 1860, 1875,
               1880, 1895, 1900, 1990, 2000, 2060, 2100, 2120, 2190, 2250};

  for (int i = 0; i < kNumResonances; i++) {
    ids_[i] = IDs[i];
    resonances_[i].SetResonanceParameters(IDs[i]);
  }
}

TComplex Manager::F1(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F1(W, costh);
  }
  return sum;
}

TComplex Manager::F2(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F2(W, costh);
  }
  return sum;
}

TComplex Manager::F3(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F3(W, costh);
  }
  return sum;
}

TComplex Manager::F4(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F4(W, costh);
  }
  return sum;
}

// for debug
std::pair<TGraph, TGraph> Manager::F_web(int i, double W) {
  const int kNum = 200;
  double x[kNum], y0[kNum], y1[kNum];
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;
    if (i == 1) {
      y0[j] = -F1(W, x[j]).Re();
      y1[j] = -F1(W, x[j]).Im();
    } else if (i == 2) {
      y0[j] = -F2(W, x[j]).Re();
      y1[j] = -F2(W, x[j]).Im();
    } else if (i == 3) {
      y0[j] = -F3(W, x[j]).Re();
      y1[j] = -F3(W, x[j]).Im();
    } else if (i == 4) {
      y0[j] = -F4(W, x[j]).Re();
      y1[j] = -F4(W, x[j]).Im();
    }
  }
  TGraph tg0(kNum, x, y0);
  tg0.SetLineColor(2);
  TGraph tg1(kNum, x, y1);
  tg1.SetLineColor(2);
  return std::make_pair(tg0, tg1);
}

void Manager::MakeMultipoleRootFile(const char *filename) {
  const int kNumPts = 10000;
  double xx[kNumPts], yy[kNumResonances][4][kNumPts];  // Re_E, Im_E, Re_M, Im_M
  for (int i = 0; i < kNumPts; i++) {
    xx[i] = 1486.2 + 0.1 * i;
    for (int j = 0; j < kNumResonances; j++) {
      TComplex EE = resonances_[j].Multipole(xx[i], true);
      yy[j][0][i] = EE.Re();
      yy[j][1][i] = EE.Im();
      TComplex MM = resonances_[j].Multipole(xx[i], false);
      yy[j][2][i] = MM.Re();
      yy[j][3][i] = MM.Im();
    }
  }

  TGraph *tg[kNumResonances][4];  // Re_E, Im_E, Re_M, Im_M
  TString title[4] = {"Re(E)", "Im(E)", "Re(M)", "Im(M)"};
  for (int i = 0; i < kNumResonances; i++) {
    for (int j = 0; j < 4; j++) {
      tg[i][j] = new TGraph(kNumPts, xx, yy[i][j]);
      tg[i][j]->SetName(Form("tg%02d_%d", i, j));
      tg[i][j]->SetTitle(Form("N(%d)  %s", ids_[i], title[j].Data()));
    }
  }

  TFile *file0 = new TFile(filename, "recreate");
  for (int i = 0; i < kNumResonances; i++) {
    for (int j = 0; j < 4; j++) tg[i][j]->Write();
  }
  file0->Close();
}
