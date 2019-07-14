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
void Manager::F_web(int W) {
  const int kNum = 200;
  double x[kNum], y[4][2][kNum];
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;
    y[0][0][j] = -F1((double)W, x[j]).Re();
    y[0][1][j] = -F1((double)W, x[j]).Im();
    y[1][0][j] = -F2((double)W, x[j]).Re();
    y[1][1][j] = -F2((double)W, x[j]).Im();
    y[2][0][j] = -F3((double)W, x[j]).Re();
    y[2][1][j] = -F3((double)W, x[j]).Im();
    y[3][0][j] = -F4((double)W, x[j]).Re();
    y[3][1][j] = -F4((double)W, x[j]).Im();
  }
  TGraph *tg[4][2];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      tg[i][j] = new TGraph(kNum, x, y[i][j]);
      tg[i][j]->SetName(Form("tg%d_%d", i, j));
      tg[i][j]->SetLineColor(kRed);
    }
  }

  TFile *file0 = new TFile(Form("rt/%04d.root", W), "recreate");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) tg[i][j]->Write();
  }
  file0->Close();
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
      tg[i][j]->SetLineColor(kRed);
    }
  }

  TFile *file0 = new TFile(filename, "recreate");
  for (int i = 0; i < kNumResonances; i++) {
    for (int j = 0; j < 4; j++) tg[i][j]->Write();
  }
  file0->Close();
}
