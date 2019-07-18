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
  sum += born_.F1(W, costh);
  sum += regge_.F1(W, costh);
  return sum;
}

TComplex Manager::F2(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F2(W, costh);
  }
  sum += born_.F2(W, costh);
  sum += regge_.F2(W, costh);
  return sum;
}

TComplex Manager::F3(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F3(W, costh);
  }
  sum += born_.F3(W, costh);
  sum += regge_.F3(W, costh);
  return sum;
}

TComplex Manager::F4(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F4(W, costh);
  }
  sum += born_.F4(W, costh);
  sum += regge_.F4(W, costh);
  return sum;
}

TComplex Manager::F1_res_born(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F1(W, costh);
  }
  sum += born_.F1(W, costh);
  return sum;
}

TComplex Manager::F2_res_born(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F2(W, costh);
  }
  sum += born_.F2(W, costh);
  return sum;
}

TComplex Manager::F3_res_born(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F3(W, costh);
  }
  sum += born_.F3(W, costh);
  return sum;
}

TComplex Manager::F4_res_born(double W, double costh) {
  TComplex sum(0.0, 0.0);
  for (int i = 0; i < kNumResonances; i++) {
    sum += resonances_[i].F4(W, costh);
  }
  sum += born_.F4(W, costh);
  return sum;
}

// for debug
void Manager::ResonanceMultipole() {
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

  TFile *file0 = new TFile("rt/resonance_multipole.root", "recreate");
  for (int i = 0; i < kNumResonances; i++) {
    for (int j = 0; j < 4; j++) tg[i][j]->Write();
  }
  file0->Close();
}

void Manager::ResonanceCGLN(int W) {
  const int kNum = 200;
  double x[kNum], y[4][2][kNum];
  TComplex sum;
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;

    // F1
    sum = TComplex(0.0, 0.0);
    for (int k = 0; k < kNumResonances; k++) {
      sum += resonances_[k].F1((double)W, x[j]);
    }
    y[0][0][j] = -sum.Re();
    y[0][1][j] = -sum.Im();

    // F2
    sum = TComplex(0.0, 0.0);
    for (int k = 0; k < kNumResonances; k++) {
      sum += resonances_[k].F2((double)W, x[j]);
    }
    y[1][0][j] = -sum.Re();
    y[1][1][j] = -sum.Im();

    // F3
    sum = TComplex(0.0, 0.0);
    for (int k = 0; k < kNumResonances; k++) {
      sum += resonances_[k].F3((double)W, x[j]);
    }
    y[2][0][j] = -sum.Re();
    y[2][1][j] = -sum.Im();

    // F4
    sum = TComplex(0.0, 0.0);
    for (int k = 0; k < kNumResonances; k++) {
      sum += resonances_[k].F4((double)W, x[j]);
    }
    y[3][0][j] = -sum.Re();
    y[3][1][j] = -sum.Im();
  }
  TGraph *tg[4][2];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      tg[i][j] = new TGraph(kNum, x, y[i][j]);
      tg[i][j]->SetName(Form("tg%d_%d", i, j));
      tg[i][j]->SetLineColor(kRed);
    }
  }

  TFile *file0 = new TFile(Form("rt/resonance_cgln%04d.root", W), "recreate");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) tg[i][j]->Write();
  }
  file0->Close();
}

void Manager::BornCGLN(int W) {
  const int kNum = 200;
  double x[kNum], y[4][kNum];
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;
    y[0][j] = -born_.F1((double)W, x[j]);
    y[1][j] = -born_.F2((double)W, x[j]);
    y[2][j] = -born_.F3((double)W, x[j]);
    y[3][j] = -born_.F4((double)W, x[j]);
  }
  TGraph *tg[4];
  for (int i = 0; i < 4; i++) {
    tg[i] = new TGraph(kNum, x, y[i]);
    tg[i]->SetName(Form("tg%d", i));
    tg[i]->SetLineColor(kRed);
  }

  TFile *file0 = new TFile(Form("rt/born_cgln%04d.root", W), "recreate");
  for (int i = 0; i < 4; i++) {
    tg[i]->Write();
  }
  file0->Close();
}

void Manager::ResBornCGLN(int W) {
  const int kNum = 200;
  double x[kNum], y[4][2][kNum];
  TComplex tmp;
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;

    tmp = F1_res_born((double)W, x[j]);
    y[0][0][j] = -tmp.Re();
    y[0][1][j] = -tmp.Im();

    tmp = F2_res_born((double)W, x[j]);
    y[1][0][j] = -tmp.Re();
    y[1][1][j] = -tmp.Im();

    tmp = F3_res_born((double)W, x[j]);
    y[2][0][j] = -tmp.Re();
    y[2][1][j] = -tmp.Im();

    tmp = F4_res_born((double)W, x[j]);
    y[3][0][j] = -tmp.Re();
    y[3][1][j] = -tmp.Im();
  }
  TGraph *tg[4][2];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      tg[i][j] = new TGraph(kNum, x, y[i][j]);
      tg[i][j]->SetName(Form("tg%d_%d", i, j));
      tg[i][j]->SetLineColor(kRed);
    }
  }

  TFile *file0 = new TFile(Form("rt/res_born_cgln%04d.root", W), "recreate");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) tg[i][j]->Write();
  }
  file0->Close();
}

void Manager::FullCGLN(int W) {
  const int kNum = 200;
  double x[kNum], y[4][2][kNum];
  TComplex tmp;
  for (int j = 0; j < kNum; j++) {
    x[j] = -1.0 + 0.01 * j;

    tmp = F1((double)W, x[j]);
    y[0][0][j] = -tmp.Re();
    y[0][1][j] = -tmp.Im();

    tmp = F2((double)W, x[j]);
    y[1][0][j] = -tmp.Re();
    y[1][1][j] = -tmp.Im();

    tmp = F3((double)W, x[j]);
    y[2][0][j] = -tmp.Re();
    y[2][1][j] = -tmp.Im();

    tmp = F4((double)W, x[j]);
    y[3][0][j] = -tmp.Re();
    y[3][1][j] = -tmp.Im();
  }
  TGraph *tg[4][2];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      tg[i][j] = new TGraph(kNum, x, y[i][j]);
      tg[i][j]->SetName(Form("tg%d_%d", i, j));
      tg[i][j]->SetLineColor(kRed);
    }
  }

  TFile *file0 = new TFile(Form("rt/full_cgln%04d.root", W), "recreate");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) tg[i][j]->Write();
  }
  file0->Close();
}
