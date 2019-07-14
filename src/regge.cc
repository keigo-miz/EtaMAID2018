#include "regge.h"

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)

TMatrixD Regge::ConvMatrix() {
  TMatrixD mat(4, 4);
  mat.Zero();

  double k = PDK(W_, 0.0, mN);
  double q = PDK(W_, mh, mN);
  double Ei = TMath::Sqrt(k * k + mN * mN);
  double Ef = TMath::Sqrt(q * q + mN * mN);
  double nu_B = (t() - mh * mh) / 4.0 / mN;
  double piW8 = 8.0 * TMath::Pi() * W_;
  double cF1 = (W_ - mN) * TMath::Sqrt((Ei + mN) * (Ef + mN)) / piW8;
  double cF2 = (W_ + mN) * q * TMath::Sqrt((Ei - mN) / (Ef + mN)) / piW8;
  double cF3 = (W_ + mN) * q * TMath::Sqrt((Ei - mN) * (Ef + mN)) / piW8;
  double cF4 = (W_ - mN) * q * q * TMath::Sqrt((Ei + mN) / (Ef + mN)) / piW8;

  mat[0][0] = cF1;
  mat[0][2] = -cF1 * 2.0 * mN * nu_B / (W_ - mN);
  mat[0][3] = cF1 * (W_ - mN + 2.0 * mN * nu_B / (W_ - mN));
  mat[1][0] = -cF2;
  mat[1][2] = -cF2 * 2.0 * mN * nu_B / (W_ + mN);
  mat[1][3] = cF2 * (W_ + mN + 2.0 * mN * nu_B / (W_ + mN));
  mat[2][1] = cF3 * (W_ - mN);
  mat[2][2] = cF3;
  mat[2][3] = -cF3;
  mat[3][1] = -cF4 * (W_ + mN);
  mat[3][2] = cF4;
  mat[3][3] = -cF4;

  return mat;
}

// damping factor
// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
// Eq. (29)
double Regge::Fd(double W) {
  const double LambdaR = 0.974;
  return (1.0 - TMath::Exp(-(W - Wthr) / LambdaR));
}

TComplex Regge::A1(double t) {
  double s = W_ * W_;
  TComplex Drho = D(s, t, a0_rho, ap_rho, -1.0);
  TComplex DrhoP = Dcut(s, t, a0_pP, ap_pP, dc_pP);
  TComplex Drhof = Dcut(s, t, a0_pf, ap_pf, dc_pf);
  TComplex rho_tmp = Drho + c_pP * DrhoP + c_pf * Drhof;
  rho_tmp = l_rho * gt_rho * rho_tmp;
  if (IsProton_ == 0) rho_tmp = -1.0 * rho_tmp;
  return (Fd(W_) * t * e / 2.0 / mh / mN * rho_tmp);
}

TComplex Regge::A2p(double t) {
  double s = W_ * W_;
  TComplex Db1 = D(s, t, a0_b1, ap_b1, -1.0);
  TComplex Drhof = Dcut(s, t, a0_pf, ap_pf, dc_pf);

  TComplex first = l_b1 * gt_b1 * Db1;
  TComplex second = l_rho * gt_rho * ct_pf * Drhof;

  return (-Fd(W_) * t * e / 2.0 / mh / mN * (first + second));
}

TComplex Regge::A2(double t) { return (1.0 / t) * (A2p(t) - A1(t)); }

TComplex Regge::A3(double t) {
  double s = W_ * W_;
  TComplex Drhof = Dcut(s, t, a0_pf, ap_pf, dc_pf);
  TComplex Domegaf = Dcut(s, t, a0_wf, ap_wf, dc_wf);

  TComplex first = l_rho * gv_rho * ct_pf * Drhof;
  if (IsProton_ == 0) first = -1.0 * first;
  TComplex second = l_omega * gv_omega * ct_wf * Domegaf;

  return (Fd(W_) * e / mh * (first + second));
}

TComplex Regge::A4(double t) {
  double s = W_ * W_;
  TComplex Drho = D(s, t, a0_rho, ap_rho, -1.0);
  TComplex DrhoP = Dcut(s, t, a0_pP, ap_pP, dc_pP);
  TComplex Drhof = Dcut(s, t, a0_pf, ap_pf, dc_pf);
  TComplex rho_tmp = Drho + c_pP * DrhoP + c_pf * Drhof;
  rho_tmp = l_rho * gv_rho * rho_tmp;
  if (IsProton_ == 0) rho_tmp = -1.0 * rho_tmp;

  TComplex Domega = D(s, t, a0_omega, ap_omega, -1.0);
  TComplex DomegaP = Dcut(s, t, a0_wP, ap_wP, dc_wP);
  TComplex Domegaf = Dcut(s, t, a0_wf, ap_wf, dc_wf);
  TComplex omega_tmp = Domega + c_wP * DomegaP + c_wf * Domegaf;
  omega_tmp = l_omega * gv_omega * omega_tmp;

  return (-Fd(W_) * e / mh * (rho_tmp + omega_tmp));
}

double Regge::ForFit(double *par) {
  double s = W_ * W_;
  double tt = t();
  TComplex Drhof = Dcut(s, tt, a0_pf, ap_pf, par[0]);
  TComplex Domegaf = Dcut(s, tt, a0_wf, ap_wf, par[1]);

  TComplex first = l_rho * gv_rho * par[2] * Drhof;
  TComplex second = l_omega * gv_omega * par[3] * Domegaf;

  TComplex tmp = Fd(W_) * e / mh * (first + second);
  return -tmp.Im();
}

TComplex Regge::D(double s, double t, double a0, double ap, double S) {
  const double s0 = 1.0;  // GeV^2
  double alpha = a0 + ap * t;
  return TMath::Power(s / s0, alpha - 1.0) * ap * TMath::Gamma(1.0 - alpha) /
         2.0 * (TComplex(S, 0.0) +
                TComplex::Exp(TComplex(0.0, -TMath::Pi() * alpha)));
}

TComplex Regge::Dcut(double s, double t, double a0, double ap, double d_c) {
  const double s0 = 1.0;  // GeV^2
  double alpha = a0 + ap * t;
  return TMath::Power(s / s0, alpha - 1.0) *
         TComplex::Exp(TComplex(d_c * t, -TMath::Pi() * alpha / 2.0));
}

double Regge::t() {
  double k = PDK(W_, 0.0, mN);
  double q = PDK(W_, mh, mN);
  return (mh * mh + 2 * k * q * costh_ - 2 * k * TMath::Sqrt(q * q + mh * mh));
}

// PDK function to calculate k, q
double Regge::PDK(double W, double m1, double m2) {
  return (TMath::Sqrt((W * W - (m1 + m2) * (m1 + m2)) *
                      (W * W - (m1 - m2) * (m1 - m2))) /
          2.0 / W);
}
