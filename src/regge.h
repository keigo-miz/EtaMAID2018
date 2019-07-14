#ifndef REGGE_H_
#define REGGE_H_

#include <TComplex.h>
#include <TMatrixD.h>
#include <TVectorD.h>

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
class Regge {
 public:
  Regge() {}

  void set_W(double W) { W_ = W; }
  void set_costh(double costh) { costh_ = costh; }
  void set_IsProton(int IsProton) { IsProton_ = IsProton; }

  static const double mh = 0.547862;
  static const double mN = 0.938272081;
  static const double Wthr = 0.547862 + 0.938272081;
  static const double GeVmfm_inv = 197.3269788;
  static const double e = 0.30282211986;  // sqrt(4 * pi / 137.036)

  // Table 7
  // electromagnetic couplings (\lambda)
  static const double l_rho = 0.910;
  static const double l_omega = 0.246;
  static const double l_b1 = 0.1;
  // vector couplings to the nucleon (gv)
  static const double gv_rho = 2.71;
  static const double gv_omega = 14.2;
  // tensor couplings to the nucleon (gt)
  static const double gt_rho = 4.20;
  static const double gt_omega = 0.0;
  static const double gt_b1 = -7.0;
  // coefficients for natural parity Regge cuts (c)
  static const double c_pP = 4.64;   // rho-Pomeron
  static const double c_pf = 3.10;   // rho-f2
  static const double c_wP = -5.00;  // omega-Pomeron
  static const double c_wf = 1.11;   // omega-f2
  // coefficients for unnatural parity Regge cuts (c~)
  static const double ct_pP = 0.0;     // rho-Pomeron
  static const double ct_pf = 0.245;   // rho-f2
  static const double ct_wP = 0.0;     // omega-Pomeron
  static const double ct_wf = -0.122;  // omega-f2
  // Regge cut parameters, d_c
  static const double dc_pP = 12.1;  // rho-Pomeron
  static const double dc_pf = 2.09;  // rho-f2 [Tiator2018 is wrong!]
  static const double dc_wP = 12.1;  // omega-Pomeron [Tiator2018 is wrong!]
  static const double dc_wf = 2.09;  // omega-f2

  // Reggeon and cut trajectories
  // Kashevarov 2017, TABLE II
  static const double a0_rho = 0.477;
  static const double ap_rho = 0.885;
  static const double a0_omega = 0.434;
  static const double ap_omega = 0.923;
  static const double a0_b1 = -0.013;
  static const double ap_b1 = 0.664;
  static const double a0_h1 = -0.013;
  static const double ap_h1 = 0.664;
  static const double a0_rho2 = -0.235;
  static const double ap_rho2 = 0.774;
  static const double a0_omega2 = -0.235;
  static const double ap_omega2 = 0.774;
  static const double a0_f2 = 0.671;
  static const double ap_f2 = 0.817;
  static const double a0_P = 1.08;
  static const double ap_P = 0.25;
  static const double a0_pf = 0.148;
  static const double ap_pf = 0.425;
  static const double a0_wf = 0.106;
  static const double ap_wf = 0.436;
  static const double a0_pP = 0.557;
  static const double ap_pP = 0.195;
  static const double a0_wP = 0.514;
  static const double ap_wP = 0.197;

  TComplex A1(double t);
  TComplex A2(double t);
  TComplex A2p(double t);
  TComplex A3(double t);
  TComplex A4(double t);
  double ForFit(double *par);

  TComplex D(double s, double t, double a0, double ap, double S);
  TComplex Dcut(double s, double t, double a0, double ap, double d_c);

  // L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
  // Appendix B
  // V.L. Kashevarov et al., Phys. Rev. C 96, 035207 (2017)
  // APPENDIX B
  TMatrixD ConvMatrix();

  // damping factor
  // L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
  // Eq. (29)
  double Fd(double W);

  double t();
  double PDK(double W, double m1, double m2);

 private:
  double W_;
  double costh_;

  int IsProton_;
};

#endif  // REGGE_H_
