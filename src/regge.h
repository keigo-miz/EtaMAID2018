#ifndef REGGE_H_
#define REGGE_H_

#include <TComplex.h>
#include <TMatrixD.h>
#include <TVectorD.h>

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
class Regge {
 public:
  Regge() : IsProton_(1) {}

  void set_W(double W) { W_ = W; }
  void set_costh(double costh) { costh_ = costh; }
  void set_IsProton(int IsProton) { IsProton_ = IsProton; }

  constexpr static double mh = 0.547862;
  constexpr static double mN = 0.938272081;
  constexpr static double Wthr = 0.547862 + 0.938272081;
  constexpr static double GeVmfm_inv = 197.3269788;
  constexpr static double e = 0.30282211986;  // sqrt(4 * pi / 137.036)

  // Table 7
  // electromagnetic couplings (\lambda)
  constexpr static double l_rho = 0.910;
  constexpr static double l_omega = 0.246;
  constexpr static double l_b1 = 1.0;  // By Kashevarov's mail
  // vector couplings to the nucleon (gv)
  constexpr static double gv_rho = 2.71;
  constexpr static double gv_omega = 14.2;
  // tensor couplings to the nucleon (gt)
  constexpr static double gt_rho = 11.4;  // By Kashevarov's mail
  constexpr static double gt_omega = 0.0;
  constexpr static double gt_b1 = -7.0;
  // coefficients for natural parity Regge cuts (c)
  constexpr static double c_pP = 4.64;   // rho-Pomeron
  constexpr static double c_pf = 3.10;   // rho-f2
  constexpr static double c_wP = -5.00;  // omega-Pomeron
  constexpr static double c_wf = 1.11;   // omega-f2
  // coefficients for unnatural parity Regge cuts (c~)
  constexpr static double ct_pP = 0.0;     // rho-Pomeron
  constexpr static double ct_pf = 0.245;   // rho-f2
  constexpr static double ct_wP = 0.0;     // omega-Pomeron
  constexpr static double ct_wf = -0.122;  // omega-f2
  // Regge cut parameters, d_c
  constexpr static double dc_pP = 12.1;  // rho-Pomeron
  constexpr static double dc_pf = 2.09;  // rho-f2 [Tiator2018 is wrong!]
  constexpr static double dc_wP = 12.1;  // omega-Pomeron [Tiator2018 is wrong!]
  constexpr static double dc_wf = 2.09;  // omega-f2

  // Reggeon and cut trajectories
  // Kashevarov 2017, TABLE II
  constexpr static double a0_rho = 0.477;
  constexpr static double ap_rho = 0.885;
  constexpr static double a0_omega = 0.434;
  constexpr static double ap_omega = 0.923;
  constexpr static double a0_b1 = -0.013;
  constexpr static double ap_b1 = 0.664;
  constexpr static double a0_h1 = -0.013;
  constexpr static double ap_h1 = 0.664;
  constexpr static double a0_rho2 = -0.235;
  constexpr static double ap_rho2 = 0.774;
  constexpr static double a0_omega2 = -0.235;
  constexpr static double ap_omega2 = 0.774;
  constexpr static double a0_f2 = 0.671;
  constexpr static double ap_f2 = 0.817;
  constexpr static double a0_P = 1.08;
  constexpr static double ap_P = 0.25;
  constexpr static double a0_pf = 0.148;
  constexpr static double ap_pf = 0.425;
  constexpr static double a0_wf = 0.106;
  constexpr static double ap_wf = 0.436;
  constexpr static double a0_pP = 0.557;
  constexpr static double ap_pP = 0.195;
  constexpr static double a0_wP = 0.514;
  constexpr static double ap_wP = 0.197;

  TComplex A1(double t);
  TComplex A2p(double t);
  TComplex A3(double t);
  TComplex A4(double t);
  TComplex F1(double W, double costh);
  TComplex F2(double W, double costh);
  TComplex F3(double W, double costh);
  TComplex F4(double W, double costh);

  // Converts invariant amplitudes {A_i} to CGLN amplitudes {F_i}.
  void A2F(double t);

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

  TComplex A1_;
  TComplex A2p_;
  TComplex A3_;
  TComplex A4_;
  TComplex F1_;
  TComplex F2_;
  TComplex F3_;
  TComplex F4_;

  int IsProton_;
};

#endif  // REGGE_H_
