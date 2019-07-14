#ifndef BORN_H_
#define BORN_H_

#include <TComplex.h>
#include <TGraph.h>

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
class Born {
 public:
  Born() {}

  void set_W(double W) { W_ = W; }
  void set_costh(double costh) { costh_ = costh; }

  static const double mh = 547.862;
  static const double mN = 938.272081;
  static const double Wthr = 547.862 + 938.272081;
  static const double MeVfm_inv = 197.3269788;
  static const double e = 0.30282211986;  // sqrt(4 * pi / 137.036)
  static const double kappa_p = 1.79;
  static const double kappa_n = -1.91;
  static const int e_p = 1;
  static const int e_n = 0;

  double F1(double W, double costh);  // [mfm]
  double F2(double W, double costh);  // [mfm]
  double F3(double W, double costh);  // [mfm]
  double F4(double W, double costh);  // [mfm]
  double g();
  double t();
  double u();
  double q();
  double E2();
  double C();
  double D();
  double PDK(double W, double m1, double m2);

 private:
  double W_;
  double costh_;
};

#endif  // BORN_H_
