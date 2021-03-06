#ifndef BORN_H_
#define BORN_H_

#include <TComplex.h>

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
class Born {
 public:
  Born() {}

  void set_W(double W) { W_ = W; }
  void set_costh(double costh) { costh_ = costh; }

  constexpr static double mh = 547.862;
  constexpr static double mN = 938.272081;
  constexpr static double Wthr = 547.862 + 938.272081;
  constexpr static double MeVfm_inv = 197.3269788;
  constexpr static double e = 0.30282211986;  // sqrt(4 * pi / 137.036)
  constexpr static double kappa_p = 1.79;
  constexpr static double kappa_n = -1.91;
  constexpr static int e_p = 1;
  constexpr static int e_n = 0;

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
