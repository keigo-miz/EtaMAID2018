#include "born.h"

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)

// [mfm]
double Born::F1() {
  double F1 =
      (-e_p + (W_ - mN) * kappa_p / 2.0 / mN) * D() +
      (t() - mh * mh) * kappa_p / (2.0 * mN * (W_ - mN) * (u() - mN * mN));
  F1 *= g() * C();
  return (F1 * MeVfm_inv * 1000.0);
}

// [mfm]
double Born::F2() {
  double F2 =
      (e_p + (W_ + mN) * kappa_p / 2.0 / mN) * D() +
      (t() - mh * mh) * kappa_p / (2.0 * mN * (W_ + mN) * (u() - mN * mN));
  F2 *= g() * C() * q() / (E2() + mN);
  return (F2 * MeVfm_inv * 1000.0);
}

// [mfm]
double Born::F3() {
  double F3 = 2 * e_p * (W_ - mN) * D() / (t() - mh * mh) -
              kappa_p / mN / (u() - mN * mN);
  F3 *= g() * C() * q();
  return (F3 * MeVfm_inv * 1000.0);
}

// [mfm]
double Born::F4() {
  double F4 = -2 * e_p * (W_ + mN) * D() / (t() - mh * mh) -
              kappa_p / mN / (u() - mN * mN);
  F4 *= g() * C() * q() * q() / (E2() + mN);
  return (F4 * MeVfm_inv * 1000.0);
}

double Born::g() {
  return (TMath::Sqrt(0.063 * 4 * TMath::Pi()) * TMath::Power(Wthr / W_, 4.51));
}

double Born::t() {
  double k = PDK(W_, 0.0, mN);
  return (mh * mh + 2 * k * q() * costh_ -
          2 * k * TMath::Sqrt(q() * q() + mh * mh));
}

double Born::u() { return (2 * mN * mN + mh * mh - W_ * W_ - t()); }

double Born::q() { return PDK(W_, mh, mN); }

double Born::E2() { return TMath::Sqrt(q() * q() + mN * mN); }
double Born::C() {
  double k = PDK(W_, 0.0, mN);
  double E1 = TMath::Sqrt(k * k + mN * mN);
  return (-(W_ - mN) * e * TMath::Sqrt((E1 + mN) * (E2() + mN)) /
          (8.0 * TMath::Pi() * W_));
}

double Born::D() { return (1.0 / (W_ * W_ - mN * mN) + 1.0 / (u() - mN * mN)); }

// PDK function to calculate k, q
double Born::PDK(double W, double m1, double m2) {
  return (TMath::Sqrt((W * W - (m1 + m2) * (m1 + m2)) *
                      (W * W - (m1 - m2) * (m1 - m2))) /
          2.0 / W);
}
