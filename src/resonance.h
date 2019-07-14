#ifndef RESONANCE_H_
#define RESONANCE_H_

#include <TComplex.h>
#include <Math/SpecFunc.h>

// L. Tiator et al., Eur. Phys. J. A 54, 210 (2018)
class Resonance {
 public:
  Resonance() {}

  void set_W(double W) { W_ = W; }

  TComplex Multipole(double W, bool IsE);  // [mfm]

  // Eq. (33)
  TComplex t(TComplex tBW);

  // Eq. (34)
  TComplex tBW(double Mbar);

  // Eq. (36)
  std::pair<double, double> GammaTotal();

  // Eqs. (37)-(39)
  double PartialWidth(double beta, double m_meson, double m_baryon, bool pi2);

  // Eq. (45)
  double PartialWidthBelowThr(double g, double m_meson, double m_baryon);

  // Eq. (35)
  // widths.first: G_tot, widths.second: G_hN
  double BreitWignerFactor(std::pair<double, double> widths);

  // Eq. (40)
  double PhotonVertex();

  // PDK function to calculate k, q
  double PDK(double W, double m1, double m2);

  void SetResonanceParameters(int W);
  void A2M(double A12, double A32);

  // CGLN amplitudes
  TComplex F1(double W, double x);             // [mfm], x = costh
  TComplex F2(double W, double x);             // [mfm], x = costh
  TComplex F3(double W, double x);             // [mfm], x = costh
  TComplex F4(double W, double x);             // [mfm], x = costh
  double dP(unsigned int l, double x) const;   // P'_{l}(x)
  double ddP(unsigned int l, double x) const;  // P''_{l}(x)

  static const double mpi = 139.57061;  // charged
  static const double mK = 493.677;
  static const double mh = 547.862;
  static const double mw = 782.65;
  static const double mhp = 958.78;
  static const double mN = 938.272081;
  static const double mL = 1115.683;
  static const double mS = 1192.642;
  static const double MeVfm_inv = 197.3269788;

 private:
  double W_;
  // Resonance parameters
  unsigned int l_;
  unsigned int J21_;  // 2J+1
  int zeta_hN_;
  int zeta_hpN_;
  double M_R_;
  double G_R_;
  double b_piN_;
  double b_pipiN_;
  double b_hN_;
  double b_KL_;
  double b_KS_;
  double b_wN_;
  double b_hpN_;

  // Eq. (45), Table 5
  double g_hN_;
  double g_KS_;
  double g_wN_;
  double g_hpN_;

  double phi_;

  double Ebar_;
  double Mbar_;
};

#endif  // RESONANCE_H_
