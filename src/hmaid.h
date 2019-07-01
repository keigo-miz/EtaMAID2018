#ifndef HMAID_H_
#define HMAID_H_

#include <TComplex.h>
#include <Math/SpecFunc.h>

// L. Tiator et al.
// Eur. Phys. J. A (2018) 54: 210
class EtaMaid {
 public:
  EtaMaid() {}

  void set_W(double W) { W_ = W; }

  TComplex Multipole(double W, bool IsE);

  // Eq. (33)
  TComplex t(TComplex tBW);

  // Eq. (34)
  TComplex tBW(double Mbar);

  // Eq. (36)
  double GammaTotal();

  // Eqs. (37)-(39)
  double PartialWidth(double beta, double m_meson, double m_baryon, bool pi2);

  // Eq. (35)
  double BreitWignerFactor(double Gtot);

  // Eq. (40)
  double PhotonVertex();

  // PDK function to calculate k, q
  double PDK(double W, double m1, double m2) {
    return (TMath::Sqrt((W * W - (m1 + m2) * (m1 + m2)) *
                        (W * W - (m1 - m2) * (m1 - m2))) /
            2.0 / W);
  }

  void SetN1535Parameters();
  void SetN1880Parameters();

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
  double M_R_;
  double G_R_;
  double b_piN_;
  double b_pipiN_;
  double b_hN_;
  double b_KL_;
  double b_KS_;
  double b_wN_;
  double b_hpN_;
  double phi_;

  double Ebar_;
  double Mbar_;
};

#endif  // HMAID_H_
