#include "hmaid.h"

// L. Tiator et al.
// Eur. Phys. J. A (2018) 54: 210

// Eq. (33)
TComplex EtaMaid::t(TComplex tBW) {
  return (tBW * TComplex::Exp(TComplex(0.0, phi_ * TMath::Pi() / 180.0)));
}

// Eq. (34)
TComplex EtaMaid::tBW(double Mbar) {
  // C_{\eta N} = -1
  double Gtot = GammaTotal();
  double f_hN = BreitWignerFactor(Gtot);
  double f_gN = PhotonVertex();
  return (TComplex(-Mbar * M_R_ * Gtot * f_hN * f_gN, 0.0) / TComplex(M_R_ * M_R_ - W_ * W_, -M_R_ * Gtot));
}


// Eq. (36)
double EtaMaid::GammaTotal() {
  double Gamma = 0.0;
  if (W_ > 1077.84) {  // piN threshold
    Gamma += PartialWidth(b_piN_, mpi, mN, false);
  }
  if (W_ > 1217.41) {  // pipiN threshold
    Gamma += PartialWidth(b_pipiN_, 2 * mpi, mN, true);
  }
  if (W_ > 1486.13) {  // etaN threshold
    Gamma += PartialWidth(b_hN_, mh, mN, false);
  }
  if (W_ > 1609.36) {  // K Lambda threshold
    Gamma += PartialWidth(b_KL_, mK, mL, false);
  }
  if (W_ > 1686.32) {  // K Sigma threshold
    Gamma += PartialWidth(b_KS_, mK, mS, false);
  }
  if (W_ > 1720.92) {  // omega N threshold
    Gamma += PartialWidth(b_wN_, mw, mN, false);
  }
  if (W_ > 1896.05) {  // eta' N threshold
    Gamma += PartialWidth(b_hpN_, mhp, mN, false);
  }

  return Gamma;
}

// Eqs. (37)-(39)
double EtaMaid::PartialWidth(double beta, double m_meson, double m_baryon, bool pi2) {
  if (beta < 1.0e-6) return 0.0;
  const double X = 450.0;  // cut-off [MeV]

  double q = PDK(W_, m_meson, m_baryon);
  double qR = PDK(M_R_, m_meson, m_baryon);

  double Gamma = beta * G_R_;
  unsigned int effective_l = (pi2 ? l_ + 2 : l_);
  for (unsigned int i = 0; i < 2 * effective_l + 1; i++) {
    Gamma *= q / qR;
  }
  for (unsigned int i = 0; i < effective_l; i++) {
    Gamma *= (X * X + qR * qR) / (X * X + q * q);
  }

  return Gamma;
}

// Eq. (35)
double EtaMaid::BreitWignerFactor(double Gtot) {
  double k = PDK(W_, 0.0, mN);
  double q = PDK(W_, mh, mN);

  return (zeta_hN_ * TMath::Sqrt(1.0 / J21_ / TMath::Pi() * k / q * mN / M_R_ * PartialWidth(b_hN_, mh, mN, false) / Gtot / Gtot));
}

// Eq. (40)
double EtaMaid::PhotonVertex() {
  // X_{\gamma} = 0
  double k = PDK(W_, 0.0, mN);
  double kR = PDK(M_R_, 0.0, mN);

  return (kR * kR / k / k);
}

// Table 5 and 6
void EtaMaid::SetN1535Parameters() {
  l_ = 0;         // S-wave
  J21_ = 2;       // 2J+1 = 1
  zeta_hN_ = +1;
  M_R_ = 1521.7;  // [MeV]
  G_R_ = 174.7;   // [MeV]
  b_piN_ = 0.520;
  b_pipiN_ = 0.136;
  b_hN_ = 0.344;
  b_KL_ = 0.0;
  b_KS_ = 0.0;
  b_wN_ = 0.0;
  b_hpN_ = 0.0;
  phi_ = 29.0;  // [deg]
}
