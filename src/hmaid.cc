#include "hmaid.h"

// L. Tiator et al.
// Eur. Phys. J. A (2018) 54: 210

TComplex EtaMaid::Multipole(double W, bool IsE) {
  set_W(W);
  double tmp_bar = (IsE ? Ebar_ : Mbar_);
  return (t(tBW(tmp_bar)) * MeVfm_inv / TMath::Sqrt(1000.0));
}

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
  return (TComplex(-Mbar * M_R_ * Gtot * f_hN * f_gN, 0.0) /
          TComplex(M_R_ * M_R_ - W_ * W_, -M_R_ * Gtot));
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

    // Eq. (45) for below threshold
    if (g_hN_ > 0.001) Gamma += PartialWidthBelowThr(g_hN_, mh, mN);
  }
  if (W_ > 1609.36) {  // K Lambda threshold
    Gamma += PartialWidth(b_KL_, mK, mL, false);
  }
  if (W_ > 1686.32) {  // K Sigma threshold
    Gamma += PartialWidth(b_KS_, mK, mS, false);

    // Eq. (45) for below threshold
    if (g_KS_ > 0.001) Gamma += PartialWidthBelowThr(g_KS_, mK, mS);
  }
  if (W_ > 1720.92) {  // omega N threshold
    Gamma += PartialWidth(b_wN_, mw, mN, false);

    // Eq. (45) for below threshold
    if (g_wN_ > 0.001) Gamma += PartialWidthBelowThr(g_wN_, mw, mN);
  }
  if (W_ > 1896.05) {  // eta' N threshold
    Gamma += PartialWidth(b_hpN_, mhp, mN, false);

    // Eq. (45) for below threshold
    if (g_hpN_ > 0.001) Gamma += PartialWidthBelowThr(g_hpN_, mhp, mN);
  }

  return Gamma;
}

// Eqs. (37)-(39)
double EtaMaid::PartialWidth(double beta, double m_meson, double m_baryon,
                             bool pi2) {
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

// Eq. (45)
double EtaMaid::PartialWidthBelowThr(double g, double m_meson, double m_baryon) {
  const double X = 450.0;  // cut-off [MeV]
  const double q = PDK(W_, m_meson, m_baryon);

  double Gamma = g * g * q;
  for (unsigned int i = 0; i < l_; i++) {
    Gamma *= q * q / (X * X + q * q);
  }

  return Gamma;
}

// Eq. (35)
double EtaMaid::BreitWignerFactor(double Gtot) {
  double k = PDK(W_, 0.0, mN);
  double q = PDK(W_, mh, mN);

  return (zeta_hN_ *
          TMath::Sqrt(1.0 / J21_ / TMath::Pi() * k / q * mN / M_R_ *
                      PartialWidth(b_hN_, mh, mN, false) / Gtot / Gtot));
}

// Eq. (40)
double EtaMaid::PhotonVertex() {
  // X_{\gamma} = 0
  double k = PDK(W_, 0.0, mN);
  double kR = PDK(M_R_, 0.0, mN);

  return (kR * kR / k / k);
}

// PDK function to calculate k, q
double EtaMaid::PDK(double W, double m1, double m2) {
  return (TMath::Sqrt((W * W - (m1 + m2) * (m1 + m2)) * (W * W - (m1 - m2) * (m1 - m2))) / 2.0 / W);
}

// Table 5 and 6
void EtaMaid::SetResonanceParameters(int W) {
  // Reads multipoles.
  FILE *fp = fopen("dat/Table5.dat", "r");
  char buf[256];
  int M0, J0, zhN0, zhpN0;
  unsigned int l0;
  double MR0, GR0, bpiN0, bpi2N0, bhN0, bKL0, bKS0, bwN0, bhpN0, ghpN0;
  while (fgets(buf, sizeof(buf), fp) != nullptr) {
    if (buf[0] == '#') continue;
    sscanf(buf, "%d %d %u %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &M0, &J0, &l0, &zhN0, &zhpN0, &MR0, &GR0, &bpiN0, &bpi2N0, &bhN0, &bKL0, &bKS0, &bwN0, &bhpN0, &ghpN0);
    if (M0 == W) {
      J21_ = J0;
      l_ = l0;
      zeta_hN_ = zhN0;
      zeta_hpN_ = zhpN0;
      M_R_ = MR0;
      G_R_ = GR0;
      b_piN_ = bpiN0 * 0.01;
      b_pipiN_ = bpi2N0 * 0.01;
      b_hN_ = bhN0 * 0.01;
      b_KL_ = bKL0 * 0.01;
      b_KS_ = bKS0 * 0.01;
      b_wN_ = bwN0 * 0.01;
      b_hpN_ = bhpN0 * 0.01;
      g_hpN_ = ghpN0;
      break;
    }
  }
  fclose(fp);
  fp = fopen("dat/Table6.dat", "r");
  double pA12, pA32, nA12, nA32, phhp0, phhn0, phhpp0, phhpn0;
  while (fgets(buf, sizeof(buf), fp) != nullptr) {
    if (buf[0] == '#') continue;
    sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &M0, &pA12, &pA32, &nA12, &nA32, &phhp0, &phhn0, &phhpp0, &phhpn0);
    if (M0 == W) {
      phi_ = phhp0;
      A2M(pA12, pA32);
      break;
    }
  }
  fclose(fp);

  // Eq. (45), Table 5
  g_hN_ = 0.0;
  g_KS_ = 0.0;
  g_wN_ = 0.0;

  if (W == 1880) {
    // P-wave, 1-
    l_ = 1;    // P-wave
    J21_ = 2;  // 2J+1 = 1
    zeta_hN_ = +1;
    M_R_ = 1882.1;  // [MeV]
    G_R_ = 90.0;   // [MeV]
    b_piN_ = 0.060;
    b_pipiN_ = 0.746;
    b_hN_ = 0.0043;
    b_KL_ = 0.020;
    b_KS_ = 0.170;
    g_hpN_ = 0.400;
    phi_ = 84.9;  // [deg]

    // Table 2
    // Ebar = 0.0
    // Mbar = A_{1/2}
    double A12 = 60.4;
    Ebar_ = 0.0;
    Mbar_ = A12;
  } else if (W == 1520) {
    // D-wave, 2-
    l_ = 2;    // D-wave
    J21_ = 4;  // 2J+1 = 1
    zeta_hN_ = +1;
    M_R_ = 1520.0;  // [MeV]
    G_R_ = 100.0;   // [MeV]
    b_piN_ = 0.61;
    b_pipiN_ = 0.389;
    b_hN_ = 0.0008;
    phi_ = 55.3;  // [deg]

    // Table 2
    // Ebar = -(1/2)(sqrt(3)*A32 + A12)
    // Mbar = -(1/2)(A32 / sqrt(3) - A12)
    double A12 = -39.7;
    double A32 = 116.8;
    Ebar_ = -(TMath::Sqrt(3.0) * A32 + A12) / 2.0;
    Mbar_ = -(A32 / TMath::Sqrt(3.0) - A12) / 2.0;
  } else if (W == 1895) {
    // S-wave, 0+
    l_ = 0;    // S-wave
    J21_ = 2;  // 2J+1 = 1
    zeta_hN_ = +1;
    M_R_ = 1894.4;  // [MeV]
    G_R_ = 70.7;   // [MeV]
    b_piN_ = 0.025;
    b_pipiN_ = 0.632;
    b_hN_ = 0.0327;
    b_KL_ = 0.180;
    b_KS_ = 0.130;
    g_hpN_ = 0.405;
    phi_ = 51.5;  // [deg]

    // Table 2
    // Ebar = -A_{1/2}
    // Mbar = 0.0
    double A12 = -32.0;
    Ebar_ = -A12;
    Mbar_ = 0.0;
  }
}

// Table 2
void EtaMaid::A2M(double A12, double A32) {
  if (l_ == 0) {      // S11 wave
    Ebar_ = -A12;
    Mbar_ = 0.0;
  } else if (l_ == 1) {
    if (J21_ == 2) {  // P11 wave
      Ebar_ = 0.0;
      Mbar_ = A12;
    } else {          // P13 wave
      Ebar_ = 0.5 * (A32 / TMath::Sqrt(3.0) - A12);
      Mbar_ = -0.5 * (TMath::Sqrt(3.0) * A32 + A12);
    }
  } else if (l_ == 2) {
    if (J21_ == 4) {  // D13 wave
      Ebar_ = -0.5 * (TMath::Sqrt(3.0) * A32 + A12);
      Mbar_ = -0.5 * (A32 / TMath::Sqrt(3.0) - A12);
    } else {          // D15 wave
      Ebar_ = (1.0 / 3.0) * (A32 / TMath::Sqrt2() - A12);
      Mbar_ = -(1.0 / 3.0) * (TMath::Sqrt2() * A32 + A12);
    }
  } else if (l_ == 3) {
    if (J21_ == 6) {  // F15 wave
      Ebar_ = -(1.0 / 3.0) * (TMath::Sqrt2() * A32 + A12);
      Mbar_ = -(1.0 / 3.0) * (A32 / TMath::Sqrt2() - A12);
    } else {          // F17 wave
      Ebar_ = (1.0 / 4.0) * (TMath::Sqrt(3.0 / 5.0) * A32 - A12);
      Mbar_ = -(1.0 / 4.0) * (TMath::Sqrt(5.0 / 3.0) * A32 + A12);
    }
  } else if (l_ == 4) {
    if (J21_ == 8) {  // G17 wave
      Ebar_ = -(1.0 / 4.0) * (TMath::Sqrt(5.0 / 3.0) * A32 + A12);
      Mbar_ = -(1.0 / 4.0) * (TMath::Sqrt(3.0 / 5.0) * A32 - A12);
    } else {          // G19 wave
      Ebar_ = (1.0 / 5.0) * (TMath::Sqrt(2.0 / 3.0) * A32 - A12);
      Mbar_ = -(1.0 / 5.0) * (TMath::Sqrt(3.0 / 2.0) * A32 + A12);
    }
  }
}
