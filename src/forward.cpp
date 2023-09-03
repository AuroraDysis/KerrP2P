#include "forward.h"


//  Real x3(Real r) const {
//    return (A * (r - r1) - B * (r - r2)) / (A * (r - r1) + B * (r - r2));
//  }
//
//  Real F3(Real r) const {
//    return 1 / mp::sqrt(A * B) * boost::math::ellint_1(mp::acos(x3(r)), k3);
//  }
//

//
//  Real Pi13(Real r) const {
//    return (2 * (r2 - r1) * mp::sqrt(A * B)) / (B * B - A * A) * R1(alpha_0, mp::acos(x3(r)), k3);
//  }
//
//  Real Pi23(Real r) const {
//    return square(((2 * (r2 - r1) * mp::sqrt(A * B)) / (B * B - A * A))) *
//           R2(alpha_0, mp::acos(x3(r)), k3);
//  }
//
//  Real I_0(Real r) const {
//    return F3(r);
//  }
//
//  Real I_1(Real r) const {
//    return ((B * r2 + A * r1) / (B + A)) * F3(r) + Pi13(r);
//  }
//
//  Real I_2(Real r) const {
//    return square(((B * r2 + A * r1) / (B + A))) * F3(r) +
//           2 * ((B * r2 + A * r1) / (B + A)) * Pi13(r) +
//           mp::sqrt(A * B) * Pi23(r);
//  }
//
//  Real I_p(Real r) const {
//    return -1 / (B * (rp - r2) + A * (rp - r1)) * ((B + A) * F3(r) +
//                                                   (2 * (r2 - r1) * mp::sqrt(A * B)) /
//                                                   (B * (rp - r2) - A * (rp - r1)) *
//                                                   R1(alpha_p, mp::acos(x3(r)), k3));
//  }
//
//  Real I_m(Real r) const {
//    return -1 / (B * (rm - r2) + A * (rm - r1)) * ((B + A) * F3(r) +
//                                                   (2 * (r2 - r1) * mp::sqrt(A * B)) /
//                                                   (B * (rm - r2) - A * (rm - r1)) *
//                                                   R1(alpha_m, mp::acos(x3(r)), k3));
//  }


//
//  Real I_r(Real r) const {
//    return I_0(r);
//  }
//
//  Real I_phi(Real r) const {
//    return (2 * a) / (rp - rm) * ((rp - (a * lambda) / 2) * I_p(r) - (rm - (a * lambda) / 2) * I_m(r));
//  }
//
//  Real I_t(Real r) const {
//    return 4 / (rp - rm) * (rp * (rp - (a * lambda) / 2) * I_p(r) - rm * (rm - (a * lambda) / 2) * I_m(r))
//           + 4 * I_0(r) + 2 * I_1(r) + I_2(r);
//  }

#define square(x) ((x) * (x))

void IIntegral2::pre_calc() {
  const Real &rp = data.rp;
  const Real &rm = data.rm;

  const Real &r1 = data.r1;
  const Real &r2 = data.r2;
  const Real &r3 = data.r3;
  const Real &r4 = data.r4;
  const Real &r_s = data.r_s;

  asin_x2_rs = mp::sqrt(((r1 - r3) * (r_s - r4)) / ((r_s - r3) * (r1 - r4)));
  asin_x2_inf = mp::sqrt((r1 - r3) / (r1 - r4));
  k = ((-r2 + r3) * (-r1 + r4)) / ((r1 - r3) * (r2 - r4));

  E2_Numerator = mp::sqrt((r1 - r3) * (r2 - r4));
  F2_Denominator = 1 / E2_Numerator;
  Pi_p2_T1 = ((r1 - r4) * (r3 - rp)) / ((r1 - r3) * (r4 - rp));
  Pi_m2_T1 = ((r1 - r4) * (r3 - rm)) / ((r1 - r3) * (r4 - rm));
  Pi_p2_m2_Numerator = 2 * (-r3 + r4);
  Pi_p2_Denominator = 1 / (E2_Numerator * (-r3 + rp) * (-r4 + rp));
  Pi_m2_Denominator = 1 / (E2_Numerator * (-r3 + rm) * (-r4 + rm));
}

void IIntegral2::calc_x(std::array<Real, 3>& integral, const Real &asin_x2) {
  const Real &a = data.a;
  const Real &lambda = data.lambda;
  const Real &rp = data.rp;
  const Real &rm = data.rm;
  const Real &r3 = data.r3;

  F2 = (2 * boost::math::ellint_1(asin_x2, k)) * F2_Denominator;
  E2 = E2_Numerator * boost::math::ellint_2(asin_x2, k);
  Pi_p2 = (Pi_p2_m2_Numerator * boost::math::ellint_3(Pi_p2_T1, asin_x2, k)) * Pi_p2_Denominator;
  Pi_m2 = (Pi_p2_m2_Numerator * boost::math::ellint_3(Pi_m2_T1, asin_x2, k)) * Pi_m2_Denominator;
  I_p = F2 / (r3 - rp) - Pi_p2;
  I_m = F2 / (r3 - rm) - Pi_m2;

  // I_r
  integral[0] = F2;
  // I_phi
  integral[1] = (a * (-2 * rp * I_p + a * I_p * lambda + (2 * rm - a * lambda) * I_m)) / (rm - rp);
}

void IIntegral2::calc(bool is_plus) {
  pre_calc();
  calc_x(integral_rs, asin_x2_rs);
  calc_x(integral_inf, asin_x2_inf);

  auto &radial_integrals = data.radial_integrals;

  for (int i = 0; i < 2; ++i) {
    if (is_plus) {
      radial_integrals[i] = integral_inf[i] + integral_rs[i];
    } else {
      radial_integrals[i] = integral_inf[i] - integral_rs[i];
    }
  }
}

void IIntegral3::pre_calc() {
  const Real &rp = data.rp;
  const Real &rm = data.rm;
  const Real &r_s = data.r_s;
  // two real roots, both inside horizon, r_1 < r_2 < r_- < r_+ and r_3 = conj(r_4)
  const Real &r1 = data.r1;
  const Real &r2 = data.r2;
  const Complex &r3 = data.r3_c;
  // const Complex &r4 = data.r4_c;

  r34_re = mp::real(r3);
  r34_im = mp::imag(r3);

  // radial coeffs
  A = mp::sqrt(square(r34_im) + square(r34_re - r2));
  B = mp::sqrt(square(r34_im) + square(r34_re - r1));

  k3 = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
  // alpha_0 = (B + A) / (B - A);
  alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
  alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));
  alpha_p2 = square(alpha_p);
  alpha_m2 = square(alpha_m);

  acos_x3_rs = mp::acos(-1 + (2 * A * (r_s - r1)) / (A * (r_s - r1) + B * (r_s - r2)));
  acos_x3_inf = mp::acos((A - B) / (A + B));
}


Real IIntegral3::R1(const Real &acos_x3, const Real &alpha, const Real &alpha2) const {
  return 1 / (1 - alpha2) *
         (boost::math::ellint_3(alpha2 / (alpha2 - 1), acos_x3, k3) -
          alpha * ((mp::sqrt((-1 + alpha2) /
                             (alpha2 + k3 - alpha2 * k3)) *
                    mp::log(mp::abs((mp::sin(acos_x3) +
                                     mp::sqrt((-1 + alpha2) /
                                              (alpha2 + k3 - alpha2 * k3)) *
                                     mp::sqrt(1 - k3 * square(mp::sin(acos_x3)))) /
                                    (-mp::sin(acos_x3) + mp::sqrt((-1 + alpha2) /
                                                                  (alpha2 + k3 - alpha2 * k3)) *
                                                         mp::sqrt(1 - k3 * square(mp::sin(acos_x3))))))) *
                   half<Real>()));
}

void IIntegral3::calc_x(std::array<Real, 3>& integral, const Real &acos_x3) {
  const Real &a = data.a;
  const Real &lambda = data.lambda;
  const Real &rp = data.rp;
  const Real &rm = data.rm;
  const Real &r1 = data.r1;
  const Real &r2 = data.r2;

  R1_alpha_p = R1(acos_x3, alpha_p, alpha_p2);
  R1_alpha_m = R1(acos_x3, alpha_m, alpha_m2);
  F3 = boost::math::ellint_1(acos_x3, k3) / mp::sqrt(A * B);
  Ip = -(((A + B) * F3 + (2 * mp::sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                         (A * (r1 - rp) + B * (-r2 + rp))) /
         (-(A * r1) - B * r2 + (A + B) * rp));
  Im = -(((A + B) * F3 + (2 * mp::sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                         (A * (r1 - rm) + B * (-r2 + rm))) /
         (-(A * r1) - B * r2 + (A + B) * rm));
  integral[0] = F3;
  integral[1] = (a * (Im * (-(a * lambda) + 2 * rm) + Ip * (a * lambda - 2 * rp))) / (rm - rp);
}

void IIntegral3::calc(bool is_plus) {
  pre_calc();
  calc_x(integral_rs, acos_x3_rs);
  calc_x(integral_inf, acos_x3_inf);

  auto &radial_integrals = data.radial_integrals;

  for (int i = 0; i < 2; ++i) {
    if (is_plus) {
      radial_integrals[i] = integral_inf[i] + integral_rs[i];
    } else {
      radial_integrals[i] = integral_inf[i] - integral_rs[i];
    }
  }
}

void GIntegral::reinit() {
  const Real &a = data.a;
  const Real &up = data.up;
  const Real &um = data.um;

  up_over_um = up / um;
  one_over_sqrt_up = mp::sqrt(up);
  one_over_umaa_sqrt = 1 / mp::sqrt(-um * a * a);

  G_theta_p[0] = -boost::math::ellint_1(up_over_um) * one_over_umaa_sqrt;
  G_theta_p[1] = -boost::math::ellint_3(up, up_over_um) * one_over_umaa_sqrt;
  G_theta_p[2] = um * (-boost::math::ellint_2(up_over_um) * one_over_umaa_sqrt + G_theta_p[0]);
}

void GIntegral::G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta) {
  const Real &up = data.up;
  const Real &um = data.um;
  asin_up_cos_theta = mp::asin(mp::cos(theta) * one_over_sqrt_up);
  G_arr[0] = -boost::math::ellint_1(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
  G_arr[1] = -boost::math::ellint_3(up, asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
  G_arr[2] = um * (boost::math::ellint_2(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt + G_arr[0]);
}

void GIntegral::calc() {
  const Real &a = data.a;
  const Real &up = data.up;
  const Real &um = data.um;
  auto &tau_o = data.tau_o;
  const Real &theta_s = data.theta_s;
  Sign nu_theta = data.nu_theta;

  G_theta_phi_t(G_theta_s, theta_s);

  const Real &G_theta_theta_s = G_theta_s[0];
  const Real &G_theta_theta_p = G_theta_p[0];

  data.theta_inf = mp::acos(-mp::sqrt(up) * to_integral(nu_theta) *
                            boost::math::jacobi_sn(
                                (tau_o + to_integral(nu_theta) * G_theta_theta_s) / one_over_umaa_sqrt, up_over_um));

  // Angular integrals
  Real m_Real = 1 + mp::floor(mp::real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                       (2 * G_theta_theta_p)));

  using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
      boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;

  // floor
  data.m = RealToInt::convert(m_Real);

  // Number of half-orbits
  data.n_half = tau_o / (2 * G_theta_theta_p);

  G_theta_phi_t(G_theta_inf, data.theta_inf);

  auto &angular_integrals = data.angular_integrals;
  // (-1)^m
  int m1_m = 1 - ((data.m & 1) << 1);

  for (int i = 0; i < 3; ++i) {
    angular_integrals[i] =
        (2 * data.m) * G_theta_p[i] + to_integral(nu_theta) * (m1_m * G_theta_inf[i] - G_theta_s[i]);
  }
}
