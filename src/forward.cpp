#include "forward.h"


//  Real x3(Real r) const {
//    return (A * (r - r1) - B * (r - r2)) / (A * (r - r1) + B * (r - r2));
//  }
//
//  Real F3(Real r) const {
//    return 1 / mp::sqrt(A * B) * boost::math::ellint_1(mp::acos(x3(r)), k3);
//  }
//
//  Real f1(Real alpha, Real curlyPhi, Real j) const {
//    Real temp1 = mp::sqrt((alpha * alpha - 1) / (j + (1 - j) * alpha * alpha));
//    Real temp2 = temp1 * mp::sqrt(1 - j * mp::sin(curlyPhi) * mp::sin(curlyPhi));
//    return temp1 / 2 * mp::log(mp::abs((temp2 + mp::sin(curlyPhi)) / (temp2 - mp::sin(curlyPhi))));
//  }
//
//  Real R1(Real alpha, Real curlyPhi, Real j) const {
//    return 1 / (1 - alpha * alpha) *
//           (boost::math::ellint_3(alpha * alpha / (alpha * alpha - 1), curlyPhi, j) -
//            alpha * f1(alpha, curlyPhi, j));
//  }
//
//  Real R2(Real alpha, Real curlyPhi, Real j) const {
//    return 1 / (alpha * alpha - 1) * (boost::math::ellint_1(curlyPhi, j) -
//                                      mp::pow(alpha, 2) / (j + (1 - j) * alpha * alpha) *
//                                      (boost::math::ellint_2(curlyPhi, j) -
//                                       (alpha * mp::sin(curlyPhi) *
//                                        mp::sqrt(1 - j * mp::pow(mp::sin(curlyPhi), 2))) /
//                                       (1 + alpha * mp::cos(curlyPhi)))) +
//           1 / (j + (1 - j) * alpha * alpha) * (2 * j - mp::pow(alpha, 2) / (alpha * alpha - 1)) *
//           R1(alpha, curlyPhi, j);
//  }
//
//  Real Pi13(Real r) const {
//    return (2 * (r2 - r1) * mp::sqrt(A * B)) / (B * B - A * A) * R1(alpha_0, mp::acos(x3(r)), k3);
//  }
//
//  Real Pi23(Real r) const {
//    return mp::pow(((2 * (r2 - r1) * mp::sqrt(A * B)) / (B * B - A * A)), 2) *
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
//    return mp::pow(((B * r2 + A * r1) / (B + A)), 2) * F3(r) +
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

//void calc_I() {
//  // radial coeffs
//  Real A, B, alpha_0, alpha_p, alpha_m, k3;
//  A = mp::sqrt((r3 - r2) * (r4 - r2));
//  B = mp::sqrt((r3 - r1) * (r4 - r1));
//
//  alpha_0 = (B + A) / (B - A);
//  alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
//  alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));
//
//  k3 = (mp::pow(A + B, 2) - mp::pow(r2 - r1, 2)) / (4 * A * B);
//
//  std::array<std::function<Real(Real)>, 3> anti_ders = {
//      [&](Real r) -> Real { return I_r(r); },
//      [&](Real r) -> Real { return I_phi(r); },
//      [&](Real r) -> Real { return I_t(r); }
//  };
//
//  std::array<Real, 3> results = {};
//  std::fill(results.begin(), results.end(), std::numeric_limits<Real>::quiet_NaN());
//  bool radial_turning = r34_is_real && r4 > rp;
//
//  // if there is a radial turning point (i.e. r4 is a real number)
//  if (radial_turning && r_s <= r4) {
//    ray_status = RayStatus::CONFINED;
//    return;
//  }
//
//  if (radial_turning && r_s > r4 && nu_r == Sign::POSITIVE) {
//    for (int i = 0; i < 3; i++) {
//      results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) - anti_ders[i](r_s);
//    }
//    return;
//  }
//
//  if (radial_turning && r_s > r4 && nu_r == Sign::NEGATIVE) {
//    for (int i = 0; i < 3; i++) {
//      results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) + anti_ders[i](r_s) - 2 * anti_ders[i](r4);
//    }
//    return;
//  }
//
//  if (!radial_turning && nu_r == Sign::NEGATIVE) {
//    ray_status = RayStatus::FALLS_IN;
//    return;
//  }
//
//  if (!radial_turning && nu_r == Sign::POSITIVE) {
//    for (int i = 0; i < 3; i++) {
//      results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) - anti_ders[i](r_s);
//    }
//    return;
//  }
//
//  ray_status = RayStatus::UNKOWN_ERROR;
//}


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
