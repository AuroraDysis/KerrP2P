#pragma once

#include "Common.h"

template<typename Real, typename Complex>
class GIntegral {
#ifdef TESTS
public:
#else
  private:
#endif
  std::array<Real, 3> G_theta_p = {};
  std::array<Real, 3> G_theta_s = {};
  std::array<Real, 3> G_theta_inf = {};

  Real one_over_sqrt_up;
  Real ellint_m;
  Real ellint_k;
  Real ellint1_coeff, ellint2_coeff, ellint3_coeff, ellint3_n;
  Real one_over_umaa_sqrt;
  Real ellint_phi; // ArcCsc[sqrt[up] Sec[\[Theta]]]
  Real ellint_theta;

  ForwardRayTracing<Real, Complex> &data;

  void G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta) {
    const Real &up = data.up;
    const Real &um = data.um;
    ellint_phi = asin(cos(theta) * one_over_sqrt_up);
    ellint_theta = asin(
        (sqrt(1 + square(ellint_k)) * sin(ellint_phi)) / sqrt(1 + square(ellint_k) * square(sin(ellint_phi))));
    G_arr[0] = -one_over_umaa_sqrt * ellint1_coeff * boost::math::ellint_1(ellint_k, ellint_theta);
    G_arr[1] = -one_over_umaa_sqrt * ellint3_coeff * boost::math::ellint_3(ellint_k, ellint3_n, ellint_theta);
    G_arr[2] = um * (-one_over_umaa_sqrt * ellint2_coeff * boost::math::ellint_2(ellint_k, ellint_theta) + G_arr[0]);
  }

public:
  explicit GIntegral(ForwardRayTracing<Real, Complex> &parent) : data(parent) {
  }

  void calc() {
    const Real &a = data.a;
    const Real &up = data.up;
    const Real &um = data.um;
    auto &tau_o = data.tau_o;
    const Real &theta_s = data.theta_s;
    Sign nu_theta = data.nu_theta;

    // https://dlmf.nist.gov/19.7#ii
    // https://analyticphysics.com/Special%20Functions/A%20Miscellany%20of%20Elliptic%20Integrals.htm
    ellint_m = -up / um;
    ellint_k = sqrt(ellint_m / (ellint_m + 1));
    ellint1_coeff = 1 / sqrt(ellint_m + 1);
    ellint2_coeff = sqrt(ellint_m + 1);
    ellint3_coeff = ellint1_coeff / (1 - up);
    ellint3_n = up / (up - 1);

    one_over_sqrt_up = sqrt(up);
    one_over_umaa_sqrt = 1 / sqrt(-um * a * a);

    G_theta_p[0] = ellint1_coeff * boost::math::ellint_1(ellint_k) * one_over_umaa_sqrt;
    G_theta_p[1] = ellint3_coeff * boost::math::ellint_3(ellint_k, ellint3_n) * one_over_umaa_sqrt;
    if (data.calc_t_f) {
      G_theta_p[2] = um * (-ellint2_coeff * boost::math::ellint_2(ellint_k) * one_over_umaa_sqrt + G_theta_p[0]);
    } else {
      G_theta_p[2] = std::numeric_limits<Real>::quiet_NaN();
    }

    G_theta_phi_t(G_theta_s, theta_s);

    const Real &G_theta_theta_s = G_theta_s[0];
    const Real &G_theta_theta_p = G_theta_p[0];

    data.theta_f = acos(-sqrt(up) * to_integral(nu_theta) *
                        boost::math::jacobi_sn(
                            (tau_o + to_integral(nu_theta) * G_theta_theta_s) / one_over_umaa_sqrt, ellint_k));

    // Angular integrals
    Real m_Real = 1 + floor(real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                 (2 * G_theta_theta_p)));

    using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
        boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;

    // floor
    data.m = RealToInt::convert(m_Real);

    // Number of half-orbits
    data.n_half = tau_o / (2 * G_theta_theta_p);

    G_theta_phi_t(G_theta_inf, data.theta_f);

    auto &angular_integrals = data.angular_integrals;
    // (-1)^m
    int m1_m = 1 - ((data.m & 1) << 1);

    int ix = data.calc_t_f ? 3 : 2;
    for (int i = 0; i < ix; ++i) {
      angular_integrals[i] =
          (2 * data.m) * G_theta_p[i] + to_integral(nu_theta) * (m1_m * G_theta_inf[i] - G_theta_s[i]);
    }
    if (!data.calc_t_f) {
      angular_integrals[2] = std::numeric_limits<Real>::quiet_NaN();
    }

#ifdef PRINT_DEBUG
    fmt::println("ellint_k: {}", ellint_k);
    fmt::println("G_theta_p: {}, {}, {}", G_theta_p[0], G_theta_p[1], G_theta_p[2]);
    fmt::println("G_theta_s: {}, {}, {}", G_theta_s[0], G_theta_s[1], G_theta_s[2]);
    fmt::println("G: {}, {}, {}", angular_integrals[0], angular_integrals[1], angular_integrals[2]);
#endif
  }
};
