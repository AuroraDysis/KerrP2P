#pragma once

#include "Common.h"

template<typename Real, typename Complex>
class GIntegral {
private:
  std::array<Real, 3> G_theta_p = {};
  std::array<Real, 3> G_theta_s = {};
  std::array<Real, 3> G_theta_inf = {};

  Real one_over_sqrt_up;
  Real up_over_um;
  Real one_over_umaa_sqrt;
  Real asin_up_cos_theta; // ArcCsc[Sqrt[up] Sec[\[Theta]]]

  ForwardRayTracing<Real, Complex> &data;

  void G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta) {
    const Real &up = data.up;
    const Real &um = data.um;
    asin_up_cos_theta = asin(cos(theta) * one_over_sqrt_up);
    G_arr[0] = -boost::math::ellint_1(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
    G_arr[1] = -boost::math::ellint_3(up, asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
    G_arr[2] = um * (boost::math::ellint_2(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt + G_arr[0]);
  }

public:
  explicit GIntegral(ForwardRayTracing<Real, Complex> &parent) : data(parent) {
  }

  void pre_calc() {
    const Real &a = data.a;
    const Real &up = data.up;
    const Real &um = data.um;

    up_over_um = up / um;
    one_over_sqrt_up = sqrt(up);
    one_over_umaa_sqrt = 1 / sqrt(-um * a * a);

    G_theta_p[0] = -boost::math::ellint_1(up_over_um) * one_over_umaa_sqrt;
    G_theta_p[1] = -boost::math::ellint_3(up, up_over_um) * one_over_umaa_sqrt;
    if (data.calc_t_f) {
      G_theta_p[2] = um * (-boost::math::ellint_2(up_over_um) * one_over_umaa_sqrt + G_theta_p[0]);
    } else {
      G_theta_p[2] = std::numeric_limits<Real>::quiet_NaN();
    }
  }

  void calc() {
    const Real &a = data.a;
    const Real &up = data.up;
    const Real &um = data.um;
    auto &tau_o = data.tau_o;
    const Real &theta_s = data.theta_s;
    Sign nu_theta = data.nu_theta;

    G_theta_phi_t(G_theta_s, theta_s);

    const Real &G_theta_theta_s = G_theta_s[0];
    const Real &G_theta_theta_p = G_theta_p[0];

    data.theta_f = acos(-sqrt(up) * to_integral(nu_theta) *
                        boost::math::jacobi_sn(
                            (tau_o + to_integral(nu_theta) * G_theta_theta_s) / one_over_umaa_sqrt, up_over_um));

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

    for (int i = 0; i < 3; ++i) {
      angular_integrals[i] =
          (2 * data.m) * G_theta_p[i] + to_integral(nu_theta) * (m1_m * G_theta_inf[i] - G_theta_s[i]);
    }
  }
};
