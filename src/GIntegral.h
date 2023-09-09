#pragma once

#include "Common.h"
#include "Integral.h"

#include <boost/numeric/conversion/converter.hpp>

template<typename Real, typename Complex>
class GIntegral : public Integral<Real, Complex> {
#ifdef TESTS
  public:
#else
private:
#endif
  std::array<Real, 3> G_theta_p = {};
  std::array<Real, 3> G_theta_s = {};
  std::array<Real, 3> G_theta_f = {};

  Real one_over_sqrt_up;
  Real ellint_m, ellint_k;
  Real ellint_kappa, ellint_kappa2, ellint_kappa_prime, ellint_one_over_kappa_prime, ellint3_n, ellint_alpha1_2;
  Real jacobi_sn_k1, jacobi_sn_k1_prime;
  Real one_over_umaa_sqrt;

  Real ellint_sin_phi, ellint_cos_theta, ellint_sin_theta, ellint_cos_theta2, ellint_sin_theta2, ellint_y;
  Real ellint_1_phi, ellint_2_phi, ellint_3_phi;

  void G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta) {
    const Real &up = this->data.up;
    const Real &um = this->data.um;

    // https://dlmf.nist.gov/19.7.E5
    ellint_sin_phi = cos(theta) * one_over_sqrt_up;
    CHECK_VAR_INT_RANGE(ellint_sin_phi, -1, 1);
    ellint_sin_theta = (sqrt(1 + ellint_m) * ellint_sin_phi) / sqrt(1 + ellint_m * MY_SQUARE(ellint_sin_phi));
    ellint_sin_theta2 = MY_SQUARE(ellint_sin_theta);
    CHECK_VAR_INT_RANGE(ellint_sin_theta, -1, 1);
    // arcsin gives -pi/2 to pi/2, so cos(theta) is always positive
    ellint_cos_theta2 = 1 - ellint_sin_theta2;
    ellint_cos_theta = sqrt(ellint_cos_theta2);
    CHECK_VAR_INT_RANGE(ellint_cos_theta, 0, 1);

    ellint_y = 1 - ellint_kappa2 * ellint_sin_theta2;
    // ellint_1_phi = boost::math::ellint_1(ellint_kappa, ellint_theta);
    ellint_1_phi = ellint_sin_theta * boost::math::ellint_rf(ellint_cos_theta2, ellint_y, 1);
    // ellint_2_phi = boost::math::ellint_2(ellint_kappa, ellint_theta);
    ellint_2_phi = ellint_1_phi - third<Real>() * ellint_kappa2 * ellint_sin_theta2 * ellint_sin_theta *
                                  boost::math::ellint_rd(ellint_cos_theta2, ellint_y, 1);
    // ellint_3_phi = boost::math::ellint_3(ellint_kappa, ellint_alpha1_2, ellint_theta);
    ellint_3_phi = ellint_1_phi + third<Real>() * ellint_alpha1_2 * ellint_sin_theta2 * ellint_sin_theta *
                                  boost::math::ellint_rj(ellint_cos_theta2, ellint_y, 1,
                                                         1 - ellint_alpha1_2 * ellint_sin_theta2);
    G_arr[0] = -one_over_umaa_sqrt * ellint_kappa_prime * ellint_1_phi;
    G_arr[1] = -one_over_umaa_sqrt * ellint_kappa_prime / ellint_alpha1_2 *
               (MY_SQUARE(ellint_kappa_prime) * up *
                ellint_3_phi +
                ellint_kappa2 * ellint_1_phi);
    G_arr[2] = um *
               (one_over_umaa_sqrt * ellint_one_over_kappa_prime * (ellint_2_phi -
                                                                    ellint_kappa2 * ellint_sin_theta *
                                                                    ellint_cos_theta / sqrt(ellint_y)) +
                G_arr[0]);
  }

public:
  explicit GIntegral(ForwardRayTracing<Real, Complex> &data_) : Integral<Real, Complex>(data_, typeid(*this).name()) {
  }

  void calc() {
    const Real &a = this->data.a;
    const Real &up = this->data.up;
    const Real &um = this->data.um;
    Real &tau_o = this->data.tau_o;
    const Real &theta_s = this->data.theta_s;
    Sign nu_theta = this->data.nu_theta;

    // up / um < 0, Imaginary-Modulus Transformation
    // https://dlmf.nist.gov/19.7#ii
    // https://analyticphysics.com/Special%20Functions/A%20Miscellany%20of%20Elliptic%20Integrals.htm
    ellint_m = -up / um;
    ellint_k = sqrt(ellint_m);
    ellint_kappa2 = ellint_m / (ellint_m + 1);
    ellint_kappa = sqrt(ellint_kappa2);
    ellint_one_over_kappa_prime = sqrt(ellint_m + 1);
    ellint_kappa_prime = 1 / ellint_one_over_kappa_prime;
    ellint3_n = up / (up - 1);
    ellint_alpha1_2 = (ellint_m + up) / (1 + ellint_m);

    one_over_sqrt_up = 1 / sqrt(up);
    one_over_umaa_sqrt = 1 / sqrt(-um * a * a);

    G_theta_p[0] = ellint_kappa_prime * boost::math::ellint_1(ellint_kappa) * one_over_umaa_sqrt;
    G_theta_p[1] =
        (ellint_kappa_prime / (1 - up)) * boost::math::ellint_3(ellint_kappa, ellint3_n) * one_over_umaa_sqrt;
    if (this->data.calc_t_f) {
      G_theta_p[2] =
          um * (-ellint_one_over_kappa_prime * boost::math::ellint_2(ellint_kappa) * one_over_umaa_sqrt + G_theta_p[0]);
    } else {
      G_theta_p[2] = std::numeric_limits<Real>::quiet_NaN();
    }

    G_theta_phi_t(G_theta_s, theta_s);
    if (this->data.ray_status != RayStatus::NORMAL) {
      return;
    }

    const Real &G_theta_theta_s = G_theta_s[0];
    const Real &G_theta_theta_p = G_theta_p[0];

    // https://dlmf.nist.gov/22.17
    jacobi_sn_k1_prime = 1 / sqrt(1 + ellint_m);
    jacobi_sn_k1 = ellint_k * jacobi_sn_k1_prime;
    this->data.theta_f = acos(-sqrt(up) * to_integral(nu_theta) *
                              jacobi_sn_k1_prime *
                              boost::math::jacobi_sd(jacobi_sn_k1, (tau_o + to_integral(nu_theta) * G_theta_theta_s) /
                                                                   (one_over_umaa_sqrt * jacobi_sn_k1_prime)));

    // Angular integrals
    Real m_Real = 1 + floor(real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                 (2 * G_theta_theta_p)));

    // floor
    this->data.m = MY_FLOOR<Real>::convert(m_Real);

    // Number of half-orbits
    this->data.n_half = tau_o / (2 * G_theta_theta_p);

    G_theta_phi_t(G_theta_f, this->data.theta_f);
    if (this->data.ray_status != RayStatus::NORMAL) {
      return;
    }

    auto &angular_integrals = this->data.angular_integrals;
    // (-1)^m
    int m1_m = 1 - ((this->data.m & 1) << 1);

    int ix = this->data.calc_t_f ? 3 : 2;
    for (int i = 0; i < ix; ++i) {
      angular_integrals[i] =
          (2 * this->data.m) * G_theta_p[i] + to_integral(nu_theta) * (m1_m * G_theta_f[i] - G_theta_s[i]);
    }
    if (!this->data.calc_t_f) {
      angular_integrals[2] = std::numeric_limits<Real>::quiet_NaN();
    }

#ifdef PRINT_DEBUG
    fmt::println("ellint_kappa: {}", ellint_kappa);
    fmt::println("G_theta_p: {}, {}, {}", G_theta_p[0], G_theta_p[1], G_theta_p[2]);
    fmt::println("G_theta_s: {}, {}, {}", G_theta_s[0], G_theta_s[1], G_theta_s[2]);
    fmt::println("G_theta_f: {}, {}, {}", G_theta_f[0], G_theta_f[1], G_theta_f[2]);
    fmt::println("G: {}, {}, {}", angular_integrals[0], angular_integrals[1], angular_integrals[2]);
#endif
  }
};
