#pragma once

#include "Common.h"

// Radial Antiderivatives for case (3)
template<typename Real, typename Complex>
class IIntegral3 {
private:
  ForwardRayTracing<Real, Complex> &data;

  Real acos_x3_rs, acos_x3_inf;
  Real r34_re, r34_im;
  Real A, B, alpha_p, alpha_p2, alpha_m, alpha_m2, k3;
  Real F3, R1_alpha_p, R1_alpha_m, Ip, Im;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_inf;


  Real R1(const Real &acos_x3, const Real &alpha, const Real &alpha2) const {
    return 1 / (1 - alpha2) *
           (boost::math::ellint_3(alpha2 / (alpha2 - 1), acos_x3, k3) -
            alpha * ((sqrt((-1 + alpha2) /
                           (alpha2 + k3 - alpha2 * k3)) *
                      log(abs((sin(acos_x3) +
                               sqrt((-1 + alpha2) /
                                    (alpha2 + k3 - alpha2 * k3)) *
                               sqrt(1 - k3 * square(sin(acos_x3)))) /
                              (-sin(acos_x3) + sqrt((-1 + alpha2) /
                                                    (alpha2 + k3 - alpha2 * k3)) *
                                               sqrt(1 - k3 * square(sin(acos_x3))))))) *
                     half<Real>()));
  }
public:
  explicit IIntegral3(ForwardRayTracing<Real, Complex> &parent) : data(parent) {
  }

  void pre_calc() {
    const Real &rp = data.rp;
    const Real &rm = data.rm;
    const Real &r_s = data.r_s;
    // two real roots, both inside horizon, r_1 < r_2 < r_- < r_+ and r_3 = conj(r_4)
    const Real &r1 = data.r1;
    const Real &r2 = data.r2;
    const Complex &r3 = data.r3_c;
    // const Complex &r4 = data.r4_c;

    r34_re = real(r3);
    r34_im = imag(r3);

    // radial coeffs
    A = sqrt(square(r34_im) + square(r34_re - r2));
    B = sqrt(square(r34_im) + square(r34_re - r1));

    k3 = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
    // alpha_0 = (B + A) / (B - A);
    alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
    alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));
    alpha_p2 = square(alpha_p);
    alpha_m2 = square(alpha_m);

    acos_x3_rs = acos(-1 + (2 * A * (r_s - r1)) / (A * (r_s - r1) + B * (r_s - r2)));
    acos_x3_inf = acos((A - B) / (A + B));
  }

  void calc_x(std::array<Real, 3> &integral, const Real &acos_x3) {
    const Real &a = data.a;
    const Real &lambda = data.lambda;
    const Real &rp = data.rp;
    const Real &rm = data.rm;
    const Real &r1 = data.r1;
    const Real &r2 = data.r2;

    R1_alpha_p = R1(acos_x3, alpha_p, alpha_p2);
    R1_alpha_m = R1(acos_x3, alpha_m, alpha_m2);
    F3 = boost::math::ellint_1(acos_x3, k3) / sqrt(A * B);
    Ip = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                           (A * (r1 - rp) + B * (-r2 + rp))) /
           (-(A * r1) - B * r2 + (A + B) * rp));
    Im = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                           (A * (r1 - rm) + B * (-r2 + rm))) /
           (-(A * r1) - B * r2 + (A + B) * rm));
    integral[0] = F3;
    integral[1] = (a * (Im * (-(a * lambda) + 2 * rm) + Ip * (a * lambda - 2 * rp))) / (rm - rp);
  }

  void calc(bool is_plus) {
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
};
