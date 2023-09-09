#pragma once

#include "Common.h"
#include "Integral.h"

#include <boost/math/special_functions/ellint_rc.hpp>

// Radial Antiderivatives for case (3)
template<typename Real, typename Complex>
class IIntegral3 : public Integral<Real, Complex> {
#ifdef TESTS
public:
#else
  private:
#endif
  Real ellint_phi, ellint_cos_phi, ellint_sin_phi;
  Real r34_re, r34_im;
  Real A, B, alpha_p, alpha_m, ellint_k, ellint_m, alpha2, ellint1_phi;
  Real alpha_0, R1_alpha_0, R2_alpha_0, Pi_13, Pi_23, I_1, I_2;
  Real ellint3_n, ellint_c, ellint_3_tmp, f1, ellint3_n1;
  Real F3, R1_alpha_p, R1_alpha_m, I_p, I_m;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_ro;

  Real R1(const Real &alpha) {
    alpha2 = MY_SQUARE(alpha);
    ellint3_n = alpha2 / (alpha2 - 1);
    ellint3_n1 = ellint_m / ellint3_n;

    // \Pi\left(\phi, \alpha^2, k\right)-F(\phi, k)=\frac{1}{3} \alpha^2 R_J\left(c-1, c-k^2, c, c-\alpha^2\right)
    using boost::math::ellint_rc;
    using boost::math::ellint_rj;

    // https://dlmf.nist.gov/19.25.E14
    // https://dlmf.nist.gov/19.25.E16
    using boost::math::constants::half_pi;
    using boost::math::constants::two_thirds;
    ellint_3_tmp = -third<Real>() * ellint3_n1 * ellint_rj(ellint_c - 1, ellint_c - ellint_m, ellint_c, ellint_c - ellint3_n1)
                   + sqrt((ellint_c - 1) * (ellint_c - ellint_m) / ((ellint3_n - 1) * (1 - ellint3_n1))) *
                     ellint_rc(ellint_c * (ellint3_n - 1) * (1 - ellint3_n1), (ellint3_n - ellint_c) * (ellint_c - ellint3_n1));
    if (ellint_phi >= half_pi<Real>()) {
      ellint_3_tmp +=
          two_thirds<Real>() * ellint_m / ellint3_n * ellint_rj(0, 1 - ellint_m, 1, 1 - ellint_m / ellint3_n);
    }

    f1 = half<Real>() * sqrt((-1 + alpha2) / (alpha2 + ellint_m - alpha2 * ellint_m)) *
         log(abs((ellint_sin_phi + sqrt((-1 + alpha2) /
                                         (alpha2 + ellint_m - alpha2 * ellint_m)) *
                                    sqrt(1 - ellint_m * MY_SQUARE(ellint_sin_phi))) /
                 (-ellint_sin_phi + sqrt((-1 + alpha2) /
                                          (alpha2 + ellint_m - alpha2 * ellint_m)) *
                                     sqrt(1 - ellint_m * MY_SQUARE(ellint_sin_phi)))));
#ifdef PRINT_DEBUG
    fmt::println("R1 - k: {}, n: {}, phi: {}", ellint_k, ellint3_n, ellint_phi);
    fmt::println("R1 - ellint_3: {}, f1: {}", ellint_3_tmp, f1);
#endif
    // sign is different from Mathematica
    return (ellint_3_tmp + alpha * f1) / (alpha2 - 1);
  }

public:
  explicit IIntegral3(ForwardRayTracing<Real, Complex> &data_) : Integral<Real, Complex>(data_, typeid(*this).name()) {
  }

  void pre_calc() {
    const Real &rp = this->data.rp;
    const Real &rm = this->data.rm;
    const Real &r_s = this->data.r_s;
    // two real roots, both inside horizon, r_1 < r_2 < r_- < r_+ and r_3 = conj(r_4)
    const Real &r1 = this->data.r1;
    const Real &r2 = this->data.r2;
    const Complex &r3 = this->data.r3_c;
    const Real &r_o = this->data.r_o;

    r34_re = real(r3);
    r34_im = imag(r3);

    // radial coeffs
    A = sqrt(MY_SQUARE(r34_im) + MY_SQUARE(r34_re - r2));
    B = sqrt(MY_SQUARE(r34_im) + MY_SQUARE(r34_re - r1));

    ellint_m = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
    ellint_k = sqrt(ellint_m);

    alpha_0 = (B + A) / (B - A);
    alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
    alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));

#ifdef PRINT_DEBUG
    fmt::println("I3 - A: {}, B: {}, ellint_k: {}", A, B, ellint_k);
    fmt::println("I3 - alpha_p: {}, alpha_m: {}", alpha_p, alpha_m);
#endif
  }

  void calc_x(std::array<Real, 3> &integral, const Real &r) {
    const Real &a = this->data.a;
    const Real &lambda = this->data.lambda;
    const Real &rp = this->data.rp;
    const Real &rm = this->data.rm;
    const Real &r1 = this->data.r1;
    const Real &r2 = this->data.r2;
    ellint_cos_phi = -1 + (2 * A * (r - r1)) / (A * (r - r1) + B * (r - r2));
    CHECK_VAR_INT_RANGE(ellint_cos_phi, -1, 1);
    ellint_sin_phi = sqrt(1 - MY_SQUARE(ellint_cos_phi));
    CHECK_VAR_INT_RANGE(ellint_sin_phi, -1, 1);
    ellint_c = 1 / (1 - MY_SQUARE(ellint_cos_phi));
    ellint_phi = acos(ellint_cos_phi);

    R1_alpha_p = R1(alpha_p);
    R1_alpha_m = R1(alpha_m);
    ellint1_phi = ellint_1(ellint_k, ellint_phi);
    F3 = ellint1_phi / sqrt(A * B);
    I_p = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                            (A * (r1 - rp) + B * (-r2 + rp))) /
            (-(A * r1) - B * r2 + (A + B) * rp));
    I_m = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                            (A * (r1 - rm) + B * (-r2 + rm))) /
            (-(A * r1) - B * r2 + (A + B) * rm));
    integral[0] = F3;
    integral[1] = (a * (I_m * (-(a * lambda) + 2 * rm) + I_p * (a * lambda - 2 * rp))) / (rm - rp);
    if (this->data.calc_t_f && !isinf(this->data.r_o)) {
      alpha2 = MY_SQUARE(alpha_0);
      R1_alpha_0 = R1(alpha_0);
      R2_alpha_0 = ((-1 + alpha2) * ellint_m * (1 + alpha_0 * ellint_cos_phi) *
                    (ellint1_phi - 2 * R1_alpha_0) +
                    alpha2 * (1 + alpha_0 * ellint_cos_phi) *
                    (ellint_2(ellint_k, ellint_phi) - ellint1_phi + R1_alpha_0) -
                    MY_CUBE(alpha_0) * ellint_sin_phi * sqrt(1 - ellint_m * MY_SQUARE(ellint_sin_phi))) /
                   ((-1 + alpha2) * (alpha2 * (-1 + ellint_m) - ellint_m) *
                    (1 + alpha_0 * ellint_cos_phi));
      Pi_13 = (2 * sqrt(A * B) * R1_alpha_0 * (r1 - r2)) / ((A - B) * (A + B));
      Pi_23 = (4 * (A * B) * R2_alpha_0 * MY_SQUARE(r1 - r2)) / MY_SQUARE((A - B) * (A + B));
      I_1 = (A * r1 + B * r2) * F3 / (A + B) + Pi_13;
      I_2 =
          sqrt(A * B) * Pi_23 + ((A * r1 + B * r2) * (2 * (A + B) * Pi_13 + F3 * (A * r1 + B * r2))) / MY_SQUARE(A + B);
      integral[2] =
          4 * F3 + I_2 + (2 * I_m * rm * (-(a * lambda) + 2 * rm) + 2 * I_p * (a * lambda - 2 * rp) * rp) / (rm - rp) +
          2 * I_1;

#ifdef PRINT_DEBUG
      fmt::println("I3 - R1_alpha_0: {}, R2_alpha_0: {}, Pi_13: {}, Pi_23: {}, I_1: {}, I_2: {}", R1_alpha_0,
                   R2_alpha_0, Pi_13, Pi_23, I_1, I_2);
#endif
    } else {
      integral[2] = std::numeric_limits<Real>::quiet_NaN();
    }

#ifdef PRINT_DEBUG
    fmt::println("I3 - F3: {}, R1_alpha_p: {}, R1_alpha_m: {}, I_p: {}, I_m: {}", F3, R1_alpha_p, R1_alpha_m, I_p, I_m);
    fmt::println("I3 - I_r: {}, I_phi: {}, I_t: {}", integral[0], integral[1], integral[2]);
#endif
  }

  void calc(bool is_plus) {
    pre_calc();

    auto &r_s = this->data.r_s;
    auto &r_o = this->data.r_o;

    calc_x(integral_rs, r_s);
    calc_x(integral_ro, r_o);

    auto &radial_integrals = this->data.radial_integrals;

    for (int i = 0; i < 3; ++i) {
      if (is_plus) {
        radial_integrals[i] = integral_ro[i] + integral_rs[i];
      } else {
        radial_integrals[i] = integral_ro[i] - integral_rs[i];
      }
    }

#ifdef PRINT_DEBUG
    fmt::println("I3: {}, {}, {}", radial_integrals[0], radial_integrals[1], radial_integrals[2]);
#endif
  }
};
