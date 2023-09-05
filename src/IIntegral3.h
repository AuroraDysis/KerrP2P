#pragma once

#include "Common.h"

#include <boost/math/special_functions/ellint_rc.hpp>

// Radial Antiderivatives for case (3)
template<typename Real, typename Complex>
class IIntegral3 {
#ifdef TESTS
public:
#else
  private:
#endif
  ForwardRayTracing<Real, Complex> &data;

  Real ellint_phi_rs, ellint_phi_ro;
  Real r34_re, r34_im;
  Real A, B, alpha_p, alpha_p2, alpha_m, alpha_m2, ellint_k, ellint_m;
  Real ellint3_n, ellint3_n_new, ellint3_c, ellint_3_tmp, f1;
  Real F3, R1_alpha_p, R1_alpha_m, Ip, Im;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_ro;

  Real R1(const Real &ellint_phi, const Real &alpha, const Real &alpha2) {
    ellint3_n = alpha2 / (alpha2 - 1);
    ellint3_n_new = ellint_m / ellint3_n;
    ellint3_c = 1 / square(sin(ellint_phi));


    // \Pi\left(\phi, \alpha^2, k\right)-F(\phi, k)=\frac{1}{3} \alpha^2 R_J\left(c-1, c-k^2, c, c-\alpha^2\right)
    using boost::math::ellint_rc;
    using boost::math::ellint_rj;
    Real c = 1 / square(sin(ellint_phi));;
    Real k2 = ellint_m;
    Real omega2 = k2 / ellint3_n;

    using boost::math::constants::half_pi;
    using boost::math::constants::two_thirds;

    // https://dlmf.nist.gov/19.25.E14
    // https://dlmf.nist.gov/19.25.E16
    ellint_3_tmp = -third<Real>() * omega2 * ellint_rj(c - 1, c - ellint_m, c, c - omega2)
                   + sqrt((c - 1) * (c - ellint_m) / ((ellint3_n - 1) * (1 - omega2))) *
                     ellint_rc(c * (ellint3_n - 1) * (1 - omega2), (ellint3_n - c) * (c - omega2));
    if (ellint_phi >= half_pi<Real>()) {
      ellint_3_tmp += two_thirds<Real>() * ellint_m / ellint3_n * ellint_rj(0, 1 - ellint_m, 1, 1 - ellint_m / ellint3_n);
    }
//    using boost::math::ellint_rj;
//    ellint_3_tmp = ellint_1(ellint_k, ellint_phi) + third<Real>() * alpha2 * ellint_rj(ellint3_c - 1, ellint3_c - ellint_m,
//                                                                             ellint3_c, ellint3_c - alpha2);

    f1 = half<Real>() * sqrt((-1 + alpha2) / (alpha2 + ellint_m - alpha2 * ellint_m)) *
         log(abs((sin(ellint_phi) + sqrt((-1 + alpha2) /
                                         (alpha2 + ellint_m - alpha2 * ellint_m)) *
                                    sqrt(1 - ellint_m * square(sin(ellint_phi)))) /
                 (-sin(ellint_phi) + sqrt((-1 + alpha2) /
                                          (alpha2 + ellint_m - alpha2 * ellint_m)) *
                                     sqrt(1 - ellint_m * square(sin(ellint_phi))))));
#ifdef PRINT_DEBUG
    fmt::println("R1 - k: {}, n: {}, n1: {}, phi: {}", ellint_k, ellint3_n, ellint3_n_new, ellint_phi);
    fmt::println("R1 - ellint_3: {}, f1: {}", ellint_3_tmp, f1);
#endif
    // sign is different from Mathematica
    return (ellint_3_tmp + alpha * f1) / (alpha2 - 1);
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
    const Real &r_o = data.r_o;

    r34_re = real(r3);
    r34_im = imag(r3);

    // radial coeffs
    A = sqrt(square(r34_im) + square(r34_re - r2));
    B = sqrt(square(r34_im) + square(r34_re - r1));

    ellint_m = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
    ellint_k = sqrt(ellint_m);
    // alpha_0 = (B + A) / (B - A);
    alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
    alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));
    alpha_p2 = square(alpha_p);
    alpha_m2 = square(alpha_m);

    ellint_phi_rs = acos(-1 + (2 * A * (r_s - r1)) / (A * (r_s - r1) + B * (r_s - r2)));
    ellint_phi_ro = acos(-1 + (2 * A * (r_o - r1)) / (A * (r_o - r1) + B * (r_o - r2)));

#ifdef PRINT_DEBUG
    fmt::println("I3 - A: {}, B: {}, ellint_k: {}", A, B, ellint_k);
    fmt::println("I3 - alpha_p: {}, alpha_m: {}, phi_rs: {}, phi_ro: {}", alpha_p, alpha_m, ellint_phi_rs,
                 ellint_phi_ro);
#endif
  }

  void calc_x(std::array<Real, 3> &integral, const Real &ellint_phi) {
    const Real &a = data.a;
    const Real &lambda = data.lambda;
    const Real &rp = data.rp;
    const Real &rm = data.rm;
    const Real &r1 = data.r1;
    const Real &r2 = data.r2;

    R1_alpha_p = R1(ellint_phi, alpha_p, alpha_p2);
    R1_alpha_m = R1(ellint_phi, alpha_m, alpha_m2);
    F3 = ellint_1(ellint_k, ellint_phi) / sqrt(A * B);
    Ip = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                           (A * (r1 - rp) + B * (-r2 + rp))) /
           (-(A * r1) - B * r2 + (A + B) * rp));
    Im = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                           (A * (r1 - rm) + B * (-r2 + rm))) /
           (-(A * r1) - B * r2 + (A + B) * rm));
    integral[0] = F3;
    integral[1] = (a * (Im * (-(a * lambda) + 2 * rm) + Ip * (a * lambda - 2 * rp))) / (rm - rp);

#ifdef PRINT_DEBUG
    fmt::println("I3 - F3: {}, R1_alpha_p: {}, R1_alpha_m: {}, Ip: {}, Im: {}", F3, R1_alpha_p, R1_alpha_m, Ip, Im);
    fmt::println("I3 - I_r: {}, I_phi: {}, T_t: {}", integral[0], integral[1], integral[2]);
#endif
  }

  void calc(bool is_plus) {
    pre_calc();
    calc_x(integral_rs, ellint_phi_rs);
    calc_x(integral_ro, ellint_phi_ro);

    auto &radial_integrals = data.radial_integrals;

    for (int i = 0; i < 2; ++i) {
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
