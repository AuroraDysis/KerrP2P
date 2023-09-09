#pragma once

#include "Common.h"
#include "Integral.h"

// Radial Antiderivatives for case (2)
template<typename Real, typename Complex>
class IIntegral2 : public Integral<Real, Complex> {
#ifdef TESTS
public:
#else
private:
#endif
  Real ellint_phi_rs, ellint_phi_ro, ellint_k;
  Real E2_coeff, F2_coeff, Pi_p2_coeff, Pi_m2_coeff, Pi_p2_ellint_n, Pi_m2_ellint_n;

  Real F2, Pi_p2, Pi_m2, I_p, I_m;
  Real E2, Pi_12, I1, I2, Pi_12_ellint_n;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_ro;
public:
  explicit IIntegral2(ForwardRayTracing<Real, Complex> &data_) : Integral<Real, Complex>(data_) {
  }

  void pre_calc() {
    const Real &rp = this->data.rp;
    const Real &rm = this->data.rm;

    const Real &r1 = this->data.r1;
    const Real &r2 = this->data.r2;
    const Real &r3 = this->data.r3;
    const Real &r4 = this->data.r4;
    const Real &r_s = this->data.r_s;
    const Real &r_o = this->data.r_o;

    ellint_phi_rs = asin(sqrt(((r1 - r3) * (r_s - r4)) / ((r_s - r3) * (r1 - r4))));
    if (isinf(r_o)) {
      ellint_phi_ro = asin(sqrt((r1 - r3) / (r1 - r4)));
    } else {
      ellint_phi_ro = asin(sqrt(((r1 - r3) * (r_o - r4)) / ((r_o - r3) * (r1 - r4))));
    }
    ellint_k = sqrt(((-r2 + r3) * (-r1 + r4)) / ((r1 - r3) * (r2 - r4)));

    E2_coeff = sqrt((r1 - r3) * (r2 - r4));
    F2_coeff = 2 / E2_coeff;
    Pi_p2_coeff = 2 * (-r3 + r4) / (E2_coeff * (-r3 + rp) * (-r4 + rp));
    Pi_m2_coeff = 2 * (-r3 + r4) / (E2_coeff * (-r3 + rm) * (-r4 + rm));
    Pi_p2_ellint_n = ((r1 - r4) * (r3 - rp)) / ((r1 - r3) * (r4 - rp));
    Pi_m2_ellint_n = ((r1 - r4) * (r3 - rm)) / ((r1 - r3) * (r4 - rm));
    Pi_12_ellint_n = (r1 - r4) / (r1 - r3);
  }

  void calc_x(std::array<Real, 3>& integral, const Real &ellint_phi, const Real &r) {
    const Real &a = this->data.a;
    const Real &lambda = this->data.lambda;
    const Real &rp = this->data.rp;
    const Real &rm = this->data.rm;
    const Real &r3 = this->data.r3;

    // ellint_k, phi
    F2 = F2_coeff * ellint_1(ellint_k, ellint_phi);
    Pi_p2 = Pi_p2_coeff * ellint_3(ellint_k, Pi_p2_ellint_n, ellint_phi);
    Pi_m2 = Pi_m2_coeff * ellint_3(ellint_k, Pi_m2_ellint_n, ellint_phi);
    I_p = F2 / (r3 - rp) - Pi_p2;
    I_m = F2 / (r3 - rm) - Pi_m2;

#ifdef PRINT_DEBUG
    fmt::println("I2 - ellint_phi: {}, ellint_k: {}", ellint_phi, ellint_k);
    fmt::println("I2 - F2: {}, E2: {}, Pi_p2: {}, Pi_m2: {}, I_p: {}, I_m: {}", F2, E2, Pi_p2, Pi_m2, I_p, I_m);
#endif

    // I_r
    integral[0] = F2;
    // I_phi
    integral[1] = (a * (-2 * rp * I_p + a * I_p * lambda + (2 * rm - a * lambda) * I_m)) / (rm - rp);
    if (this->data.calc_t_f && !isinf(this->data.r_o)) {
      // I_t
      const Real &r4 = this->data.r4;
      const Real &r1 = this->data.r1;
      const Real &r2 = this->data.r2;
      const Real &eta = this->data.eta;

      E2 = E2_coeff * boost::math::ellint_2(ellint_k, ellint_phi);
      Pi_12 = F2_coeff * boost::math::ellint_3(ellint_k, Pi_12_ellint_n, ellint_phi);
      I1 = r3 * F2 + (r4 - r3) * Pi_12;
      I2 = -E2 + sqrt(-((eta + MY_SQUARE(a - lambda)) * (MY_SQUARE(a) + (-2 + r) * r)) +
                      MY_SQUARE(MY_SQUARE(a) - a * lambda + MY_SQUARE(r))) / (r - r3) -
           (F2 * (r2 * r3 + r1 * r4)) * half<Real>();
      integral[2] = 4 * F2 + 2 * I1 + I2 +
                    (-2 * a * I_m * lambda * rm + 4 * I_m * MY_SQUARE(rm) + 2 * I_p * (a * lambda - 2 * rp) * rp) /
                    (rm - rp);
#ifdef PRINT_DEBUG
      fmt::println("I2 - E2: {}, Pi_12: {}, I1: {}, I2: {}", E2, Pi_12, I1, I2);
#endif
    } else {
      integral[2] = std::numeric_limits<Real>::quiet_NaN();
    }
  }

  void calc(bool is_plus) {
    pre_calc();
    const Real &r_s = this->data.r_s;
    const Real &r_o = this->data.r_o;
    calc_x(integral_rs, ellint_phi_rs, r_s);
    calc_x(integral_ro, ellint_phi_ro, r_o);

    auto &radial_integrals = this->data.radial_integrals;

    if (is_plus) {
      for (int i = 0; i < 3; ++i) {
        radial_integrals[i] = integral_ro[i] + integral_rs[i];
      }
    } else {
      for (int i = 0; i < 3; ++i) {
        radial_integrals[i] = integral_ro[i] - integral_rs[i];
      }
    }
#ifdef PRINT_DEBUG
    fmt::println("I2: {}, {}, {}", radial_integrals[0], radial_integrals[1], radial_integrals[2]);
#endif
  }
};