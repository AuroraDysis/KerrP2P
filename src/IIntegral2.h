#pragma once

#include "Common.h"

// Radial Antiderivatives for case (2)
template<typename Real, typename Complex>
class IIntegral2 {
private:
  ForwardRayTracing<Real, Complex> &data;

  Real asin_x2_rs, asin_x2_inf, k;
  Real E2_Numerator, F2_Denominator, Pi_p2_T1, Pi_m2_T1, Pi_p2_m2_Numerator, Pi_p2_Denominator, Pi_m2_Denominator;

  Real F2, E2, Pi_p2, Pi_m2, I_p, I_m;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_inf;
public:
  explicit IIntegral2(ForwardRayTracing<Real, Complex> &parent) : data(parent) {
  }

  void pre_calc() {
    const Real &rp = data.rp;
    const Real &rm = data.rm;

    const Real &r1 = data.r1;
    const Real &r2 = data.r2;
    const Real &r3 = data.r3;
    const Real &r4 = data.r4;
    const Real &r_s = data.r_s;

    asin_x2_rs = sqrt(((r1 - r3) * (r_s - r4)) / ((r_s - r3) * (r1 - r4)));
    asin_x2_inf = sqrt((r1 - r3) / (r1 - r4));
    k = ((-r2 + r3) * (-r1 + r4)) / ((r1 - r3) * (r2 - r4));

    E2_Numerator = sqrt((r1 - r3) * (r2 - r4));
    F2_Denominator = 1 / E2_Numerator;
    Pi_p2_T1 = ((r1 - r4) * (r3 - rp)) / ((r1 - r3) * (r4 - rp));
    Pi_m2_T1 = ((r1 - r4) * (r3 - rm)) / ((r1 - r3) * (r4 - rm));
    Pi_p2_m2_Numerator = 2 * (-r3 + r4);
    Pi_p2_Denominator = 1 / (E2_Numerator * (-r3 + rp) * (-r4 + rp));
    Pi_m2_Denominator = 1 / (E2_Numerator * (-r3 + rm) * (-r4 + rm));
  }

  void calc_x(std::array<Real, 3>& integral, const Real &asin_x2) {
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
    if (data.calc_t_f) {
      // I_t
      // integral[2] = (a * a * (I_p * (a * lambda - 2 * rp) + I_m * (2 * rm - a * lambda))) / (rm - rp);
    } else {
      integral[2] = std::numeric_limits<Real>::quiet_NaN();
    }
  }

  void calc(bool is_plus) {
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

#ifdef PRINT_DEBUG
    fmt::println("I2: {}, {}, {}", radial_integrals[0], radial_integrals[1], radial_integrals[2]);
#endif
  }
};