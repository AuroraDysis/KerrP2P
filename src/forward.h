#pragma once

#include <array>
#include <cmath>
#include <boost/numeric/conversion/converter.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

using boost::math::constants::half;
using boost::math::constants::third;
using boost::math::constants::sixth;

//using Real = double;
//using Complex = std::complex<Real>;
//namespace mp = std;
// using Real = long double;

//#include <boost/multiprecision/float128.hpp>
//#include <boost/multiprecision/complex128.hpp>
//using Real = boost::multiprecision::float128;
//using Complex = boost::multiprecision::complex128;
//namespace mp = boost::multiprecision;

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <utility>

using Real = boost::multiprecision::mpfr_float_50;
using Complex = boost::multiprecision::mpc_complex_50;
namespace mp = boost::multiprecision;

enum class RayStatus {
  NORMAL,
  FALLS_IN,
  CONFINED,
  ETA_OUT_OF_RANGE, // eta should be positive
  THETA_OUT_OF_RANGE, // theta should be in [theta_m, theta_p]
  UNKOWN_ERROR,
};

enum class Sign : int {
  POSITIVE = 1,
  NEGATIVE = -1,
};

template<typename E>
constexpr auto to_integral(E e) -> typename std::underlying_type<E>::type {
  return static_cast<typename std::underlying_type<E>::type>(e);
}

class ForwardRayTracing {
private:
  const Real a, rp, rm, r_s, theta_s;
  const Sign nu_r, nu_theta;
  // need to set before using other method
  Real lambda, q, eta;
  // auto initialized
  RayStatus ray_status = RayStatus::NORMAL;
  Real delta_theta, up, um, theta_p, theta_m;
  Complex r1_c, r2_c, r3_c, r4_c;
  Real r1, r2, r3, r4;
  bool r12_is_real, r34_is_real;
  // radial coeffs
  Real A, B, alpha_0, alpha_p, alpha_m, k3;
  // minor time
  Real tau_o;
  // I_theta, I_phi, I_t
  std::array<Real, 3> radial_integrals = {std::numeric_limits<Real>::quiet_NaN(),
                                          std::numeric_limits<Real>::quiet_NaN(),
                                          std::numeric_limits<Real>::quiet_NaN()};
  // G_theta, G_phi, G_t
  std::array<Real, 3> angular_integrals = {std::numeric_limits<Real>::quiet_NaN(),
                                           std::numeric_limits<Real>::quiet_NaN(),
                                           std::numeric_limits<Real>::quiet_NaN()};

  // Range of theta: For Type A (eta > 0)
  void init_theta_pm() {
    delta_theta = half<Real>() * (1 - (eta + lambda * lambda) / (a * a));
    up = delta_theta + mp::sqrt(delta_theta * delta_theta + eta / (a * a));
    // um = delta_theta - mp::sqrt(delta_theta * delta_theta + eta / (a * a));
    theta_p = mp::acos(-mp::sqrt(up));
    theta_m = mp::acos(mp::sqrt(up));
  }

  // find roots
  void init_radial_potential_roots() {
    Real AA = a * a - eta - lambda * lambda;
    Real BB = 2 * (eta + (lambda - a) * (lambda - a));
    Real CC = -a * a * eta;
    Real PP = -AA * AA / 12 - CC;
    Real QQ = (-2 * mp::pow(AA, 3) - 27 * mp::pow(BB, 2) + 72 * AA * CC) / 216;

    Complex omega_pm;
    Real omega_pm_1 = -QQ * half<Real>();
    Real omega_pm_2 = mp::pow(PP, 3) / 27 + mp::pow(QQ, 2) / 4;

    if (omega_pm_2 > 0) {
      omega_pm = mp::cbrt(omega_pm_1 + mp::sqrt(omega_pm_2)) +
                 mp::cbrt(omega_pm_1 - mp::sqrt(omega_pm_2));
    } else {
      Complex omega_pm_2_c = Complex{omega_pm_2};
      omega_pm = mp::pow(omega_pm_1 + mp::sqrt(omega_pm_2_c), third<Real>()) +
                 mp::pow(omega_pm_1 - mp::sqrt(omega_pm_2_c), third<Real>());
    }

    Real z = mp::sqrt((mp::real(omega_pm) - AA * third<Real>()) * half<Real>());

    Complex sqrt_in_1 = -(AA * half<Real>()) - z * z + BB / (4 * z);
    Complex sqrt_in_2 = -(AA * half<Real>()) - z * z - BB / (4 * z);

    r1_c = -z - mp::sqrt(sqrt_in_1);
    r2_c = -z + mp::sqrt(sqrt_in_1);
    r3_c = z - mp::sqrt(sqrt_in_2);
    r4_c = z + mp::sqrt(sqrt_in_2);

    if (mp::real(sqrt_in_1) < 0) {
      r12_is_real = false;
      r1 = r2 = std::numeric_limits<Real>::quiet_NaN();
    } else {
      r12_is_real = true;
      r1 = mp::real(r1_c);
      r2 = mp::real(r2_c);
    }

    if (mp::real(sqrt_in_2) < 0) {
      r34_is_real = false;
      r3 = r4 = std::numeric_limits<Real>::quiet_NaN();
    } else {
      r34_is_real = true;
      r3 = mp::real(r3_c);
      r4 = mp::real(r4_c);
    }
  }

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

  void reset_by_lambda_q(Real lambda_, Real q_) {
    ray_status = RayStatus::NORMAL;
    lambda = std::move(lambda_);
    q = std::move(q_);
    eta = q * q;
    init_radial_potential_roots();
    init_theta_pm();
  }

  // convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
  void reset_by_rc_d(const Real &rc, const Real &d) {
    ray_status = RayStatus::NORMAL;
    Real lambda_c = a + (rc * (2 * mp::pow(a, 2) + (-3 + rc) * rc)) / (a - a * rc);
    Real eta_c = -((mp::pow(rc, 3) * (-4 * mp::pow(a, 2) + mp::pow(-3 + rc, 2) * rc)) /
                   (mp::pow(a, 2) * mp::pow(-1 + rc, 2)));
    Real qc = mp::sqrt(eta_c);

    Real coeff = mp::sqrt(
        mp::pow(-3 + rc, 2) / (mp::pow(a, 2) * mp::pow(-1 + rc, 2)) + eta_c / mp::pow(rc, 4));

    reset_by_lambda_q(lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff),
                      qc + d * (mp::sqrt(eta_c) / mp::pow(rc, 2) / coeff));
  }

public:
  Real theta_inf = std::numeric_limits<Real>::quiet_NaN();
  Real phi_inf = std::numeric_limits<Real>::quiet_NaN();
  Real t_inf = std::numeric_limits<Real>::quiet_NaN();
  int m = 0;
  Real n_half = std::numeric_limits<Real>::quiet_NaN();

  // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
  ForwardRayTracing(Real a_, Real r_s_, Real theta_s_, Sign nu_r_, Sign nu_theta_)
      : a(std::move(a_)), r_s(std::move(r_s_)), theta_s(std::move(theta_s_)), nu_r(nu_r_), nu_theta(nu_theta_),
        rp(1 + mp::sqrt(1 - a * a)), rm(1 - mp::sqrt(1 - a * a)) {
  }

  Real I_r(Real r) const {
    return I_0(r);
  }

  Real I_phi(Real r) const {
    return (2 * a) / (rp - rm) * ((rp - (a * lambda) / 2) * I_p(r) - (rm - (a * lambda) / 2) * I_m(r));
  }

  Real I_t(Real r) const {
    return 4 / (rp - rm) * (rp * (rp - (a * lambda) / 2) * I_p(r) - rm * (rm - (a * lambda) / 2) * I_m(r))
           + 4 * I_0(r) + 2 * I_1(r) + I_2(r);
  }

  void calc_I() {
    A = mp::sqrt((r3 - r2) * (r4 - r2));
    B = mp::sqrt((r3 - r1) * (r4 - r1));

    alpha_0 = (B + A) / (B - A);
    alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
    alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));

    k3 = (mp::pow(A + B, 2) - mp::pow(r2 - r1, 2)) / (4 * A * B);

    std::array<std::function<Real(Real)>, 3> anti_ders = {
        [&](Real r) -> Real { return I_r(r); },
        [&](Real r) -> Real { return I_phi(r); },
        [&](Real r) -> Real { return I_t(r); }
    };

    std::array<Real, 3> results = {};
    std::fill(results.begin(), results.end(), std::numeric_limits<Real>::quiet_NaN());
    bool radial_turning = r34_is_real && r4 > rp;

    // if there is a radial turning point (i.e. r4 is a real number)
    if (radial_turning && r_s <= r4) {
      ray_status = RayStatus::CONFINED;
      return;
    }

    if (radial_turning && r_s > r4 && nu_r == Sign::POSITIVE) {
      for (int i = 0; i < 3; i++) {
        results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) - anti_ders[i](r_s);
      }
      return;
    }

    if (radial_turning && r_s > r4 && nu_r == Sign::NEGATIVE) {
      for (int i = 0; i < 3; i++) {
        results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) + anti_ders[i](r_s) - 2 * anti_ders[i](r4);
      }
      return;
    }

    if (!radial_turning && nu_r == Sign::NEGATIVE) {
      ray_status = RayStatus::FALLS_IN;
      return;
    }

    if (!radial_turning && nu_r == Sign::POSITIVE) {
      for (int i = 0; i < 3; i++) {
        results[i] = anti_ders[i](std::numeric_limits<Real>::infinity()) - anti_ders[i](r_s);
      }
      return;
    }

    ray_status = RayStatus::UNKOWN_ERROR;
  }

  void calc_G() const {
    std::array<Real, 3> result = {};
    result[0] = m * (G_theta(theta_p) - G_theta(theta_m)) +
                to_integral(nu_theta) * (pow(-1, m) * G_theta(theta_inf) - G_theta(theta_s));
    result[1] = m * (G_phi(theta_p) - G_phi(theta_m)) +
                to_integral(nu_theta) * (pow(-1, m) * G_phi(theta_inf) - G_phi(theta_s));
    result[2] = m * (G_t(theta_p) - G_t(theta_m)) +
                to_integral(nu_theta) * (pow(-1, m) * G_t(theta_inf) - G_t(theta_s));

    G_theta = -1 / mp::sqrt(-um * a * a) *
              boost::math::ellint_1(mp::asin(mp::cos(theta) / mp::sqrt(up)), up / um);
    G_phi = -1 / mp::sqrt(-um * a * a) *
            boost::math::ellint_3(up, mp::asin(mp::cos(theta) / mp::sqrt(up)), up / um);
    G_t = (2 * up) / mp::sqrt(-um * a * a) / (2 * up / um) *
          (boost::math::ellint_2(mp::asin(mp::cos(theta) / mp::sqrt(up)), up / um) -
           boost::math::ellint_1(theta));

    Real G_theta_theta_s = G_theta(theta_s);
    Real G_theta_theta_p = G_theta(theta_p);
    Real G_theta_theta_m = G_theta(theta_m);

    theta_inf = mp::acos(-mp::sqrt(up) * to_integral(nu_theta) *
                         boost::math::jacobi_sn(
                             mp::sqrt(-um * a * a) * (tau_o + to_integral(nu_theta) * G_theta_theta_s),
                             up / um));

    // Angular integrals
    Real m_Real = 1 + mp::floor(mp::real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                         (G_theta_theta_p - G_theta_theta_m)));

    using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
        boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;
    // floor
    int m = RealToInt::convert(m_Real);

    calc_G(theta_inf, m);

    // Number of half-orbits
    n_half = tau_o / (G_theta_theta_p - G_theta_theta_m);
  }

  RayStatus calc_ray_by_lambda_q(Real lambda_, Real q_) {
    reset_by_lambda_q(std::move(lambda_), std::move(q_));
    return calc_ray();
  }

  RayStatus calc_ray_by_rc_d(const Real &rc, const Real &d) {
    reset_by_rc_d(rc, d);
    return calc_ray();
  }

  RayStatus calc_ray() {
    if (eta <= 0) {
      return RayStatus::ETA_OUT_OF_RANGE;
    }

    if (theta_s < theta_m || theta_s > theta_p) {
      return RayStatus::THETA_OUT_OF_RANGE;
    }

    // Radial integrals
    calc_I();

    tau_o = radial_integrals[0];

    if (ray_status != RayStatus::NORMAL) {
      return ray_status;
    }

    calc_G();

    // Final values of phi and t
    phi_inf = radial_integrals[1] + lambda * angular_integrals[1];
    // t_inf = radial_integrals[2] + a * a * angular_integrals[1];

    std::cout << "theta_inf: " << theta_inf << ", phi_inf: " << phi_inf
              << ", m: "
              << m << ", nhalf: " << n_half << std::endl;
    return RayStatus::NORMAL;
  }
};
