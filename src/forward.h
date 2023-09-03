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

class ForwardRayTracing;

class IIntegral {
private:
  ForwardRayTracing &data;
public:
  explicit IIntegral(ForwardRayTracing &parent) : data(parent) {
  }

  void reinit();

  void calc();
};

class GIntegral {
private:
  std::array<Real, 3> G_theta_p = {};
  std::array<Real, 3> G_theta_s = {};
  std::array<Real, 3> G_theta_inf = {};

  Real one_over_sqrt_up;
  Real up_over_um;
  Real one_over_umaa_sqrt;
  Real asin_up_cos_theta; // ArcCsc[Sqrt[up] Sec[\[Theta]]]

  ForwardRayTracing &data;
  void G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta);
public:
  explicit GIntegral(ForwardRayTracing &parent) : data(parent) {
  }

  void reinit();
  void calc();
};

class ForwardRayTracing {
  friend class IIntegral;
  friend class GIntegral;

private:
  const Real a, rp, rm, r_s, theta_s;

  // need to set before using other method
  Sign nu_r, nu_theta;
  Real lambda, q, eta;

  // auto initialized
  RayStatus ray_status = RayStatus::NORMAL;
  Real delta_theta, up, um, theta_p, theta_m;
  Complex r1_c, r2_c, r3_c, r4_c;
  Real r1, r2, r3, r4;
  bool r12_is_real, r34_is_real;

  // minor time
  Real tau_o;
  // I_theta, I_phi, I_t
  std::array<Real, 3> radial_integrals;
  // G_theta, G_phi, G_t
  std::array<Real, 3> angular_integrals;

  void reset_variables() {
    tau_o = std::numeric_limits<Real>::quiet_NaN();
    std::fill(radial_integrals.begin(), radial_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
    std::fill(angular_integrals.begin(), angular_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
  }

  // Range of theta: For Type A (eta > 0)
  void init_theta_pm() {
    delta_theta = half<Real>() * (1 - (eta + lambda * lambda) / (a * a));
    up = delta_theta + mp::sqrt(delta_theta * delta_theta + eta / (a * a));
    um = delta_theta - mp::sqrt(delta_theta * delta_theta + eta / (a * a));
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

  void reset_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_) {
    ray_status = RayStatus::NORMAL;
    lambda = std::move(lambda_);
    q = std::move(q_);
    eta = q * q;
    nu_r = nu_r_;
    nu_theta = nu_theta_;
    init_radial_potential_roots();
    init_theta_pm();
    reset_variables();
    I_integral->reinit();
    G_integral->reinit();
  }

  // convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
  void reset_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_) {
    ray_status = RayStatus::NORMAL;
    Real lambda_c = a + (rc * (2 * mp::pow(a, 2) + (-3 + rc) * rc)) / (a - a * rc);
    Real eta_c = -((mp::pow(rc, 3) * (-4 * mp::pow(a, 2) + mp::pow(-3 + rc, 2) * rc)) /
                   (mp::pow(a, 2) * mp::pow(-1 + rc, 2)));
    Real qc = mp::sqrt(eta_c);

    Real coeff = mp::sqrt(
        mp::pow(-3 + rc, 2) / (mp::pow(a, 2) * mp::pow(-1 + rc, 2)) + eta_c / mp::pow(rc, 4));

    reset_by_lambda_q(lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff),
                      qc + d * (mp::sqrt(eta_c) / mp::pow(rc, 2) / coeff), nu_r_, nu_theta_);
    reset_variables();
  }

public:
  Real theta_inf = std::numeric_limits<Real>::quiet_NaN();
  Real phi_inf = std::numeric_limits<Real>::quiet_NaN();
  Real t_inf = std::numeric_limits<Real>::quiet_NaN();
  int m = 0;
  Real n_half = std::numeric_limits<Real>::quiet_NaN();

  // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
  ForwardRayTracing(Real a_, Real r_s_, Real theta_s_)
      : a(std::move(a_)), r_s(std::move(r_s_)), theta_s(std::move(theta_s_)),
        rp(1 + mp::sqrt(1 - a * a)), rm(1 - mp::sqrt(1 - a * a)) {
    I_integral = std::make_shared<IIntegral>(*this);
    G_integral = std::make_shared<GIntegral>(*this);
  }

  std::shared_ptr<IIntegral> I_integral;
  std::shared_ptr<GIntegral> G_integral;

  RayStatus calc_ray_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_) {
    reset_by_lambda_q(std::move(lambda_), std::move(q_), nu_r_, nu_theta_);
    return calc_ray();
  }

  RayStatus calc_ray_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_) {
    reset_by_rc_d(rc, d, nu_r_, nu_theta_);
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
    I_integral->calc();

    tau_o = radial_integrals[0];

    if (ray_status != RayStatus::NORMAL) {
      return ray_status;
    }

    G_integral->calc();

    // Final values of phi and t
    phi_inf = radial_integrals[1] + lambda * angular_integrals[1];
    // t_inf = radial_integrals[2] + a * a * angular_integrals[1];

    std::cout << "theta_inf: " << theta_inf << ", phi_inf: " << phi_inf
              << ", m: "
              << m << ", nhalf: " << n_half << std::endl;
    return RayStatus::NORMAL;
  }
};
