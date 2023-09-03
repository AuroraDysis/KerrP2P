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

#ifdef FLOAT64
using Real = double;
using Complex = std::complex<Real>;
namespace mp = std;
#endif

#ifdef FLOAT128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
using Real = boost::multiprecision::float128;
using Complex = boost::multiprecision::complex128;
namespace mp = boost::multiprecision;
#endif

#ifdef BIGFLOAT
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
using Real = boost::multiprecision::mpfr_float_50;
using Complex = boost::multiprecision::mpc_complex_50;
namespace mp = boost::multiprecision;
#endif

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

// Radial Antiderivatives for case (2)
class IIntegral2 {
private:
  ForwardRayTracing &data;

  Real asin_x2_rs, asin_x2_inf, k;
  Real E2_Numerator, F2_Denominator, Pi_p2_T1, Pi_m2_T1, Pi_p2_m2_Numerator, Pi_p2_Denominator, Pi_m2_Denominator;

  Real F2, E2, Pi_p2, Pi_m2, I_p, I_m;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_inf;
public:
  explicit IIntegral2(ForwardRayTracing &parent) : data(parent) {
  }

  void pre_calc();

  void calc_x(std::array<Real, 3>& integral, const Real &x);

  void calc(bool is_plus);
};

// // Radial Antiderivatives for case (3)
class IIntegral3 {
private:
  ForwardRayTracing &data;

  Real acos_x3_rs, acos_x3_inf;
  Real r34_re, r34_im;
  Real A, B, alpha_p, alpha_p2, alpha_m, alpha_m2, k3;
  Real F3, R1_alpha_p, R1_alpha_m, Ip, Im;

  std::array<Real, 3> integral_rs;
  std::array<Real, 3> integral_inf;

  Real R1(const Real &acos_x3, const Real &alpha, const Real &alpha2) const;
public:
  explicit IIntegral3(ForwardRayTracing &parent) : data(parent) {
  }

  void pre_calc();

  void calc_x(std::array<Real, 3> &integral, const Real &x);

  void calc(bool is_plus);
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

  void pre_calc();

  void calc();
};

class ForwardRayTracing {
  friend class IIntegral2;

  friend class IIntegral3;

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

  void reset_variables();

  // Range of theta: For Type A (eta > 0)
  void init_theta_pm();

  // find roots
  void init_radial_potential_roots();

  void reset_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_);

  // convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
  void reset_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_);

  void calcI();

public:
  Real theta_f;
  Real phi_f;
  Real t_f;
  int m;
  Real n_half;
  bool calc_t_f = false;

  // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
  ForwardRayTracing(Real a_, Real r_s_, Real theta_s_)
      : a(std::move(a_)), r_s(std::move(r_s_)), theta_s(std::move(theta_s_)),
        rp(1 + mp::sqrt(1 - a * a)), rm(1 - mp::sqrt(1 - a * a)) {
    reset_variables();
    I_integral_2 = std::make_shared<IIntegral2>(*this);
    I_integral_3 = std::make_shared<IIntegral3>(*this);
    G_integral = std::make_shared<GIntegral>(*this);
  }

  std::shared_ptr<IIntegral2> I_integral_2;
  std::shared_ptr<IIntegral3> I_integral_3;
  std::shared_ptr<GIntegral> G_integral;

  RayStatus calc_ray_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_) {
    reset_by_lambda_q(std::move(lambda_), std::move(q_), nu_r_, nu_theta_);
    return calc_ray();
  }

  RayStatus calc_ray_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_) {
    reset_by_rc_d(rc, d, nu_r_, nu_theta_);
    return calc_ray();
  }

  RayStatus calc_ray();
};
