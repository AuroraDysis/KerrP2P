#pragma once

#include "Common.h"
#include "IIntegral2.h"
#include "IIntegral3.h"
#include "GIntegral.h"

#include <boost/numeric/conversion/converter.hpp>

template<typename Real, typename Complex>
class ForwardRayTracing {
  friend class IIntegral2<Real, Complex>;

  friend class IIntegral3<Real, Complex>;

  friend class GIntegral<Real, Complex>;

#ifdef TESTS
public:
#else
private:
#endif

  const Real a, rp, rm, r_s, theta_s, r_o;

  // need to set before using other method
  Sign nu_r, nu_theta;
  Real lambda, q, eta;

  // auto initialized
  RayStatus ray_status = RayStatus::NORMAL;
  Real delta_theta, up, um, theta_p, theta_m;
  Complex r1_c, r2_c, r3_c, r4_c;
  Real r1, r2, r3, r4;
  bool r12_is_real, r34_is_real;

  void reset_variables() {
    tau_o = std::numeric_limits<Real>::quiet_NaN();
    theta_f = std::numeric_limits<Real>::quiet_NaN();
    phi_f = std::numeric_limits<Real>::quiet_NaN();
    t_f = std::numeric_limits<Real>::quiet_NaN();
    m = std::numeric_limits<int>::max();
    n_half = std::numeric_limits<Real>::quiet_NaN();
    std::fill(radial_integrals.begin(), radial_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
    std::fill(angular_integrals.begin(), angular_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
  }

  void init_radial_potential_roots() {
    Real AA = a * a - eta - lambda * lambda;
    Real BB = 2 * (eta + (lambda - a) * (lambda - a));
    Real CC = -a * a * eta;
    Real PP = -AA * AA / 12 - CC;
    Real QQ = (-2 * cube(AA) - 27 * square(BB) + 72 * AA * CC) / 216;

    Complex omega_pm;
    Real omega_pm_1 = -QQ * half<Real>();
    Real omega_pm_2 = cube(PP) / 27 + square(QQ) / 4;

    if (omega_pm_2 > 0) {
      omega_pm = cbrt(omega_pm_1 + sqrt(omega_pm_2)) +
                 cbrt(omega_pm_1 - sqrt(omega_pm_2));
    } else {
      Complex omega_pm_2_c = Complex{omega_pm_2};
      omega_pm = pow(omega_pm_1 + sqrt(omega_pm_2_c), third<Real>()) +
                 pow(omega_pm_1 - sqrt(omega_pm_2_c), third<Real>());
    }

    Real z = sqrt((real(omega_pm) - AA * third<Real>()) * half<Real>());

    Complex sqrt_in_1 = -(AA * half<Real>()) - square(z) + BB / (4 * z);
    Complex sqrt_in_2 = -(AA * half<Real>()) - square(z) - BB / (4 * z);

    r1_c = -z - sqrt(sqrt_in_1);
    r2_c = -z + sqrt(sqrt_in_1);
    r3_c = z - sqrt(sqrt_in_2);
    r4_c = z + sqrt(sqrt_in_2);

    if (real(sqrt_in_1) < 0) {
      r12_is_real = false;
      r1 = r2 = std::numeric_limits<Real>::quiet_NaN();
    } else {
      r12_is_real = true;
      r1 = real(r1_c);
      r2 = real(r2_c);
    }

    if (real(sqrt_in_2) < 0) {
      r34_is_real = false;
      r3 = r4 = std::numeric_limits<Real>::quiet_NaN();
    } else {
      r34_is_real = true;
      r3 = real(r3_c);
      r4 = real(r4_c);
    }

#ifdef PRINT_DEBUG
    fmt::println("r1: {}, r2: {}, r3: {}, r4: {}", r1, r2, r3, r4);
#endif
  }

  void init_theta_pm() {
    delta_theta = half<Real>() * (1 - (eta + lambda * lambda) / (a * a));
    up = delta_theta + sqrt(delta_theta * delta_theta + eta / (a * a));
    um = delta_theta - sqrt(delta_theta * delta_theta + eta / (a * a));
    theta_p = acos(-sqrt(up));
    theta_m = acos(sqrt(up));
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
#ifdef PRINT_DEBUG
    fmt::println("lambda: {}, q: {}, nu_r: {}, nu_theta: {}", lambda, q, to_integral(nu_r), to_integral(nu_theta));
#endif
  }

  void reset_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_) {
    ray_status = RayStatus::NORMAL;
    Real lambda_c = a + (rc * (2 * square(a) + (-3 + rc) * rc)) / (a - a * rc);
    Real eta_c = -((cube(rc) * (-4 * square(a) + square(-3 + rc) * rc)) /
                   (square(a) * square(-1 + rc)));
    Real qc = sqrt(eta_c);

    Real coeff = sqrt(
        square(-3 + rc) / (square(a) * square(-1 + rc)) + eta_c / pow(rc, 4));

    reset_by_lambda_q(lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff),
                      qc + d * (sqrt(eta_c) / square(rc) / coeff), nu_r_, nu_theta_);
    reset_variables();
  }

  void calcI() {
    bool radial_turning = r34_is_real && r4 > rp;

    // if there is a radial turning point (i.e. r4 is a real number)
    if (radial_turning && r_s <= r4) {
      ray_status = RayStatus::CONFINED;
      return;
    }

    if (radial_turning && r_s > r4 && nu_r == Sign::POSITIVE) {
      I_integral_2->calc(false);
      return;
    }

    if (radial_turning && r_s > r4 && nu_r == Sign::NEGATIVE) {
      I_integral_2->calc(true);
      return;
    }

    if (!radial_turning && nu_r == Sign::NEGATIVE) {
      ray_status = RayStatus::FALLS_IN;
      return;
    }

    if (!radial_turning && nu_r == Sign::POSITIVE) {
      I_integral_3->calc(false);
      return;
    }

    ray_status = RayStatus::UNKOWN_ERROR;
  }

public:
  // minor time
  Real tau_o;
  // I_theta, I_phi, I_t
  std::array<Real, 3> radial_integrals;
  // G_theta, G_phi, G_t
  std::array<Real, 3> angular_integrals;

  Real theta_f;
  Real phi_f;
  Real t_f;
  int m;
  Real n_half;
  bool calc_t_f = false;

  // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
  ForwardRayTracing<Real, Complex>(Real a_, Real r_s_, Real theta_s_, Real r_o_)
      : a(std::move(a_)), r_s(std::move(r_s_)), theta_s(std::move(theta_s_)), r_o(std::move(r_o_)),
        rp(1 + sqrt(1 - a * a)), rm(1 - sqrt(1 - a * a)) {
    reset_variables();
    I_integral_2 = std::make_shared<IIntegral2<Real, Complex>>(*this);
    I_integral_3 = std::make_shared<IIntegral3<Real, Complex>>(*this);
    G_integral = std::make_shared<GIntegral<Real, Complex>>(*this);
#ifdef PRINT_DEBUG
    fmt::println("a: {}, r_s: {}, theta_s: {}, r_o: {}, rp: {}, rm: {}", a, r_s, theta_s, r_o, rp, rm);
#endif
  }

  std::shared_ptr<IIntegral2<Real, Complex>> I_integral_2;
  std::shared_ptr<IIntegral3<Real, Complex>> I_integral_3;
  std::shared_ptr<GIntegral<Real, Complex>> G_integral;

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
    calcI();

    tau_o = radial_integrals[0];

    if (ray_status != RayStatus::NORMAL) {
      return ray_status;
    }

    G_integral->calc();

    // Final values of phi and t
    phi_f = radial_integrals[1] + lambda * angular_integrals[1];
    if (calc_t_f) {
      t_f = radial_integrals[2] + square(a) * angular_integrals[2];
    }

#ifdef PRINT_DEBUG
    fmt::println("theta_f, phi_f, t_f, m, nhalf: {}, {}, {}, {}, {}", theta_f, phi_f, t_f, m, n_half);
#endif
    return RayStatus::NORMAL;
  }
};
