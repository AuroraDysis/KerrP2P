#pragma once

#include "Common.h"
#include "IIntegral2.h"
#include "IIntegral3.h"
#include "GIntegral.h"
#include "ObjectPool.h"

template<typename Real>
struct ForwardRayTracingParams {
  Real a;
  Real r_s;
  Real theta_s;
  Real r_o;
  Sign nu_r;
  Sign nu_theta;

  Real rc, lgd;
  Sign lgd_sign;

  Real lambda, q;

  bool calc_t_f = false;

  ForwardRayTracingParams() = default;

  ForwardRayTracingParams(const ForwardRayTracingParams &params) {
    a = params.a;
    r_s = params.r_s;
    theta_s = params.theta_s;
    r_o = params.r_o;
    nu_r = params.nu_r;
    nu_theta = params.nu_theta;
    rc = params.rc;
    lgd = params.lgd;
    lgd_sign = params.lgd_sign;
    lambda = params.lambda;
    q = params.q;
    calc_t_f = params.calc_t_f;
  }

  void rc_d_to_lambda_q() {
    Real r_up = 4 * MY_SQUARE(cos(acos(a) * third<Real>()));
    Real r_down = 4 * MY_SQUARE(cos(acos(-a) * third<Real>()));

    if (rc < r_down || rc > r_up) {
	  fmt::println("rc out of range: rc = {}, r_down: {}, r_up: {}", rc, r_down, r_up);
      lambda = std::numeric_limits<Real>::quiet_NaN();
      q = std::numeric_limits<Real>::quiet_NaN();
      return;
	}

    Real lambda_c = a + (rc * (2 * MY_SQUARE(a) + (-3 + rc) * rc)) / (a * (1 - rc));
    Real eta_c = -((MY_CUBE(rc) * (-4 * MY_SQUARE(a) + MY_SQUARE(-3 + rc) * rc)) /
                   (MY_SQUARE(a) * MY_SQUARE(-1 + rc)));
    Real qc = sqrt(eta_c);

    Real coeff = sqrt(
        MY_SQUARE(-3 + rc) / (MY_SQUARE(a) * MY_SQUARE(-1 + rc)) + eta_c / pow(rc, 4));

    Real d = to_integral(lgd_sign) * pow(static_cast<Real>(10), lgd);
    lambda = lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff);
    q = qc + d * (sqrt(eta_c) / MY_SQUARE(rc) / coeff);
#ifdef PRINT_DEBUG
    fmt::println("rc: {}, lgd: {}", rc, lgd);
    fmt::println("lambda_c: {}, eta_c: {}, qc: {}", lambda_c, eta_c, qc);
    fmt::println("coeff: {}", coeff);
#endif
  }
};

template<typename Real, typename Complex>
struct ForwardRayTracingResult {
  Real a, rp, rm, r_s, theta_s, r_o;
  Real r1, r2, r3, r4;
  Complex r1_c, r2_c, r3_c, r4_c;
  Real t_f, theta_f, phi_f;
  int m;
  Real n_half;
  Real eta, lambda, q;
  Real rc, lgd;
  Sign lgd_sign;
  RayStatus ray_status;
};

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
  inline static ObjectPool<ForwardRayTracing<Real, Complex>> pool;

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
    Real QQ = (-2 * MY_CUBE(AA) - 27 * MY_SQUARE(BB) + 72 * AA * CC) / 216;

    Complex omega_pm;
    Real omega_pm_1 = -QQ * half<Real>();
    Real omega_pm_2 = MY_CUBE(PP) / 27 + MY_SQUARE(QQ) / 4;

    if (omega_pm_2 > 0) {
      omega_pm = cbrt(omega_pm_1 + sqrt(omega_pm_2)) +
                 cbrt(omega_pm_1 - sqrt(omega_pm_2));
    } else {
      Complex omega_pm_2_c = Complex{omega_pm_2};
      omega_pm = pow(omega_pm_1 + sqrt(omega_pm_2_c), third<Real>()) +
                 pow(omega_pm_1 - sqrt(omega_pm_2_c), third<Real>());
    }

    Real z = sqrt((real(omega_pm) - AA * third<Real>()) * half<Real>());

    Complex sqrt_in_1 = -(AA * half<Real>()) - MY_SQUARE(z) + BB / (4 * z);
    Complex sqrt_in_2 = -(AA * half<Real>()) - MY_SQUARE(z) - BB / (4 * z);

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
    fmt::println("AA: {}, BB: {}, CC: {}, PP: {}, QQ: {}", AA, BB, CC, PP, QQ);
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
  Real a, rp, rm, r_s, theta_s, r_o;

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

  Real theta_f;
  Real phi_f;
  Real t_f;
  int m;
  Real n_half;
  bool calc_t_f = false;

  ForwardRayTracing<Real, Complex>() {
    reset_variables();
    I_integral_2 = std::make_shared<IIntegral2<Real, Complex>>(*this);
    I_integral_3 = std::make_shared<IIntegral3<Real, Complex>>(*this);
    G_integral = std::make_shared<GIntegral<Real, Complex>>(*this);
  }

  static std::shared_ptr<ForwardRayTracing<Real, Complex>> get_from_cache() {
    return pool.create();
  }

  static void clear_cache() {
    pool.clear();
  }

  std::shared_ptr<IIntegral2<Real, Complex>> I_integral_2;
  std::shared_ptr<IIntegral3<Real, Complex>> I_integral_3;
  std::shared_ptr<GIntegral<Real, Complex>> G_integral;

  void calc_ray(const ForwardRayTracingParams<Real> &params) {
    a = params.a;
    r_s = params.r_s;
    theta_s = params.theta_s;
    r_o = params.r_o;

    rp = 1 + sqrt(1 - a * a);
    rm = 1 - sqrt(1 - a * a);

    calc_t_f = params.calc_t_f;

#ifdef PRINT_DEBUG
    fmt::println("a: {}, r_s: {}, theta_s: {}, r_o: {}", a, r_s, theta_s, r_o);
    fmt::println("rp: {}, rm: {}", rp, rm);
    fmt::println("lambda: {}, q: {}", params.lambda, params.q);
#endif

    if (isnan(params.lambda) || isnan(params.q) || abs(params.lambda) <= 10000 * ErrorLimit<Real>::Value) {
      ray_status = RayStatus::ARGUMENT_ERROR;
      return;
    }

    reset_by_lambda_q(params.lambda, params.q, params.nu_r, params.nu_theta);

    if (eta <= 0) {
      ray_status = RayStatus::ETA_OUT_OF_RANGE;
      return;
    }

    if (theta_s < theta_m || theta_s > theta_p) {
      ray_status = RayStatus::THETA_OUT_OF_RANGE;
      return;
    }

    // Radial integrals
    calcI();
    if (ray_status != RayStatus::NORMAL) {
      return;
    }

    tau_o = radial_integrals[0];

    if (ray_status != RayStatus::NORMAL) {
      return;
    }

    G_integral->calc();
    if (ray_status != RayStatus::NORMAL) {
      return;
    }

    // Final values of phi and t
    phi_f = radial_integrals[1] + lambda * angular_integrals[1];
    if (calc_t_f) {
      t_f = radial_integrals[2] + MY_SQUARE(a) * angular_integrals[2];
    }

#ifdef PRINT_DEBUG
    fmt::println("theta_f, phi_f, t_f, m, nhalf: {}, {}, {}, {}, {}", theta_f, phi_f, t_f, m, n_half);
#endif
  }

  ForwardRayTracingResult<Real, Complex> to_result() {
    ForwardRayTracingResult<Real, Complex> result;
    result.a = a;
    result.rp = rp;
    result.rm = rm;
    result.r_s = r_s;
    result.theta_s = theta_s;
    result.r_o = r_o;
    result.r1 = r1;
    result.r2 = r2;
    result.r3 = r3;
    result.r4 = r4;
    result.r1_c = r1_c;
    result.r2_c = r2_c;
    result.r3_c = r3_c;
    result.r4_c = r4_c;
    result.t_f = t_f;
    result.theta_f = theta_f;
    result.phi_f = phi_f;
    result.m = m;
    result.n_half = n_half;
    result.ray_status = ray_status;
    result.eta = eta;
    result.lambda = lambda;
    result.q = q;
    return result;
  }
};
