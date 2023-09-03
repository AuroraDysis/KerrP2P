#include "forward.h"

#define square(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))

void IIntegral2::pre_calc() {
  const Real &rp = data.rp;
  const Real &rm = data.rm;

  const Real &r1 = data.r1;
  const Real &r2 = data.r2;
  const Real &r3 = data.r3;
  const Real &r4 = data.r4;
  const Real &r_s = data.r_s;

  asin_x2_rs = mp::sqrt(((r1 - r3) * (r_s - r4)) / ((r_s - r3) * (r1 - r4)));
  asin_x2_inf = mp::sqrt((r1 - r3) / (r1 - r4));
  k = ((-r2 + r3) * (-r1 + r4)) / ((r1 - r3) * (r2 - r4));

  E2_Numerator = mp::sqrt((r1 - r3) * (r2 - r4));
  F2_Denominator = 1 / E2_Numerator;
  Pi_p2_T1 = ((r1 - r4) * (r3 - rp)) / ((r1 - r3) * (r4 - rp));
  Pi_m2_T1 = ((r1 - r4) * (r3 - rm)) / ((r1 - r3) * (r4 - rm));
  Pi_p2_m2_Numerator = 2 * (-r3 + r4);
  Pi_p2_Denominator = 1 / (E2_Numerator * (-r3 + rp) * (-r4 + rp));
  Pi_m2_Denominator = 1 / (E2_Numerator * (-r3 + rm) * (-r4 + rm));
}

void IIntegral2::calc_x(std::array<Real, 3>& integral, const Real &asin_x2) {
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

void IIntegral2::calc(bool is_plus) {
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
}

void IIntegral3::pre_calc() {
  const Real &rp = data.rp;
  const Real &rm = data.rm;
  const Real &r_s = data.r_s;
  // two real roots, both inside horizon, r_1 < r_2 < r_- < r_+ and r_3 = conj(r_4)
  const Real &r1 = data.r1;
  const Real &r2 = data.r2;
  const Complex &r3 = data.r3_c;
  // const Complex &r4 = data.r4_c;

  r34_re = mp::real(r3);
  r34_im = mp::imag(r3);

  // radial coeffs
  A = mp::sqrt(square(r34_im) + square(r34_re - r2));
  B = mp::sqrt(square(r34_im) + square(r34_re - r1));

  k3 = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
  // alpha_0 = (B + A) / (B - A);
  alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
  alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));
  alpha_p2 = square(alpha_p);
  alpha_m2 = square(alpha_m);

  acos_x3_rs = mp::acos(-1 + (2 * A * (r_s - r1)) / (A * (r_s - r1) + B * (r_s - r2)));
  acos_x3_inf = mp::acos((A - B) / (A + B));
}


Real IIntegral3::R1(const Real &acos_x3, const Real &alpha, const Real &alpha2) const {
  return 1 / (1 - alpha2) *
         (boost::math::ellint_3(alpha2 / (alpha2 - 1), acos_x3, k3) -
          alpha * ((mp::sqrt((-1 + alpha2) /
                             (alpha2 + k3 - alpha2 * k3)) *
                    mp::log(mp::abs((mp::sin(acos_x3) +
                                     mp::sqrt((-1 + alpha2) /
                                              (alpha2 + k3 - alpha2 * k3)) *
                                     mp::sqrt(1 - k3 * square(mp::sin(acos_x3)))) /
                                    (-mp::sin(acos_x3) + mp::sqrt((-1 + alpha2) /
                                                                  (alpha2 + k3 - alpha2 * k3)) *
                                                         mp::sqrt(1 - k3 * square(mp::sin(acos_x3))))))) *
                   half<Real>()));
}

void IIntegral3::calc_x(std::array<Real, 3>& integral, const Real &acos_x3) {
  const Real &a = data.a;
  const Real &lambda = data.lambda;
  const Real &rp = data.rp;
  const Real &rm = data.rm;
  const Real &r1 = data.r1;
  const Real &r2 = data.r2;

  R1_alpha_p = R1(acos_x3, alpha_p, alpha_p2);
  R1_alpha_m = R1(acos_x3, alpha_m, alpha_m2);
  F3 = boost::math::ellint_1(acos_x3, k3) / mp::sqrt(A * B);
  Ip = -(((A + B) * F3 + (2 * mp::sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                         (A * (r1 - rp) + B * (-r2 + rp))) /
         (-(A * r1) - B * r2 + (A + B) * rp));
  Im = -(((A + B) * F3 + (2 * mp::sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                         (A * (r1 - rm) + B * (-r2 + rm))) /
         (-(A * r1) - B * r2 + (A + B) * rm));
  integral[0] = F3;
  integral[1] = (a * (Im * (-(a * lambda) + 2 * rm) + Ip * (a * lambda - 2 * rp))) / (rm - rp);
}

void IIntegral3::calc(bool is_plus) {
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

void GIntegral::pre_calc() {
  const Real &a = data.a;
  const Real &up = data.up;
  const Real &um = data.um;

  up_over_um = up / um;
  one_over_sqrt_up = mp::sqrt(up);
  one_over_umaa_sqrt = 1 / mp::sqrt(-um * a * a);

  G_theta_p[0] = -boost::math::ellint_1(up_over_um) * one_over_umaa_sqrt;
  G_theta_p[1] = -boost::math::ellint_3(up, up_over_um) * one_over_umaa_sqrt;
  if (data.calc_t_f) {
    G_theta_p[2] = um * (-boost::math::ellint_2(up_over_um) * one_over_umaa_sqrt + G_theta_p[0]);
  } else {
    G_theta_p[2] = std::numeric_limits<Real>::quiet_NaN();
  }
}

void GIntegral::G_theta_phi_t(std::array<Real, 3> &G_arr, const Real &theta) {
  const Real &up = data.up;
  const Real &um = data.um;
  asin_up_cos_theta = mp::asin(mp::cos(theta) * one_over_sqrt_up);
  G_arr[0] = -boost::math::ellint_1(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
  G_arr[1] = -boost::math::ellint_3(up, asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt;
  G_arr[2] = um * (boost::math::ellint_2(asin_up_cos_theta, up_over_um) * one_over_umaa_sqrt + G_arr[0]);
}

void GIntegral::calc() {
  const Real &a = data.a;
  const Real &up = data.up;
  const Real &um = data.um;
  auto &tau_o = data.tau_o;
  const Real &theta_s = data.theta_s;
  Sign nu_theta = data.nu_theta;

  G_theta_phi_t(G_theta_s, theta_s);

  const Real &G_theta_theta_s = G_theta_s[0];
  const Real &G_theta_theta_p = G_theta_p[0];

  data.theta_f = mp::acos(-mp::sqrt(up) * to_integral(nu_theta) *
                          boost::math::jacobi_sn(
                                (tau_o + to_integral(nu_theta) * G_theta_theta_s) / one_over_umaa_sqrt, up_over_um));

  // Angular integrals
  Real m_Real = 1 + mp::floor(mp::real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                       (2 * G_theta_theta_p)));

  using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
      boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;

  // floor
  data.m = RealToInt::convert(m_Real);

  // Number of half-orbits
  data.n_half = tau_o / (2 * G_theta_theta_p);

  G_theta_phi_t(G_theta_inf, data.theta_f);

  auto &angular_integrals = data.angular_integrals;
  // (-1)^m
  int m1_m = 1 - ((data.m & 1) << 1);

  for (int i = 0; i < 3; ++i) {
    angular_integrals[i] =
        (2 * data.m) * G_theta_p[i] + to_integral(nu_theta) * (m1_m * G_theta_inf[i] - G_theta_s[i]);
  }
}

void ForwardRayTracing::init_radial_potential_roots() {
  Real AA = a * a - eta - lambda * lambda;
  Real BB = 2 * (eta + (lambda - a) * (lambda - a));
  Real CC = -a * a * eta;
  Real PP = -AA * AA / 12 - CC;
  Real QQ = (-2 * cube(AA) - 27 * square(BB) + 72 * AA * CC) / 216;

  Complex omega_pm;
  Real omega_pm_1 = -QQ * half<Real>();
  Real omega_pm_2 = cube(PP) / 27 + square(QQ) / 4;

  if (omega_pm_2 > 0) {
    omega_pm = mp::cbrt(omega_pm_1 + mp::sqrt(omega_pm_2)) +
               mp::cbrt(omega_pm_1 - mp::sqrt(omega_pm_2));
  } else {
    Complex omega_pm_2_c = Complex{omega_pm_2};
    omega_pm = mp::pow(omega_pm_1 + mp::sqrt(omega_pm_2_c), third<Real>()) +
               mp::pow(omega_pm_1 - mp::sqrt(omega_pm_2_c), third<Real>());
  }

  Real z = mp::sqrt((mp::real(omega_pm) - AA * third<Real>()) * half<Real>());

  Complex sqrt_in_1 = -(AA * half<Real>()) - square(z) + BB / (4 * z);
  Complex sqrt_in_2 = -(AA * half<Real>()) - square(z) - BB / (4 * z);

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

void ForwardRayTracing::calcI() {
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

RayStatus ForwardRayTracing::calc_ray() {
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

  std::cout << "theta_f: " << theta_f << ", phi_f: " << phi_f << ", t_f: " << t_f
            << ", m: "
            << m << ", nhalf: " << n_half << std::endl;
  return RayStatus::NORMAL;
}

void ForwardRayTracing::reset_variables() {
  tau_o = std::numeric_limits<Real>::quiet_NaN();
  theta_f = std::numeric_limits<Real>::quiet_NaN();
  phi_f = std::numeric_limits<Real>::quiet_NaN();
  t_f = std::numeric_limits<Real>::quiet_NaN();
  m = std::numeric_limits<int>::max();
  n_half = std::numeric_limits<Real>::quiet_NaN();
  std::fill(radial_integrals.begin(), radial_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
  std::fill(angular_integrals.begin(), angular_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
}

void ForwardRayTracing::init_theta_pm() {
  delta_theta = half<Real>() * (1 - (eta + lambda * lambda) / (a * a));
  up = delta_theta + mp::sqrt(delta_theta * delta_theta + eta / (a * a));
  um = delta_theta - mp::sqrt(delta_theta * delta_theta + eta / (a * a));
  theta_p = mp::acos(-mp::sqrt(up));
  theta_m = mp::acos(mp::sqrt(up));
}

void ForwardRayTracing::reset_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_) {
  ray_status = RayStatus::NORMAL;
  lambda = std::move(lambda_);
  q = std::move(q_);
  eta = q * q;
  nu_r = nu_r_;
  nu_theta = nu_theta_;
  init_radial_potential_roots();
  init_theta_pm();
  reset_variables();
}

void ForwardRayTracing::reset_by_rc_d(const Real &rc, const Real &d, Sign nu_r_, Sign nu_theta_) {
  ray_status = RayStatus::NORMAL;
  Real lambda_c = a + (rc * (2 * square(a) + (-3 + rc) * rc)) / (a - a * rc);
  Real eta_c = -((cube(rc) * (-4 * square(a) + square(-3 + rc) * rc)) /
                 (square(a) * square(-1 + rc)));
  Real qc = mp::sqrt(eta_c);

  Real coeff = mp::sqrt(
      square(-3 + rc) / (square(a) * square(-1 + rc)) + eta_c / mp::pow(rc, 4));

  reset_by_lambda_q(lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff),
                    qc + d * (mp::sqrt(eta_c) / square(rc) / coeff), nu_r_, nu_theta_);
  reset_variables();
}
