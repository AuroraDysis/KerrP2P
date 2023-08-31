#include "forward.h"

#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

// convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
std::array<Real, 2> lamqv2(Real a, Real rc, Real d) {
  Real lambda_c = a + (rc * (2 * std::pow(a, 2) + (-3 + rc) * rc)) / (a - a * rc);
  Real eta_c = -((std::pow(rc, 3) * (-4 * std::pow(a, 2) + std::pow(-3 + rc, 2) * rc)) /
                 (std::pow(a, 2) * std::pow(-1 + rc, 2)));
  Real qc = std::sqrt(eta_c);

  Real coeff = std::sqrt(std::pow(-3 + rc, 2) / (std::pow(a, 2) * std::pow(-1 + rc, 2)) + eta_c / std::pow(rc, 4));
  std::array<Real, 2> nLambdaQ = {(3 - rc) / (a * (-1 + rc)) / coeff, std::sqrt(eta_c) / std::pow(rc, 2) / coeff};

  // unit normal vector of the critical curve
  return {lambda_c + d * nLambdaQ[0], qc + d * nLambdaQ[1]};
}

//  convert rc, d to lambda, eta
std::array<Real, 2> lametav2(Real a, Real rc, Real d) {
  std::array<Real, 2> lambda_qc = lamqv2(a, rc, d);
  return {lambda_qc[0], lambda_qc[1] * lambda_qc[1]};
}

class RadialAnti {
private:
  Real a, lambda, eta;
  Real r1, r2, r3, r4;
  Real rp, rm;
  Real A, B;
  Real alpha_0, alpha_p, alpha_m, k3;

  Real x3(Real r) const {
    return (A * (r - r1) - B * (r - r2)) / (A * (r - r1) + B * (r - r2));
  }

  Real F3(Real r) const {
    return 1.0 / std::sqrt(A * B) * boost::math::ellint_1(std::acos(x3(r)), k3());
  }

  Real f1(Real alpha, Real curlyPhi, Real j) const {
    Real temp1 = sqrt((alpha * alpha - 1) / (j + (1 - j) * alpha * alpha));
    Real temp2 = temp1 * sqrt(1 - j * sin(curlyPhi) * sin(curlyPhi));
    return temp1 / 2 * log(fabs((temp2 + sin(curlyPhi)) / (temp2 - sin(curlyPhi))));
  }

  int H(Real x) const {
    if (x > 0) {
      return 1;
    } else {
      return 0;
    }
  }

  Real R1(Real alpha, Real curlyPhi, Real j) const {
    return 1.0 / (1.0 - alpha * alpha) *
           (boost::math::ellint_3(alpha * alpha / (alpha * alpha - 1.0), curlyPhi, j) - alpha * f1(alpha, curlyPhi, j));
  }

  Real R2(Real alpha, Real curlyPhi, Real j) const {
    return 1.0 / (alpha * alpha - 1.0) * (boost::math::ellint_1(curlyPhi, j) -
                                          pow(alpha, 2) / (j + (1.0 - j) * alpha * alpha) *
                                          (boost::math::ellint_2(curlyPhi, j) -
                                           (alpha * sin(curlyPhi) * sqrt(1.0 - j * pow(sin(curlyPhi), 2))) /
                                           (1.0 + alpha * cos(curlyPhi)))) +
           1.0 / (j + (1.0 - j) * alpha * alpha) * (2.0 * j - pow(alpha, 2) / (alpha * alpha - 1.0)) *
           R1(alpha, curlyPhi, j);
  }

  Real Pi13(Real r) const {
    // Note: A, B, r1, r2, alpha0, x3, and k3 should be globally defined or passed as parameters
    return (RCONST(2.0) * (r2 - r1) * sqrt(A * B)) / (B * B - A * A) * R1(alpha_0, acos(x3(r)), k3);
  }

  Real Pi23(Real r) const {
    // Note: A, B, r1, r2, alpha0, x3, and k3 should be globally defined or passed as parameters
    return std::pow(((RCONST(2.0) * (r2 - r1) * sqrt(A * B)) / (B * B - A * A)), 2) * R2(alpha_0, acos(x3(r)), k3);
  }

  Real I0(Real r) const {
    return F3(r);
  }

  Real I1(Real r) const {
    return ((B * r2 + A * r1) / (B + A)) * F3(r) + Pi13(r);
  }

  Real I2(Real r) const {
    return pow(((B * r2 + A * r1) / (B + A)), 2) * F3(r) + RCONST(2.0) * ((B * r2 + A * r1) / (B + A)) * Pi13(r) +
           sqrt(A * B) * Pi23(r);
  }

  Real Ip(Real r) const {
    return -1.0 / (B * (rp - r2) + A * (rp - r1)) * ((B + A) * F3(r) +
                                                     (RCONST(2.0) * (r2 - r1) * sqrt(A * B)) /
                                                     (B * (rp - r2) - A * (rp - r1)) *
                                                     R1(alpha_p, acos(x3(r)), k3));
  }

  Real Im(Real r) const {
    return -1.0 / (B * (rm - r2) + A * (rm - r1)) * ((B + A) * F3(r) +
                                                     (RCONST(2.0) * (r2 - r1) * sqrt(A * B)) /
                                                     (B * (rm - r2) - A * (rm - r1)) *
                                                     R1(alpha_m, acos(x3(r)), k3));
  }

public:
  RadialAnti(Real a_, Real lambda_, Real eta_)
      : a(a_), lambda(lambda_), eta(eta_) {
    std::tie(r1, r2, r3, r4) = roots(a, lambda, eta);

    rp = 1 + std::sqrt(1 - a * a);
    rm = 1 - std::sqrt(1 - a * a);

    A = sqrt((r3 - r2) * (r4 - r2));
    B = sqrt((r3 - r1) * (r4 - r1));

    alpha_0 = (B + A) / (B - A);
    alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
    alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));

    k3 = ((A + B) * (A + B) - (r2 - r1) * (r2 - r1)) / (4 * A * B);
  }


  Real I_r(Real r) const {
    return I0(r);
  }

  Real I_phi(Real r) const {
    return (2 * a) / (rp - rm) * ((rp - (a * lambda) / 2) * Ip(r) - (rm - (a * lambda) / 2) * Im(r));
  }

  Real I_t(Real r) const {
    return 4 / (rp - rm) * (rp * (rp - (a * lambda) / 2) * Ip(r) - rm * (rm - (a * lambda) / 2) * Im(r))
           + 4 * I0(r) + 2 * I1(r) + I2(r);
  }
};

std::pair<RayStatus, std::array<Real, 3>> calc_I(Real a, Real lambda, Real eta, Real r_s, Real r_inf, Real nu_r) {
  Real r4 = roots(a, lambda, eta).back();
  Real rp = 1 + std::sqrt(1 - a * a);

  bool radial_turning = (r4 > rp); // Assuming roots returns a real number. Otherwise, check for real numbers.

  RadialAnti radial_anti(a, lambda, eta);
  std::array<std::function<Real(Real)>, 3> anti_ders = {
      [&](Real r) -> Real { return radial_anti.I_r(r); },
      [&](Real r) -> Real { return radial_anti.I_phi(r); },
      [&](Real r) -> Real { return radial_anti.I_t(r); }
  };

  std::array<Real, 3> results = {};
  std::fill(results.begin(), results.end(), std::numeric_limits<Real>::quiet_NaN());

  if (radial_turning) {
    // if there is a radial turning point (i.e. r4 is a real number)
    if (r_s <= r4) {
      return {RayStatus::CONFINED, results};
    } else {
      if (nu_r == 1) {
        // n_ur == +1
        for (int i = 0; i < 3; i++) {
          results[i] = anti_ders[i](r_inf) - anti_ders[i](r_s);
        }
      } else {
        // nu_r == -1
        for (int i = 0; i < 3; i++) {
          results[i] = anti_ders[i](r_inf) + anti_ders[i](r_s) - 2 * anti_ders[i](r4);
        }
      }
    }
  } else {
    // if there is no radial turning point (i.e. r4 is a complex number), then, depending on the value of nu_r = \pm 1
    if (nu_r == -1) {
      // get into the black hole
      return {RayStatus::FALLS_IN, results};
    } else {
      for (int i = 0; i < 3; i++) {
        results[i] = anti_ders[i](r_inf) - anti_ders[i](r_s);
      }
    }
  }

  return {RayStatus::NORMAL, results};
}


class ForwardRayTracing {
private:
  const Real a, rp, rm, r_s, r_inf, theta_s, nu_r, nu_theta;
  // Real lambda, q;
  Real eta, delta_theta, up, um, theta_p, theta_m;
  Real r1, r2, r3, r4;

  // Range of theta: For Type A (eta > 0)
  void theta_pm() {
    delta_theta = RCONST(0.5) * (RCONST(1.0) - (eta + lambda * lambda) / (a * a));
    up = delta_theta + std::sqrt(delta_theta * delta_theta + eta / (a * a));
    um = delta_theta - std::sqrt(delta_theta * delta_theta + eta / (a * a));
    theta_p = std::acos(-std::sqrt(up));
    theta_m = std::acos(std::sqrt(up));
  }

  // find roots
  void roots() {
    Real A = a * a - eta - lambda * lambda;
    Real B = 2 * (eta + (lambda - a) * (lambda - a));
    Real CC = -a * a * eta;
    Real P = -(A * A / 12) - CC;
    Real Q = -(A / 3) * (pow(A / 6, 2) - CC) - B * B / 8;

    Real omega_p, omega_m;

    if (pow(P / 3, 3) + pow(Q / 2, 2) > 0) {
      omega_p = std::pow(-Q / 2 + std::sqrt(pow(P / 3, 3) + pow(Q / 2, 2)), RCONST(1.0) / RCONST(3.0));
      omega_m = std::pow(-Q / 2 - std::sqrt(pow(P / 3, 3) + pow(Q / 2, 2)), RCONST(1.0) / RCONST(3.0));
    } else {
      omega_p = std::cbrt(-Q / 2 + std::sqrt(pow(P / 3, 3) + pow(Q / 2, 2)));
      omega_m = std::cbrt(-Q / 2 - std::sqrt(pow(P / 3, 3) + pow(Q / 2, 2)));
    }

    Real xi0 = omega_p + omega_m - A / 3;
    Real z = std::sqrt(xi0 / 2);

    r1 = -z - std::sqrt(-(A / 2) - z * z + B / (4 * z));
    r2 = -z + std::sqrt(-(A / 2) - z * z + B / (4 * z));
    r3 = z - std::sqrt(-(A / 2) - z * z - B / (4 * z));
    r4 = z + std::sqrt(-(A / 2) - z * z - B / (4 * z));
  }

public:
  // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
  ForwardRayTracing(Real r_s_, Real r_inf_, Real theta_s_, Real nu_r_, Real nu_theta_, Real a_, Real lambda_, Real q_)
      : r_s(r_s_), r_inf(r_inf_), theta_s(theta_s_), nu_r(nu_r_), nu_theta(nu_theta_), a(a_), lambda(lambda_), q(q_),
        rp(1 + std::sqrt(1 - a * a)), rm(1 - std::sqrt(1 - a * a)), eta(q * q) {
    theta_pm();
  }

  Real G_theta(Real theta) {
    return -1.0 / std::sqrt(-um * a * a) * boost::math::ellint_1(std::asin(std::cos(theta) / std::sqrt(up)), up / um);
  }

  Real G_phi(Real theta) {
    return -1.0 / std::sqrt(-um * a * a) *
           boost::math::ellint_3(up, std::asin(std::cos(theta) / std::sqrt(up)), up / um);
  }

  Real G_t(Real theta) {
    return (2.0 * up) / std::sqrt(-um * a * a) * 1.0 / (2.0 * up / um) *
           (boost::math::ellint_2(std::asin(std::cos(theta) / std::sqrt(up)), up / um) - boost::math::ellint_1(theta));
  }

  std::array<Real, 3> calc_G(Real theta_inf, int m) {
    std::array<std::function<Real(Real)>, 3> anti_ders = {
        [&](Real theta) -> Real { return G_theta(theta); },
        [&](Real theta) -> Real { return G_phi(theta); },
        [&](Real theta) -> Real { return G_t(theta); }
    };
    std::array<Real, 3> result = {};
    for (int i = 0; i < 3; i++) {
      result[i] = (m * (anti_ders[i](theta_p) - anti_ders[i](theta_m)) +
                   nu_theta * (pow(-1, m) * anti_ders[i](theta_inf) - anti_ders[i](theta_s)));
    }
    return result;
  }

  RayStatus ray_lambda_q(Real lambda, Real q) {
    if (eta <= 0) {
      return RayStatus::ETA_OUT_OF_RANGE;
    }

    if (theta_s < theta_m || theta_s > theta_p) {
      return RayStatus::THETA_OUT_OF_RANGE;
    }

    roots();

    Real rp = 1 + sqrt(1 - a * a);
    bool radial_turning = r4 > rp;

    if (!radial_turning && nu_r < 0) return RayStatus::FALLS_IN;
    if (radial_turning && r_s <= r4) return RayStatus::CONFINED;

    // Radial integrals
    auto [radial_status, radial_integrals] = calc_I(a, lambda, eta, r_s, r_inf, nu_r);

    // Minor time
    Real tau_o = radial_integrals[0];

    // Final value of theta
    auto scriptG_theta = [&](Real theta) -> Real {
      return RCONST(-1.0) / std::sqrt(-um * a * a) *
             boost::math::ellint_1(std::asin(std::cos(theta) / std::sqrt(up)), up / um);
    };

    Real theta_inf = std::acos(-std::sqrt(up) * nu_theta *
                               boost::math::jacobi_sn(
                                   std::sqrt(-um * a * a) * (tau_o + nu_theta * scriptG_theta(theta_s)), up / um));

    // Angular integrals
    Real m = RCONST(1.0) + std::floor(std::real((tau_o - scriptG_theta(theta_p) + nu_theta * scriptG_theta(theta_s)) /
                                                (scriptG_theta(theta_p) - scriptG_theta(theta_m))));

    std::array<Real, 3> angular_integrals = calc_G(theta_inf, static_cast<int>(m));

    // Number of half-orbits
    Real n_half = tau_o / (scriptG_theta(theta_p) - scriptG_theta(theta_m));

    // Final values of phi and t
    Real phi_inf = radial_integrals[1] + lambda * angular_integrals[1];
    Real t_inf = radial_integrals[2] + a * a * angular_integrals[1];

    std::cout << "Theta_inf: " << theta_inf << ", Phi_inf: " << phi_inf << ", tinf-r_inf: " << (t_inf - r_inf)
              << ", m: "
              << m << ", nhalf: " << n_half << std::endl;
    return RayStatus::NORMAL;
  }
};



//std::array<Real, 2> ob_to_conserve(Real alpha, Real beta, Real theta_o, Real a) {
//  Real lambda, eta;
//
//  lambda = -alpha * std::sin(theta_o);
//  eta = (alpha * alpha - a * a) * std::cos(theta_o) * std::cos(theta_o) + beta * beta;
//
//  return {lambda, eta};
//}
//
//// conserve_to_ob function
//std::array<Real, 2> conserve_to_ob(Real lambda, Real eta, Real theta_o, Real a, Real nu_theta_o) {
//  Real alpha, beta;
//
//  alpha = -lambda / std::sin(theta_o);
//  beta = nu_theta_o * std::sqrt(eta + a * a * std::cos(theta_o) * std::cos(theta_o) -
//                                lambda * lambda / (1 / std::tan(theta_o)) / (1 / std::tan(theta_o)));
//
//  return {alpha, beta};
//}
//

//
//// Potentials
//
//Real R(Real r, Real a, Real lambda, Real eta) {
//  return pow(r * r + a * a - a * lambda, 2) - (r * r - 2 * r + a * a) * (eta + pow(lambda - a, 2));
//}
//
//Real Theta(Real theta, Real a, Real lambda, Real eta) {
//  return eta + a * a * std::cos(theta) * std::cos(theta) - lambda * lambda / std::tan(theta) / std::tan(theta);
//}
//
//std::array<Real, 2> Theta_pm(Real a, Real lambda, Real eta) {
//  Real deltaTheta = 0.5 * (1 - (eta + lambda * lambda) / a / a);
//  Real up = deltaTheta + std::sqrt(deltaTheta * deltaTheta + eta / a / a);
//  Real um = deltaTheta - std::sqrt(deltaTheta * deltaTheta + eta / a / a);
//  Real thetaP = std::acos(-std::sqrt(up));
//  Real thetaM = std::acos(std::sqrt(up));
//  return {thetaP, thetaM};
//}
//
////int main() {
////  // Test the functions with example values.
////  Real a = 1.0, lambda = 2.0, eta = 3.0;
////  std::array<Real, 2> res = roots(a, lambda, eta);
////  for (Real val : res) {
////    std::cout << val << " ";
////  }
////  std::cout << std::endl;
////  return 0;
////}
//

//


//Real FindImage(Real a, Real rs, Real Theta_s, Real Nu_rs, Real Nu_Theta_s,
//               Real Theta_o, Real Phi_o, std::vector<Real> parameters) {
//  Real sol;
//  Real rcsol;
//  Real dsol;
//  Real ray;
//  Real m;
//  Real Lambda;
//  Real Eta;
//  Real r4;
//  Real w0;

//  if (opt == "out") {
//    // sol = FindRoot...
//    std::cout << "{r_c, lgd} = " << sol << std::endl;
//    rcsol = sol; // Assuming a correct implementation of FindRoot
//    dsol = std::pow(10, sol); // For simplicity, this is a placeholder
//  }
//  else if (opt == "in") {
//    // sol = FindRoot...
//    std::cout << "{r_c, lgmd} = " << sol << std::endl;
//    rcsol = sol; // Assuming a correct implementation of FindRoot
//    dsol = -std::pow(10, sol); // For simplicity, this is a placeholder
//  }
//  else {
//    std::cerr << "'opt' should be either 'out' or 'in'!" << std::endl;
//    return -1;  // Error value
//  }
//
//  ray = /* Call to Rayrcd */;
//  m = ray; // Simplified for this demonstration.
//
//  std::cout << "{Theta_inf, Phi_inf, t_inf-r_inf, m, nhalf} = " << ray << std::endl;
//  std::cout << "Mod[Phi_inf, 2 * PI] = " << std::fmod(ray, 2 * M_PI) << std::endl;
//
//  std::tie(Lambda, Eta) = /* Call to lametav2 */;
//  std::cout << "{Lambda, Eta} = " << Lambda << ", " << Eta << std::endl;
//
//  std::vector<Real> results = /* Call to roots */;
//  std::sort(results.begin(), results.end());
//  r4 = results.back();
//
//  w0 = /* Call to Cos or Theta_pm */;
//  std::cout << "{r4, w0} = " << r4 << ", " << w0 << std::endl;
//
//  std::cout << "{Alpha, Beta} = " << /* Call to ConserveToOb */ << std::endl;
//
//  return /* Call to Shape */;
//}

