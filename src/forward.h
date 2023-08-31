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

// using Real = double;
using Real = long double;

enum class RayStatus {
  NORMAL,
  FALLS_IN,
  CONFINED,
  ETA_OUT_OF_RANGE, // eta should be positive
  THETA_OUT_OF_RANGE, // theta should be in [theta_m, theta_p]
};

enum class Sign : int {
    POSITIVE = 1,
    NEGATIVE = -1,
};

template<typename E>
constexpr auto to_integral(E e) -> typename std::underlying_type<E>::type
{
    return static_cast<typename std::underlying_type<E>::type>(e);
}

enum class RadialPotentialRootType {
    TYPE_I, // case I: four real roots, two outside horizon, outside critical curve, motion between r4 and infinity with one radial turning
    TYPE_II, // case II: four real roots inside horizon, inside critical curve, motion between r4 and infinity with no radial turning
    TYPE_III, // case III: two real roots inside horizon, inside critical curve, motion between r2 and infinity with no radial turning
    TYPE_IV, // case IV: no real roots, inside critical curve, motion between -z and infinity with no radial turning
};

class ForwardRayTracing {
private:
    const Real a, rp, rm, r_s, r_inf, theta_s;
    const Sign nu_r, nu_theta;
    // need to set before using other method
    Real lambda, q, eta;
    // auto initialized
    Real delta_theta, up, um, theta_p, theta_m;
    Real r1, r2, r3, r4;
    RadialPotentialRootType radial_potential_root_type;
    // radial coeffs
    Real A, B, alpha_0, alpha_p, alpha_m, k3;

    // Range of theta: For Type A (eta > 0)
    void init_theta_pm() {
        delta_theta = half<Real>() * (1 - (eta + lambda * lambda) / (a * a));
        up = delta_theta + std::sqrt(delta_theta * delta_theta + eta / (a * a));
        um = delta_theta - std::sqrt(delta_theta * delta_theta + eta / (a * a));
        theta_p = std::acos(-std::sqrt(up));
        theta_m = std::acos(std::sqrt(up));
    }

    // find roots
    void init_radial_potential_roots() {
        Real AA = a * a - eta - lambda * lambda;
        Real BB = 2 * (eta + (lambda - a) * (lambda - a));
        Real CC = -a * a * eta;
        Real PP = -(AA * AA / 12) - CC;
        Real QQ = (-2 * std::pow(A, 3) - 27 * std::pow(B, 2) + 72 * A * CC) / 216;

        Real omega_p, omega_m;

        if (std::pow(PP * third<Real>(), 3) + std::pow(QQ * half<Real>(), 2) > 0) {
            omega_p = std::pow(
                    -QQ * half<Real>() + std::sqrt(std::pow(PP * third<Real>(), 3) + std::pow(QQ * half<Real>(), 2)),
                    third<Real>());
            omega_m = std::pow(
                    -QQ * half<Real>() - std::sqrt(std::pow(PP * third<Real>(), 3) + std::pow(QQ * half<Real>(), 2)),
                    third<Real>());
        } else {
            omega_p = std::cbrt(
                    -QQ * half<Real>() + std::sqrt(std::pow(PP * third<Real>(), 3) + std::pow(QQ * half<Real>(), 2)));
            omega_m = std::cbrt(
                    -QQ * half<Real>() - std::sqrt(std::pow(PP * third<Real>(), 3) + std::pow(QQ * half<Real>(), 2)));
        }

        Real z = std::sqrt((omega_p + omega_m - AA * third<Real>()) * half<Real>());

        Real sqrt_in_1 = -(AA * half<Real>()) - z * z + BB / (4 * z);
        Real sqrt_in_2 = -(AA * half<Real>()) - z * z - BB / (4 * z);

        std::complex<Real> r1_c = -z - std::sqrt(std::complex<Real>(sqrt_in_1));
        std::complex<Real> r2_c = -z + std::sqrt(std::complex<Real>(sqrt_in_1));
        std::complex<Real> r3_c = z - std::sqrt(std::complex<Real>(sqrt_in_2));
        std::complex<Real> r4_c = z + std::sqrt(std::complex<Real>(sqrt_in_2));

        // determine the radial case
        if (sqrt_in_1 < 0) {
            // case IV: no real roots
            // inside critical curve, motion between -z and infinity with no radial turning
            radial_potential_root_type = RadialPotentialRootType::TYPE_IV;
            r1 = r2 = r3 = r4 = std::numeric_limits<Real>::quiet_NaN();
        } else {
            r1 = std::real(r1_c);
            r2 = std::real(r2_c);

            if (sqrt_in_2 < 0) {
                // case III: two real roots inside horiz.
                // inside critical curve, motion between r2 and infinity with no radial turning
                radial_potential_root_type = RadialPotentialRootType::TYPE_III;
                r3 = r4 = std::numeric_limits<Real>::quiet_NaN();
            } else {
                r3 = std::real(r3_c);
                r4 = std::real(r4_c);

                if (r4 < rp) {
                    // case II: four real roots inside horiz.
                    // inside critical curve, motion between r4 and infinity with no radial turning
                    radial_potential_root_type = RadialPotentialRootType::TYPE_II;
                } else {
                    // std::real(r4) > rp
                    // case I: four real roots, two outside horiz.
                    // outside critical curve, motion between r4 and infinity with one radial turning
                    radial_potential_root_type = RadialPotentialRootType::TYPE_I;
                }
            }
        }
    }

    void init_radial_coeffs() {
        A = std::sqrt((r3 - r2) * (r4 - r2));
        B = std::sqrt((r3 - r1) * (r4 - r1));

        alpha_0 = (B + A) / (B - A);
        alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
        alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));

        k3 = (std::pow(A + B, 2) - std::pow(r2 - r1, 2)) / (4 * A * B);
    }

    Real x3(Real r) const {
        return (A * (r - r1) - B * (r - r2)) / (A * (r - r1) + B * (r - r2));
    }

    Real F3(Real r) const {
        return 1 / std::sqrt(A * B) * boost::math::ellint_1(std::acos(x3(r)), k3);
    }

    Real f1(Real alpha, Real curlyPhi, Real j) const {
        Real temp1 = std::sqrt((alpha * alpha - 1) / (j + (1 - j) * alpha * alpha));
        Real temp2 = temp1 * std::sqrt(1 - j * std::sin(curlyPhi) * std::sin(curlyPhi));
        return temp1 / 2 * std::log(std::abs((temp2 + std::sin(curlyPhi)) / (temp2 - std::sin(curlyPhi))));
    }

    Real R1(Real alpha, Real curlyPhi, Real j) const {
        return 1 / (1 - alpha * alpha) *
               (boost::math::ellint_3(alpha * alpha / (alpha * alpha - 1), curlyPhi, j) -
                alpha * f1(alpha, curlyPhi, j));
    }

    Real R2(Real alpha, Real curlyPhi, Real j) const {
        return 1 / (alpha * alpha - 1) * (boost::math::ellint_1(curlyPhi, j) -
                                          std::pow(alpha, 2) / (j + (1 - j) * alpha * alpha) *
                                          (boost::math::ellint_2(curlyPhi, j) -
                                           (alpha * std::sin(curlyPhi) *
                                            std::sqrt(1 - j * std::pow(std::sin(curlyPhi), 2))) /
                                           (1 + alpha * std::cos(curlyPhi)))) +
               1 / (j + (1 - j) * alpha * alpha) * (2 * j - std::pow(alpha, 2) / (alpha * alpha - 1)) *
               R1(alpha, curlyPhi, j);
    }

    Real Pi13(Real r) const {
        return (2 * (r2 - r1) * std::sqrt(A * B)) / (B * B - A * A) * R1(alpha_0, std::acos(x3(r)), k3);
    }

    Real Pi23(Real r) const {
        return std::pow(((2 * (r2 - r1) * std::sqrt(A * B)) / (B * B - A * A)), 2) *
               R2(alpha_0, std::acos(x3(r)), k3);
    }

    Real I_0(Real r) const {
        return F3(r);
    }

    Real I_1(Real r) const {
        return ((B * r2 + A * r1) / (B + A)) * F3(r) + Pi13(r);
    }

    Real I_2(Real r) const {
        return std::pow(((B * r2 + A * r1) / (B + A)), 2) * F3(r) +
               2 * ((B * r2 + A * r1) / (B + A)) * Pi13(r) +
               std::sqrt(A * B) * Pi23(r);
    }

    Real I_p(Real r) const {
        return -1 / (B * (rp - r2) + A * (rp - r1)) * ((B + A) * F3(r) +
                                                       (2 * (r2 - r1) * std::sqrt(A * B)) /
                                                       (B * (rp - r2) - A * (rp - r1)) *
                                                       R1(alpha_p, std::acos(x3(r)), k3));
    }

    Real I_m(Real r) const {
        return -1 / (B * (rm - r2) + A * (rm - r1)) * ((B + A) * F3(r) +
                                                       (2 * (r2 - r1) * std::sqrt(A * B)) /
                                                       (B * (rm - r2) - A * (rm - r1)) *
                                                       R1(alpha_m, std::acos(x3(r)), k3));
    }

public:
    // 输入参数lambda, q，输出光线到无穷远处的theta、phi、传播时间、角向转折次数m、角向"半轨道"数
    ForwardRayTracing(Real a_, Real r_s_, Real theta_s_, Real r_inf_, Sign nu_r_, Sign nu_theta_)
            : a(a_), r_s(r_s_), theta_s(theta_s_), r_inf(r_inf_), nu_r(nu_r_), nu_theta(nu_theta_),
              rp(1 + std::sqrt(1 - a_ * a_)), rm(1 - std::sqrt(1 - a_ * a_)) {
    }

    void init_by_lambda_q(Real lambda_, Real q_) {
        lambda = lambda_;
        q = q_;
        eta = q_ * q_;
        init_radial_potential_roots();
        init_theta_pm();
        init_radial_coeffs();
    }

    // convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
    void init_by_rc_d(Real rc, Real d) {
        Real lambda_c = a + (rc * (2 * std::pow(a, 2) + (-3 + rc) * rc)) / (a - a * rc);
        Real eta_c = -((std::pow(rc, 3) * (-4 * std::pow(a, 2) + std::pow(-3 + rc, 2) * rc)) /
                       (std::pow(a, 2) * std::pow(-1 + rc, 2)));
        Real qc = std::sqrt(eta_c);

        Real coeff = std::sqrt(
                std::pow(-3 + rc, 2) / (std::pow(a, 2) * std::pow(-1 + rc, 2)) + eta_c / std::pow(rc, 4));

        // unit normal vector of the critical curve
        init_by_lambda_q(lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff),
                         qc + d * (std::sqrt(eta_c) / std::pow(rc, 2) / coeff));
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

    std::pair<RayStatus, std::array<Real, 3>> calc_I() const {
        bool radial_turning = r4 > rp; // Assuming roots returns a real number. Otherwise, check for real numbers.

        std::array<std::function<Real(Real)>, 3> anti_ders = {
                [&](Real r) -> Real { return I_r(r); },
                [&](Real r) -> Real { return I_phi(r); },
                [&](Real r) -> Real { return I_t(r); }
        };

        std::array<Real, 3> results = {};
        std::fill(results.begin(), results.end(), std::numeric_limits<Real>::quiet_NaN());

        if (radial_turning) {
            // if there is a radial turning point (i.e. r4 is a real number)
            if (r_s <= r4) {
                return {RayStatus::CONFINED, results};
            } else {
                if (nu_r == Sign::POSITIVE) {
                    // nu_r == +1
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
            if (nu_r == Sign::NEGATIVE) {
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

    Real G_theta(Real theta) const {
        return -1 / std::sqrt(-um * a * a) *
               boost::math::ellint_1(std::asin(std::cos(theta) / std::sqrt(up)), up / um);
    }

    Real G_phi(Real theta) const {
        return -1 / std::sqrt(-um * a * a) *
               boost::math::ellint_3(up, std::asin(std::cos(theta) / std::sqrt(up)), up / um);
    }

    Real G_t(Real theta) const {
        return (2 * up) / std::sqrt(-um * a * a) / (2 * up / um) *
               (boost::math::ellint_2(std::asin(std::cos(theta) / std::sqrt(up)), up / um) -
                boost::math::ellint_1(theta));
    }

    std::array<Real, 3> calc_G(Real theta_inf, int m) const {
        std::array<Real, 3> result = {};
        result[0] = m * (G_theta(theta_p) - G_theta(theta_m)) +
                    to_integral(nu_theta) * (pow(-1, m) * G_theta(theta_inf) - G_theta(theta_s));
        result[1] = m * (G_phi(theta_p) - G_phi(theta_m)) +
                    to_integral(nu_theta) * (pow(-1, m) * G_phi(theta_inf) - G_phi(theta_s));
        result[2] = m * (G_t(theta_p) - G_t(theta_m)) +
                    to_integral(nu_theta) * (pow(-1, m) * G_t(theta_inf) - G_t(theta_s));
        return result;
    }

    RayStatus ray_lambda_q() {
        if (eta <= 0) {
            return RayStatus::ETA_OUT_OF_RANGE;
        }

        if (theta_s < theta_m || theta_s > theta_p) {
            return RayStatus::THETA_OUT_OF_RANGE;
        }

        bool radial_turning = r4 > rp;

        if (!radial_turning && nu_r == Sign::NEGATIVE) return RayStatus::FALLS_IN;
        if (radial_turning && r_s <= r4) return RayStatus::CONFINED;

        // Radial integrals
        auto [radial_status, radial_integrals] = calc_I();

        Real G_theta_theta_s = G_theta(theta_s);
        Real G_theta_theta_p = G_theta(theta_p);
        Real G_theta_theta_m = G_theta(theta_m);

        // Minor time
        Real tau_o = radial_integrals[0];

        Real theta_inf = std::acos(-std::sqrt(up) * to_integral(nu_theta) *
                                   boost::math::jacobi_sn(
                                           std::sqrt(-um * a * a) * (tau_o + to_integral(nu_theta) * G_theta_theta_s),
                                           up / um));

        // Angular integrals
        Real m_Real = 1 + std::floor(std::real((tau_o - G_theta_theta_p + to_integral(nu_theta) * G_theta_theta_s) /
                                               (G_theta_theta_p - G_theta_theta_m)));

        using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
                boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;
        // floor
        int m = RealToInt::convert(m_Real);

        std::array<Real, 3> angular_integrals = calc_G(theta_inf, m);

        // Number of half-orbits
        Real n_half = tau_o / (G_theta_theta_p - G_theta_theta_m);

        // Final values of phi and t
        Real phi_inf = radial_integrals[1] + lambda * angular_integrals[1];
        Real t_inf = radial_integrals[2] + a * a * angular_integrals[1];

        std::cout << "Theta_inf: " << theta_inf << ", Phi_inf: " << phi_inf << ", tinf-r_inf: " << (t_inf - r_inf)
                  << ", m: "
                  << m << ", nhalf: " << n_half << std::endl;
        return RayStatus::NORMAL;
    }
};
