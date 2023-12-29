#pragma once

#include "Common.h"
#include "Integral.h"

#include <boost/math/special_functions/ellint_rc.hpp>

// Radial Antiderivatives for case (3)
template<typename Real, typename Complex>
class IIntegral3 : public Integral<Real, Complex> {
private:
    Real ellint_phi, ellint_cos_phi, ellint_sin_phi, ellint_sin_phi2;
    Real r34_re, r34_im;
    Real A, B, alpha_p, alpha_m, ellint_k, ellint_m, alpha2, ellint1_phi;
    Real alpha_0, R1_alpha_0, R2_alpha_0, Pi_13, Pi_23, I1, I2;
    Real ellint3_n, ellint_3_tmp, f1, p1, ellint3_n1;
    Real F3, R1_alpha_p, R1_alpha_m, Ip, Im;

    std::array<Real, 3> integral_rs;
    std::array<Real, 3> integral_ro;

    void R1(Real &res, const Real &alpha) {
        // alpha2 > 1
        alpha2 = MY_SQUARE(alpha);
        CHECK_VAR(alpha2, alpha2 > 1);

        // alpha2 / (alpha2 - 1) > 1
        ellint3_n = alpha2 / (alpha2 - 1);
        ellint3_n1 = ellint_m / ellint3_n;

        using boost::math::ellint_rc;
        using boost::math::ellint_rj;

        using boost::math::constants::half_pi;
        using boost::math::constants::two_thirds;

        // Cauchy principal value: https://dlmf.nist.gov/19.25.E16

        // Symmetric Integrals: https://dlmf.nist.gov/19.20.E6
        ellint_3_tmp = -third<Real>() * ellint3_n1 * ellint_sin_phi2 * ellint_sin_phi *
            ellint_rj(1 - ellint_sin_phi2, 1 - ellint_m * ellint_sin_phi2, 1, 1 - ellint3_n1 * ellint_sin_phi2);
        // Symmetric Integrals: https://dlmf.nist.gov/19.20.E1
        ellint_3_tmp += sqrt((1 - ellint_sin_phi2) * (1 - ellint_m * ellint_sin_phi2) / ((ellint3_n - 1) * (1 - ellint3_n1))) *
            ellint_rc(ellint_sin_phi2 * (ellint3_n - 1) * (1 - ellint3_n1),
                (ellint3_n * ellint_sin_phi2 - 1) * (1 - ellint3_n1 * ellint_sin_phi2));
        if (ellint_phi >= half_pi<Real>()) {
            // Cauchy principal value: https://dlmf.nist.gov/19.25.E4
            ellint_3_tmp += two_thirds<Real>() * ellint3_n1 * ellint_rj(0, 1 - ellint_m, 1, 1 - ellint3_n1);
            ellint_3_tmp = -ellint_3_tmp;
        }

        // p_1 > 0 (B65)
        p1 = sqrt((alpha2 - 1) / (ellint_m + (1 - ellint_m) * alpha2));
        f1 = half<Real>() * p1 * log(abs((ellint_sin_phi + p1 * sqrt(1 - ellint_m * ellint_sin_phi2)) /
                                         (-ellint_sin_phi + p1 * sqrt(1 - ellint_m * ellint_sin_phi2))));
#ifdef PRINT_DEBUG
        fmt::println("R1 - alpha: {}, ellint_phi: {}, ellint_m: {}, ellint3_n: {}, ellint3_n1: {}, alpha2: {}, f1: {}, ellint_sin_phi2: {}, ellint_3_tmp: {}", alpha, ellint_phi, ellint_m, ellint3_n, ellint3_n1, alpha2, f1, ellint_sin_phi2, ellint_3_tmp);
#endif

        res = (-ellint_3_tmp + alpha * f1) / (alpha2 - 1);
    }

public:
    explicit IIntegral3(ForwardRayTracing<Real, Complex> &data_) : Integral<Real, Complex>(data_,
                                                                                           TypeName<IIntegral3<Real, Complex>>::Get()) {
    }

    void pre_calc() {
        const Real &rp = this->data.rp;
        const Real &rm = this->data.rm;
        // two real roots, both inside horizon, r_1 < r_2 < r_- < r_+ and r_3 = conj(r_4)
        const Real &r1 = this->data.r1;
        const Real &r2 = this->data.r2;
        const Complex &r3 = this->data.r3_c;

        r34_re = real(r3);
        r34_im = imag(r3);

        // radial coeffs
        A = sqrt(MY_SQUARE(r34_im) + MY_SQUARE(r34_re - r2));
        B = sqrt(MY_SQUARE(r34_im) + MY_SQUARE(r34_re - r1));

#ifdef PRINT_DEBUG
        fmt::println("I3 - A: {}, B: {}", A, B);
#endif

        // k3 \in (0, 1)
        ellint_m = ((A + B + r1 - r2) * (A + B - r1 + r2)) / (4 * A * B);
        CHECK_VAR(ellint_m, ellint_m > 0 && ellint_m < 1);
        ellint_k = sqrt(ellint_m);

        alpha_0 = (B + A) / (B - A);
        alpha_p = (B * (rp - r2) + A * (rp - r1)) / (B * (rp - r2) - A * (rp - r1));
        alpha_m = (B * (rm - r2) + A * (rm - r1)) / (B * (rm - r2) - A * (rm - r1));

#ifdef PRINT_DEBUG
        fmt::println("I3 - ellint_k: {}, alpha_p: {}, alpha_m: {}, alpha_0: {}", ellint_k, alpha_p, alpha_m, alpha_0);
#endif
    }

    void calc_x(std::array<Real, 3> &integral, const Real &r) {
        const Real &a = this->data.a;
        const Real &lambda = this->data.lambda;
        const Real &rp = this->data.rp;
        const Real &rm = this->data.rm;
        const Real &r1 = this->data.r1;
        const Real &r2 = this->data.r2;
        // x_3 \in (- 1/alpha_0, 1) \in (-1, 1)
        if (isinf(r)) {
            ellint_cos_phi = (A - B) / (A + B);
        } else {
            ellint_cos_phi = -1 + (2 * A * (r - r1)) / (A * (r - r1) + B * (r - r2));
        }

        CHECK_VAR(ellint_cos_phi, ellint_cos_phi >= -1 && ellint_cos_phi <= 1);
        ellint_sin_phi2 = 1 - MY_SQUARE(ellint_cos_phi);

        ellint_sin_phi = sqrt(ellint_sin_phi2);
        ellint_phi = acos(ellint_cos_phi);

        R1(R1_alpha_p, alpha_p);
        CHECK_DATA_STATUS

        R1(R1_alpha_m, alpha_m);
        CHECK_DATA_STATUS

        ellint1_phi = ellint_1(ellint_k, ellint_phi);
        F3 = ellint1_phi / sqrt(A * B);
        Ip = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_p * (-r1 + r2)) /
                               (A * (r1 - rp) + B * (-r2 + rp))) /
               (-(A * r1) - B * r2 + (A + B) * rp));
        Im = -(((A + B) * F3 + (2 * sqrt(A * B) * R1_alpha_m * (-r1 + r2)) /
                               (A * (r1 - rm) + B * (-r2 + rm))) /
               (-(A * r1) - B * r2 + (A + B) * rm));
        integral[0] = F3;
        integral[1] = (a * (Im * (-(a * lambda) + 2 * rm) + Ip * (a * lambda - 2 * rp))) / (rm - rp);
        if (this->data.calc_t_f && !isinf(this->data.r_o)) {
            alpha2 = MY_SQUARE(alpha_0);

            R1(R1_alpha_0, alpha_0);
            CHECK_DATA_STATUS

            R2_alpha_0 = ((-1 + alpha2) * ellint_m * (1 + alpha_0 * ellint_cos_phi) *
                          (ellint1_phi - 2 * R1_alpha_0) +
                          alpha2 * (1 + alpha_0 * ellint_cos_phi) *
                          (ellint_2(ellint_k, ellint_phi) - ellint1_phi + R1_alpha_0) -
                          MY_CUBE(alpha_0) * ellint_sin_phi * sqrt(1 - ellint_m * ellint_sin_phi2)) /
                         ((-1 + alpha2) * (alpha2 * (-1 + ellint_m) - ellint_m) *
                          (1 + alpha_0 * ellint_cos_phi));
            Pi_13 = (2 * sqrt(A * B) * R1_alpha_0 * (r1 - r2)) / ((A - B) * (A + B));
            Pi_23 = (4 * (A * B) * R2_alpha_0 * MY_SQUARE(r1 - r2)) / MY_SQUARE((A - B) * (A + B));
            I1 = (A * r1 + B * r2) * F3 / (A + B) + Pi_13;
            I2 =
                    sqrt(A * B) * Pi_23 +
                    ((A * r1 + B * r2) * (2 * (A + B) * Pi_13 + F3 * (A * r1 + B * r2))) / MY_SQUARE(A + B);
            integral[2] =
                    4 * F3 + I2 +
                    (2 * Im * rm * (-(a * lambda) + 2 * rm) + 2 * Ip * (a * lambda - 2 * rp) * rp) / (rm - rp) +
                    2 * I1;

#ifdef PRINT_DEBUG
            fmt::println("I3 - R1_alpha_0: {}, R2_alpha_0: {}, Pi_13: {}, Pi_23: {}, I1: {}, I2: {}", R1_alpha_0,
                         R2_alpha_0, Pi_13, Pi_23, I1, I2);
#endif
        } else {
            integral[2] = std::numeric_limits<Real>::quiet_NaN();
        }

#ifdef PRINT_DEBUG
        fmt::println("I3 - F3: {}, R1_alpha_p: {}, R1_alpha_m: {}, Ip: {}, Im: {}", F3, R1_alpha_p, R1_alpha_m, Ip, Im);
        fmt::println("I3 - I_r: {}, I_phi: {}, I_t: {}", integral[0], integral[1], integral[2]);
#endif
    }

    void calc(bool is_plus) {
        pre_calc();

        CHECK_DATA_STATUS

        auto &r_s = this->data.r_s;
        calc_x(integral_rs, r_s);

        CHECK_DATA_STATUS

        auto &r_o = this->data.r_o;
        calc_x(integral_ro, r_o);

        CHECK_DATA_STATUS

        auto &radial_integrals = this->data.radial_integrals;
        for (int i = 0; i < 3; ++i) {
            radial_integrals[i] = is_plus ? integral_ro[i] + integral_rs[i] : integral_ro[i] - integral_rs[i];
        }
#ifdef PRINT_DEBUG
        fmt::println("I3: {}, {}, {}", radial_integrals[0], radial_integrals[1], radial_integrals[2]);
#endif
    }
};

template<typename Real, typename Complex>
struct TypeName<IIntegral3<Real, Complex>> {
    static std::string Get() {
        return fmt::format("IIntegral3<{}, {}>", TypeName<Real>::Get(), TypeName<Complex>::Get());
    }
};
