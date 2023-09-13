#pragma once

#include "Common.h"
#include "Integral.h"

// radial anti-derivatives for case (2)
template<typename Real, typename Complex>
class IIntegral2 : public Integral<Real, Complex> {
private:
    Real ellint_sin_phi_rs2, ellint_sin_phi_ro2, ellint_sin_phi, ellint_sin_phi3, ellint_cos_phi2, ellint_c, ellint_k, ellint_phi;
    Real E2_coeff, F2_coeff, Pi_p2_coeff, Pi_m2_coeff, Pi_p2_ellint_n, Pi_m2_ellint_n;

    Real F2, Pi_p2, Pi_m2, Ip, Im;
    Real E2, Pi_12, I1, I2, Pi_12_ellint_n;

    std::array<Real, 3> integral_rs;
    std::array<Real, 3> integral_ro;
public:
    explicit IIntegral2(ForwardRayTracing<Real, Complex> &data_) : Integral<Real, Complex>(data_,
                                                                                           TypeName<IIntegral2<Real, Complex>>::Get()) {
    }

    void pre_calc() {
        const Real &rp = this->data.rp;
        const Real &rm = this->data.rm;

        const Real &r1 = this->data.r1;
        const Real &r2 = this->data.r2;
        const Real &r3 = this->data.r3;
        const Real &r4 = this->data.r4;
        const Real &r_s = this->data.r_s;
        const Real &r_o = this->data.r_o;

        ellint_sin_phi_rs2 = ((r1 - r3) * (r_s - r4)) / ((r_s - r3) * (r1 - r4));
        if (isinf(r_o)) {
            ellint_sin_phi_ro2 = (r1 - r3) / (r1 - r4);
        } else {
            ellint_sin_phi_ro2 = ((r1 - r3) * (r_o - r4)) / ((r_o - r3) * (r1 - r4));
        }
        CHECK_VAR(ellint_sin_phi_rs2, ellint_sin_phi_rs2 >= 0 && ellint_sin_phi_rs2 <= 1);
        CHECK_VAR(ellint_sin_phi_ro2, ellint_sin_phi_ro2 >= 0 && ellint_sin_phi_ro2 <= 1);

        ellint_k = sqrt(((-r2 + r3) * (-r1 + r4)) / ((r1 - r3) * (r2 - r4)));

        E2_coeff = sqrt((r1 - r3) * (r2 - r4));
        F2_coeff = 2 / E2_coeff;
        Pi_p2_coeff = 2 * (-r3 + r4) / (E2_coeff * (-r3 + rp) * (-r4 + rp));
        Pi_m2_coeff = 2 * (-r3 + r4) / (E2_coeff * (-r3 + rm) * (-r4 + rm));
        Pi_p2_ellint_n = ((r1 - r4) * (r3 - rp)) / ((r1 - r3) * (r4 - rp));
        Pi_m2_ellint_n = ((r1 - r4) * (r3 - rm)) / ((r1 - r3) * (r4 - rm));
        Pi_12_ellint_n = (r1 - r4) / (r1 - r3);
    }

    void calc_x(std::array<Real, 3> &integral, const Real &ellint_sin_phi2, const Real &r) {
        const Real &a = this->data.a;
        const Real &lambda = this->data.lambda;
        const Real &rp = this->data.rp;
        const Real &rm = this->data.rm;
        const Real &r3 = this->data.r3;

        ellint_cos_phi2 = 1 - ellint_sin_phi2;
        ellint_sin_phi = sqrt(ellint_sin_phi2);
        ellint_sin_phi3 = ellint_sin_phi2 * ellint_sin_phi;

        ellint_phi = asin(ellint_sin_phi);
        F2 = F2_coeff * ellint_1(ellint_k, ellint_phi);
        Pi_p2 = Pi_p2_coeff * ellint_3(ellint_k, Pi_p2_ellint_n, ellint_phi);
        Pi_m2 = Pi_m2_coeff * ellint_3(ellint_k, Pi_m2_ellint_n, ellint_phi);

        Ip = F2 / (r3 - rp) - Pi_p2;
        Im = F2 / (r3 - rm) - Pi_m2;

#ifdef PRINT_DEBUG
        fmt::println("I2 - ellint_sin_phi2: {}, ellint_k: {}", ellint_sin_phi2, ellint_k);
        fmt::println("I2 - F2: {}, E2: {}, Pi_p2: {}, Pi_m2: {}, Ip: {}, Im: {}", F2, E2, Pi_p2, Pi_m2, Ip, Im);
#endif

        // I_r
        integral[0] = F2;
        // I_phi
        integral[1] = (a * (-2 * rp * Ip + a * Ip * lambda + (2 * rm - a * lambda) * Im)) / (rm - rp);
        if (this->data.calc_t_f && !isinf(this->data.r_o)) {
            // I_t
            const Real &r1 = this->data.r1;
            const Real &r2 = this->data.r2;
            const Real &r4 = this->data.r4;
            const Real &eta = this->data.eta;

            using boost::math::ellint_rd;
            E2 = E2_coeff * boost::math::ellint_2(ellint_k, ellint_phi);
            Pi_12 = F2_coeff * boost::math::ellint_3(ellint_k, Pi_12_ellint_n, ellint_phi);
            I1 = r3 * F2 + (r4 - r3) * Pi_12;
            I2 = -E2 + sqrt(-((eta + MY_SQUARE(a - lambda)) * (MY_SQUARE(a) + (-2 + r) * r)) +
                            MY_SQUARE(MY_SQUARE(a) - a * lambda + MY_SQUARE(r))) / (r - r3) -
                 (F2 * (r2 * r3 + r1 * r4)) * half<Real>();
            integral[2] = 4 * F2 + 2 * I1 + I2 +
                          (-2 * a * Im * lambda * rm + 4 * Im * MY_SQUARE(rm) +
                           2 * Ip * (a * lambda - 2 * rp) * rp) /
                          (rm - rp);
#ifdef PRINT_DEBUG
            fmt::println("I2 - E2: {}, Pi_12: {}, I1: {}, I2: {}", E2, Pi_12, I1, I2);
#endif
        } else {
            integral[2] = std::numeric_limits<Real>::quiet_NaN();
        }
    }

    void calc(bool is_plus) {
        pre_calc();

        CHECK_DATA_STATUS

        const Real &r_s = this->data.r_s;
        calc_x(integral_rs, ellint_sin_phi_rs2, r_s);

        CHECK_DATA_STATUS

        const Real &r_o = this->data.r_o;
        calc_x(integral_ro, ellint_sin_phi_ro2, r_o);

        CHECK_DATA_STATUS

        auto &radial_integrals = this->data.radial_integrals;
        for (int i = 0; i < 3; ++i) {
            radial_integrals[i] = is_plus ? integral_ro[i] + integral_rs[i] : integral_ro[i] - integral_rs[i];
        }
#ifdef PRINT_DEBUG
        fmt::println("I2: {}, {}, {}", radial_integrals[0], radial_integrals[1], radial_integrals[2]);
#endif
    }
};

template<typename Real, typename Complex>
struct TypeName<IIntegral2<Real, Complex>> {
    static std::string Get() {
        return fmt::format("IIntegral2<{}, {}>", TypeName<Real>::Get(), TypeName<Complex>::Get());
    }
};
