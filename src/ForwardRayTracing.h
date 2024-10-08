#pragma once

#include "Common.h"
#include "IIntegral2.h"
#include "IIntegral3.h"
#include "GIntegral.h"
#include "ObjectPool.h"

#define CHECK_STATUS                     \
    if (ray_status != RayStatus::NORMAL) \
        return;

template <typename Real>
std::pair<Real, Real> get_rc_range(const Real &a)
{
    Real r_up = 2 * cos(acos(a) * third<Real>());
    Real r_down = 2 * cos(acos(-a) * third<Real>());
    r_up = MY_SQUARE(r_up);
    r_down = MY_SQUARE(r_down);
    return {r_down, r_up};
}

template <typename Real>
struct ForwardRayTracingParams
{
    Real a;
    Real r_s;
    Real theta_s;
    Real r_o;
    Sign nu_r;
    Sign nu_theta;

    Real rc, log_abs_d;
    Sign d_sign;

    Real lambda, q;

    bool calc_t_f = false;
    bool print_args_error = true;

    ForwardRayTracingParams() = default;

    ForwardRayTracingParams(const ForwardRayTracingParams &params)
    {
        a = params.a;
        r_s = params.r_s;
        theta_s = params.theta_s;
        r_o = params.r_o;
        nu_r = params.nu_r;
        nu_theta = params.nu_theta;
        rc = params.rc;
        log_abs_d = params.log_abs_d;
        d_sign = params.d_sign;
        lambda = params.lambda;
        q = params.q;
        calc_t_f = params.calc_t_f;
        print_args_error = params.print_args_error;
    }

    template <typename TH>
    ForwardRayTracingParams<TH> get_high_prec() const
    {
        ForwardRayTracingParams<TH> params;
        params.a = a;
        params.r_s = r_s;
        params.theta_s = theta_s;
        params.r_o = r_o;
        params.nu_r = nu_r;
        params.nu_theta = nu_theta;
        params.rc = rc;
        params.log_abs_d = log_abs_d;
        params.d_sign = d_sign;
        params.lambda = lambda;
        params.q = q;
        params.calc_t_f = calc_t_f;
        return params;
    }

    bool rc_d_to_lambda_q()
    {
        auto [r_down, r_up] = get_rc_range(a);

        if (rc < r_down || rc > r_up)
        {
            if (print_args_error)
                fmt::println("rc out of range: rc = {}, r_down: {}, r_up: {}", rc, r_down, r_up);
            lambda = std::numeric_limits<Real>::quiet_NaN();
            q = std::numeric_limits<Real>::quiet_NaN();
            return false;
        }

        Real lambda_c = a + (rc * (2 * MY_SQUARE(a) + (-3 + rc) * rc)) / (a * (1 - rc));
        Real eta_c = -((MY_CUBE(rc) * (-4 * MY_SQUARE(a) + MY_SQUARE(-3 + rc) * rc)) /
                       (MY_SQUARE(a) * MY_SQUARE(-1 + rc)));
        Real qc = sqrt(eta_c);

        Real coeff = sqrt(
            MY_SQUARE(-3 + rc) / (MY_SQUARE(a) * MY_SQUARE(-1 + rc)) + eta_c / pow(rc, 4));

        Real d = GET_SIGN(d_sign) * pow(static_cast<Real>(10), log_abs_d);
        lambda = lambda_c + d * ((3 - rc) / (a * (-1 + rc)) / coeff);
        q = qc + d * (sqrt(eta_c) / MY_SQUARE(rc) / coeff);

        if (d_sign == Sign::NEGATIVE && q < 0)
        {
            if (print_args_error)
                fmt::println("q out of range, which should be positive when d_sign is NEGATIVE: q = {}", q);
            lambda = std::numeric_limits<Real>::quiet_NaN();
            q = std::numeric_limits<Real>::quiet_NaN();
            return false;
        }
#ifdef PRINT_DEBUG
        fmt::println("rc: {}, log_abs_d: {}", rc, log_abs_d);
        fmt::println("lambda_c: {}, eta_c: {}, qc: {}", lambda_c, eta_c, qc);
        fmt::println("coeff: {}", coeff);
#endif
        return true;
    }
};

template <typename Real, typename Complex>
struct ForwardRayTracingResult
{
    Real a, rp, rm, r_s, theta_s, r_o;
    Real r1, r2, r3, r4;
    Complex r1_c, r2_c, r3_c, r4_c;
    Real t_f, theta_f, phi_f;
    int m;
    Real n_half;
    Real eta, lambda, q;
    Real rc, log_abs_d;
    Sign d_sign;
    RayStatus ray_status;
};

template <typename LReal, typename LComplex, typename Real, typename Complex>
ForwardRayTracingResult<LReal, LComplex> get_low_prec(const ForwardRayTracingResult<Real, Complex> &x)
{
    ForwardRayTracingResult<LReal, LComplex> result;
    result.a = x.a.template convert_to<LReal>();
    result.rp = x.rp.template convert_to<LReal>();
    result.rm = x.rm.template convert_to<LReal>();
    result.r_s = x.r_s.template convert_to<LReal>();
    result.theta_s = x.theta_s.template convert_to<LReal>();
    result.r_o = x.r_o.template convert_to<LReal>();
    result.r1 = x.r1.template convert_to<LReal>();
    result.r2 = x.r2.template convert_to<LReal>();
    result.r3 = x.r3.template convert_to<LReal>();
    result.r4 = x.r4.template convert_to<LReal>();
    result.r1_c = x.r1_c.template convert_to<LComplex>();
    result.r2_c = x.r2_c.template convert_to<LComplex>();
    result.r3_c = x.r3_c.template convert_to<LComplex>();
    result.r4_c = x.r4_c.template convert_to<LComplex>();
    result.t_f = x.t_f.template convert_to<LReal>();
    result.theta_f = x.theta_f.template convert_to<LReal>();
    result.phi_f = x.phi_f.template convert_to<LReal>();
    result.m = x.m;
    result.n_half = x.n_half.template convert_to<LReal>();
    result.eta = x.eta.template convert_to<LReal>();
    result.lambda = x.lambda.template convert_to<LReal>();
    result.q = x.q.template convert_to<LReal>();
    result.rc = x.rc.template convert_to<LReal>();
    result.log_abs_d = x.log_abs_d.template convert_to<LReal>();
    result.d_sign = x.d_sign;
    result.ray_status = x.ray_status;
    return result;
}

template <typename Real, typename Complex>
class ForwardRayTracing
{
    friend class IIntegral2<Real, Complex>;

    friend class IIntegral3<Real, Complex>;

    friend class GIntegral<Real, Complex>;

private:
    inline static ObjectPool<ForwardRayTracing<Real, Complex>> pool;

    void reset_variables()
    {
        ray_status = RayStatus::NORMAL;
        r1 = r2 = r3 = r4 = tau_o = t_f = theta_f = phi_f = n_half = eta = lambda = q = std::numeric_limits<Real>::quiet_NaN();
        r1_c = r2_c = r3_c = r4_c = Complex{std::numeric_limits<Real>::quiet_NaN()};
        m = std::numeric_limits<int>::max();
        std::fill(radial_integrals.begin(), radial_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
        std::fill(angular_integrals.begin(), angular_integrals.end(), std::numeric_limits<Real>::quiet_NaN());
    }

    void init_radial_potential_roots()
    {
        Real AA = a * a - eta - lambda * lambda;
        Real BB = 2 * (eta + (lambda - a) * (lambda - a));
        Real CC = -a * a * eta;
        Real PP = -AA * AA / 12 - CC;
        Real QQ = (-2 * MY_CUBE(AA) - 27 * MY_SQUARE(BB) + 72 * AA * CC) / 216;

        Complex omega_pm;
        Real omega_pm_1 = -QQ * half<Real>();
        Real omega_pm_2 = MY_CUBE(PP) / 27 + MY_SQUARE(QQ) / 4;

        if (omega_pm_2 > 0)
        {
            omega_pm = cbrt(omega_pm_1 + sqrt(omega_pm_2)) +
                       cbrt(omega_pm_1 - sqrt(omega_pm_2));
        }
        else
        {
            Complex omega_pm_2_c{omega_pm_2};
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

        if (real(sqrt_in_1) < 0)
        {
            r12_is_real = false;
            r1 = r2 = std::numeric_limits<Real>::quiet_NaN();
        }
        else
        {
            r12_is_real = true;
            r1 = real(r1_c);
            r2 = real(r2_c);
        }

        if (real(sqrt_in_2) < 0)
        {
            r34_is_real = false;
            r3 = r4 = std::numeric_limits<Real>::quiet_NaN();
        }
        else
        {
            r34_is_real = true;
            r3 = real(r3_c);
            r4 = real(r4_c);
        }

#ifdef PRINT_DEBUG
        fmt::println("AA: {}, BB: {}, CC: {}, PP: {}, QQ: {}", AA, BB, CC, PP, QQ);
        fmt::println("r1: {}, r2: {}, r3: {}, r4: {}", r1, r2, r3, r4);
#endif
    }

    void init_theta_pm()
    {
        delta_theta = half<Real>() * (1 - (eta + MY_SQUARE(lambda)) / MY_SQUARE(a));
        up = delta_theta + sqrt(MY_SQUARE(delta_theta) + eta / MY_SQUARE(a));
        um = delta_theta - sqrt(MY_SQUARE(delta_theta) + eta / MY_SQUARE(a));
        // (28), theta_p \in (pi/2, pi), theta_m \in (0, pi/2)
        theta_p = acos(-sqrt(up));
        theta_m = acos(sqrt(up));
    }

    void reset_by_lambda_q(Real lambda_, Real q_, Sign nu_r_, Sign nu_theta_)
    {
        lambda = std::move(lambda_);
        q = std::move(q_);
        eta = q * q;
        nu_r = nu_r_;
        nu_theta = nu_theta_;
        init_radial_potential_roots();
        init_theta_pm();
#ifdef PRINT_DEBUG
        fmt::println("lambda: {}, q: {}, nu_r: {}, nu_theta: {}", lambda, q, GET_SIGN(nu_r), GET_SIGN(nu_theta));
#endif
    }

    void calcI()
    {
        bool radial_turning = r34_is_real && r4 > rp;

        // if there is a radial turning point (i.e. r4 is a real number)
        if (radial_turning && r_s <= r4)
        {
            ray_status = RayStatus::CONFINED;
            return;
        }

        if (radial_turning && r_s > r4 && nu_r == Sign::POSITIVE)
        {
            I_integral_2->calc(false);
            return;
        }

        if (radial_turning && r_s > r4 && nu_r == Sign::NEGATIVE)
        {
            I_integral_2->calc(true);
            return;
        }

        if (!radial_turning && nu_r == Sign::NEGATIVE)
        {
            ray_status = RayStatus::FALLS_IN;
            return;
        }

        if (!radial_turning && nu_r == Sign::POSITIVE)
        {
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
    RayStatus ray_status;
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

    ForwardRayTracing<Real, Complex>()
    {
        reset_variables();
        I_integral_2 = std::make_shared<IIntegral2<Real, Complex>>(*this);
        I_integral_3 = std::make_shared<IIntegral3<Real, Complex>>(*this);
        G_integral = std::make_shared<GIntegral<Real, Complex>>(*this);
    }

    static std::shared_ptr<ForwardRayTracing<Real, Complex>> get_from_cache()
    {
        return pool.create();
    }

    static void clear_cache()
    {
        pool.clear();
    }

    std::shared_ptr<IIntegral2<Real, Complex>> I_integral_2;
    std::shared_ptr<IIntegral3<Real, Complex>> I_integral_3;
    std::shared_ptr<GIntegral<Real, Complex>> G_integral;

    void calc_ray(const ForwardRayTracingParams<Real> &params)
    {
        reset_variables();

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
        if (isnan(params.lambda) || isnan(params.q) || abs(params.lambda) <= 10000 * ErrorLimit<Real>::Value)
        {
            ray_status = RayStatus::ARGUMENT_ERROR;
            return;
        }

        reset_by_lambda_q(params.lambda, params.q, params.nu_r, params.nu_theta);

        if (eta <= 0)
        {
            ray_status = RayStatus::ETA_OUT_OF_RANGE;
            return;
        }

        if (theta_s < theta_m || theta_s > theta_p)
        {
            ray_status = RayStatus::THETA_OUT_OF_RANGE;
            return;
        }

        if (r34_is_real && (r_s < r4 || r_o < r4))
        {
            ray_status = RayStatus::R_OUT_OF_RANGE;
            return;
        }

        // Radial integrals
        calcI();

        CHECK_STATUS

        tau_o = radial_integrals[0];

        G_integral->calc();

        CHECK_STATUS

        // Final values of phi and t
        phi_f = radial_integrals[1] + lambda * angular_integrals[1];
        if (calc_t_f)
        {
            t_f = radial_integrals[2] + MY_SQUARE(a) * angular_integrals[2];
        }

#ifdef PRINT_DEBUG
        fmt::println("theta_f, phi_f, t_f, m, nhalf: {}, {}, {}, {}, {}", theta_f, phi_f, t_f, m, n_half);
#endif
    }

    ForwardRayTracingResult<Real, Complex> to_result()
    {
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
        result.rc = std::numeric_limits<Real>::quiet_NaN();
        result.log_abs_d = std::numeric_limits<Real>::quiet_NaN();
        return result;
    }
};
