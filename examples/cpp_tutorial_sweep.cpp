#include <iostream>

#include "ForwardRayTracing.h"
#include "Utils.h"

using std::string;

int main(int argc, char *argv[]) {
    using Real = double;
    using Complex = std::complex<double>;
    ForwardRayTracingParams<Real> params;

    const auto &pi = boost::math::constants::pi<Real>();
    params.a = boost::lexical_cast<Real>("0.8");
    params.r_s = boost::lexical_cast<Real>("10.0");
    params.theta_s = boost::lexical_cast<Real>("85") * pi / 180;
    params.r_o = 1000;
    params.nu_r = Sign::NEGATIVE;
    params.nu_theta = Sign::NEGATIVE;
    params.d_sign = Sign::POSITIVE;

    auto [rc_down, rc_up] = get_rc_range(params.a);
    rc_down += 0.05;
    rc_up -= 0.05;

    std::cout << "rc_down: " << rc_down << ", rc_up: " << rc_up << std::endl;

    std::vector<Real> rc_list(1000);
    std::vector<Real> lgd_list(2000);
    for (int i = 0; i < rc_list.size(); i++) {
        rc_list[i] = rc_down + (rc_up - rc_down) * i / (rc_list.size() - 1.);
    }
    for (int i = 0; i < lgd_list.size(); i++) {
        lgd_list[i] = -10 + 12 * i / (lgd_list.size() - 1.);
    }

    double theta_o = 17 * pi / 180;
    double phi_o = pi / 4;
    int cut_off = 50;
    double tol = 1e-6;
    auto data = ForwardRayTracingUtils<Real, Complex>::sweep_rc_d(params, theta_o, phi_o, rc_list, lgd_list, cut_off,
                                                                  tol);

    // print data.results
    std::cout << "data.results: " << data.results.size() << std::endl;
    for (auto &item: data.results) {
        std::cout << item.rc << ", " << item.log_abs_d << std::endl;
    }
}
