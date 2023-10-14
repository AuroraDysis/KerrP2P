#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "ForwardRayTracing.h"
#include "Utils.h"

using std::string;

int main(int argc, char *argv[]) {
  using Real = double;
  using Complex = std::complex<Real>;
  ForwardRayTracingParams<Real> params;

  auto pi = boost::math::constants::pi<Real>();
  params.a = 0.8;
  params.r_s = 10;
  params.theta_s = 85 * pi / 180;
  auto theta_o = 17 * pi / 180;
  auto phi_o = pi / 4;
  params.r_o = 1000;
  params.nu_r = Sign::POSITIVE;
  params.nu_theta = Sign::NEGATIVE;
  params.rc = 1.8117208808167675;
  params.log_abs_d = 0.34917458729364625;
  params.log_abs_d_sign = Sign::NEGATIVE;

  params.rc_d_to_lambda_q();

  auto forward = ForwardRayTracing<Real, Complex>::get_from_cache();
  forward->calc_ray(params);
  fmt::println("ray status: {}", ray_status_to_str(forward->ray_status));
  fmt::println("theta_f: {}, phi_f: {}", forward->theta_f, forward->phi_f);
}

//int main(int argc, char *argv[]) {
////    namespace po = boost::program_options;
////
////    po::options_description desc("Allowed options");
////    desc.add_options()
////            ("help", "produce help message")
////            ("a", po::value<std::string>(), "black hole spin");
////
////    po::variables_map vm;
////    po::store(po::parse_command_line(argc, argv, desc), vm);
////    po::notify(vm);
////
////    if (vm.count("help")) {
////        std::cout << desc << "\n";
////        return 1;
////    }
//
//    using Real = Float256;
//    using Complex = Complex256;
//    ForwardRayTracingParams<Real> params;
//    auto forward = ForwardRayTracing<Real, Complex>::get_from_cache();
//
//    std::vector<std::array<std::string, 9>> data;
//    std::ifstream ifs(argv[1]);
//    std::string line;
//    while (std::getline(ifs, line)) {
//        std::stringstream ss(line);
//        std::array<std::string, 9> row;
//        std::string cell;
//        int i = 0;
//        while (std::getline(ss, cell, ',')) {
//            row[i++] = cell;
//        }
//        data.push_back(row);
//    }
//    fmt::println("data size: {}", data.size());
//
//    size_t i_start = 26062;
//    size_t i_end = 26100;
//    for (size_t i = i_start; i < i_end; i++) {
//        fmt::println("i: {}", i);
//        const auto &item = data[i];
//        params.a = boost::lexical_cast<Real>(item[0]);
//        params.r_s = boost::lexical_cast<Real>(item[1]);
//        params.theta_s = boost::lexical_cast<Real>(item[2]);
//        params.r_o = boost::lexical_cast<Real>(item[3]);
//        params.nu_r = Sign::POSITIVE;
//        params.nu_theta = Sign::NEGATIVE;
//
//        params.lambda = boost::lexical_cast<Real>(item[4]);
//        params.q = sqrt(boost::lexical_cast<Real>(item[5]));
//        fmt::println("lambda = {}, q = {}", item[4], item[5]);
//        fmt::println("a: {}, r_s: {}, theta_s: {}, r_o: {}, lambda: {}, q: {}", params.a, params.r_s, params.theta_s,
//                     params.r_o, params.lambda, params.q);
//        params.calc_t_f = true;
//
//        forward->calc_ray(params);
//
//        Real delta_t_f = forward->t_f - boost::lexical_cast<Real>(item[6]);
//        Real delta_theta_f = forward->theta_f - boost::lexical_cast<Real>(item[7]);
//        Real delta_phi_f = forward->phi_f - boost::lexical_cast<Real>(item[8]);
//
//        if (delta_t_f > 1e-3) {
//            fmt::println("delta_t_f: {}, delta_theta_f: {}, delta_phi_f: {}", delta_t_f, delta_theta_f, delta_phi_f);
//        }
//    }
//}
