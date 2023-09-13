#include <string>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "ForwardRayTracing.h"
#include "Utils.h"

using std::string;

int main(int argc, char *argv[]) {
//    namespace po = boost::program_options;
//
//    po::options_description desc("Allowed options");
//    desc.add_options()
//            ("help", "produce help message")
//            ("a", po::value<std::string>(), "black hole spin");
//
//    po::variables_map vm;
//    po::store(po::parse_command_line(argc, argv, desc), vm);
//    po::notify(vm);
//
//    if (vm.count("help")) {
//        std::cout << desc << "\n";
//        return 1;
//    }

    using Real = Float256;
    using Complex = Complex256;
    ForwardRayTracingParams<Real> params;
    auto forward = ForwardRayTracing<Real, Complex>::get_from_cache();

    std::vector<std::array<std::string, 9>> data;
    std::ifstream ifs(argv[1]);
    std::string line;
    std::array<std::string, 9> row;
    while (std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string cell;
        int i = 0;
        while (std::getline(ss, cell, ',')) {
            row[i++] = cell;
        }
        data.emplace_back(row);
    }

    size_t i_start = 26062;
    size_t i_end = 26100;
    for (size_t i = i_start; i < i_end; i++) {
        const auto &item = data[i];
        params.a = boost::lexical_cast<Real>(item[0]);
        params.r_s = boost::lexical_cast<Real>(item[1]);
        params.theta_s = boost::lexical_cast<Real>(item[2]);
        params.r_o = boost::lexical_cast<Real>(item[3]);
        params.nu_r = Sign::POSITIVE;
        params.nu_theta = Sign::NEGATIVE;

        params.lambda = boost::lexical_cast<Real>(item[4]);
        params.q = sqrt(boost::lexical_cast<Real>(item[5]));
        params.calc_t_f = true;

        forward->calc_ray(params);

        Real delta_t_f = forward->t_f - boost::lexical_cast<Real>(item[6]);
        Real delta_theta_f = forward->theta_f - boost::lexical_cast<Real>(item[7]);
        Real delta_phi_f = forward->phi_f - boost::lexical_cast<Real>(item[8]);

        if (delta_t_f > 1e-3) {
            fmt::println("delta_t_f: {}, delta_theta_f: {}, delta_phi_f: {}", delta_t_f, delta_theta_f, delta_phi_f);
        }
    }
}
