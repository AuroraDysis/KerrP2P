#include <tuple>

#include "TestData.h"
#include <oneapi/tbb.h>

using std::string;

template<typename Real, typename Complex>
void test_case(std::vector<std::array<std::string, 9>> &test_data, Sign nu_r, Sign nu_theta) {
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    if (test_data.empty()) {
        fmt::println(std::cerr,
                     "No test data found, nu_r = {}, nu_theta = {}. Please set the data path with -d or --data_path",
                     GET_SIGN(nu_r), GET_SIGN(nu_theta));
        return;
    }

    Vector t_f_vec = Vector::Zero(test_data.size());
    Vector theta_f_vec = Vector::Zero(test_data.size());
    Vector phi_f_vec = Vector::Zero(test_data.size());

    std::mutex mtx;
    std::vector<size_t> error_indices;

    oneapi::tbb::parallel_for(tbb::blocked_range<size_t>(0u, test_data.size()), [&](const auto &range) {
        ForwardRayTracingParams<Real> params;
        auto forward = ForwardRayTracing<Real, Complex>::get_from_cache();
        for (size_t i = range.begin(); i < range.end(); ++i) {
            const auto &item = test_data[i];
            params.a = boost::lexical_cast<Real>(item[0]);
            params.r_s = boost::lexical_cast<Real>(item[1]);
            params.theta_s = boost::lexical_cast<Real>(item[2]);
            params.r_o = boost::lexical_cast<Real>(item[3]);
            params.nu_r = nu_r;
            params.nu_theta = nu_theta;

            params.lambda = boost::lexical_cast<Real>(item[4]);
            params.q = sqrt(boost::lexical_cast<Real>(item[5]));
            params.calc_t_f = true;

            try {
                forward->calc_ray(params);
                // CHECK(forward->ray_status == RayStatus::NORMAL);

                if (forward->ray_status != RayStatus::NORMAL) {
                    fmt::println(std::cerr, "[{}, {}, {}] Ray status: {}", GET_SIGN(nu_r), GET_SIGN(nu_theta), i,
                                 ray_status_to_str(forward->ray_status));
                    {
                        std::lock_guard<std::mutex> lock(mtx);
                        error_indices.push_back(i);
                    }
                    t_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
                    theta_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
                    phi_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
                } else {
                    t_f_vec[i] = forward->t_f - boost::lexical_cast<Real>(item[6]);
                    theta_f_vec[i] = forward->theta_f - boost::lexical_cast<Real>(item[7]);
                    phi_f_vec[i] = forward->phi_f - boost::lexical_cast<Real>(item[8]);
                }
            } catch (std::exception &ex) {
                fmt::println("[{}, {}, {}] Exception: {}", GET_SIGN(nu_r), GET_SIGN(nu_theta), i, ex.what());

                t_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
                theta_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
                phi_f_vec[i] = std::numeric_limits<Real>::quiet_NaN();
            }
        }
    });

    // sort error_indices
    std::sort(error_indices.begin(), error_indices.end());

    // how many elements smaller than error limit
    fmt::println("nu_r: {}, nu_theta: {}", GET_SIGN(nu_r), GET_SIGN(nu_theta));
    fmt::println("error indices: {}", error_indices.size());
    if (!error_indices.empty()) {
        fmt::println("{}", fmt::join(error_indices, ", "));
    }
    fmt::println("t_f: {} / {}, max error: {}", (t_f_vec.array().abs() < ErrorLimit<Real>::Value).count(),
                 t_f_vec.size(), t_f_vec.cwiseAbs().maxCoeff());
    fmt::println("theta_f: {} / {}, max error: {}", (theta_f_vec.array().abs() < ErrorLimit<Real>::Value).count(),
                 theta_f_vec.size(), theta_f_vec.cwiseAbs().maxCoeff());
    fmt::println("phi_f: {} / {}", (phi_f_vec.array().abs() < ErrorLimit<Real>::Value).count(), phi_f_vec.size(),
                 phi_f_vec.cwiseAbs().maxCoeff());

    CHECK(t_f_vec.array().abs().maxCoeff() < ErrorLimit<Real>::Value * 100000000);
    CHECK(theta_f_vec.array().abs().maxCoeff() < ErrorLimit<Real>::Value * 100000000);
    CHECK(phi_f_vec.array().abs().maxCoeff() < ErrorLimit<Real>::Value * 100000000);
}

TEMPLATE_TEST_CASE("Forward Function", "[forward]", TEST_TYPES) {
    using Real = std::tuple_element_t<0u, TestType>;
    using Complex = std::tuple_element_t<1u, TestType>;

    fmt::println("[{}, {}] Error limit: {}", TypeName<Real>::Get(), TypeName<Complex>::Get(), ErrorLimit<Real>::Value);

    SECTION("nu_r = POSITIVE, nu_theta = POSITIVE") {
        test_case<Real, Complex>(TEST_DATA_PP, Sign::POSITIVE, Sign::POSITIVE);
    }
    SECTION("nu_r = POSITIVE, nu_theta = NEGATIVE") {
        test_case<Real, Complex>(TEST_DATA_PM, Sign::POSITIVE, Sign::NEGATIVE);
    }
    SECTION("nu_r = NEGATIVE, nu_theta = POSITIVE") {
        test_case<Real, Complex>(TEST_DATA_MP, Sign::NEGATIVE, Sign::POSITIVE);
    }
    SECTION("nu_r = NEGATIVE, nu_theta = NEGATIVE") {
        test_case<Real, Complex>(TEST_DATA_MM, Sign::NEGATIVE, Sign::NEGATIVE);
    }
}

//TEMPLATE_TEST_CASE("Find Root Function", "[root]", TEST_TYPES) {
//  using Real = std::tuple_element_t<0u, TestType>;
//  using Complex = std::tuple_element_t<1u, TestType>;
//
//  for (const auto &data : TEST_DATA) {
//    ForwardRayTracingParams<Real> params;
//    params.a = boost::lexical_cast<Real>(data.get<string>("a"));
//    params.r_s = boost::lexical_cast<Real>(data.get<string>("r_s"));
//    params.theta_s = boost::lexical_cast<Real>(data.get<string>("theta_s"));
//    params.r_o = boost::lexical_cast<Real>(data.get<string>("r_o"));
//    params.nu_r = data.get<int>("nu_r") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
//    params.nu_theta = data.get<int>("nu_theta") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
//
//    Real lambda = boost::lexical_cast<Real>(data.get<string>("lambda"));
//    Real eta = boost::lexical_cast<Real>(data.get<string>("eta"));
//
//    Real rc = boost::lexical_cast<Real>(data.get<string>("rc"));
//    Real d = boost::lexical_cast<Real>(data.get<string>("d"));
//    params.log_abs_d_sign = d > 0 ? Sign::POSITIVE : Sign::NEGATIVE;
//    Real log_abs_d = log10(abs(d));
//    params.calc_t_f = false;
//
//    Real theta_o = boost::lexical_cast<Real>(data.get<string>("theta_f"));
//    Real phi_o = boost::lexical_cast<Real>(data.get<string>("phi_f"));
//
//    params.rc = rc + 0.001;
//    params.log_abs_d = log_abs_d - 0.001;
//    int period = MY_FLOOR<Real>::convert(theta_o / boost::math::constants::two_pi<Real>());
//    auto root_res = ForwardRayTracingUtils<Real, Complex>::find_root_period(params, period, theta_o, phi_o);
//    REQUIRE(root_res.success);
//    auto &res = *(root_res.root);
//    CHECK(res.ray_status == RayStatus::NORMAL);
//
//    Real ERROR_LIMIT = ErrorLimit<Real>::Value * 1000;
//
//    fmt::println("eta: {}, lambda: {}", eta, lambda);
//    fmt::println("eta: {}, lambda: {}", res.eta, res.lambda);
//    fmt::println("delta eta: {}, delta lambda: {}", abs(eta - res.eta), abs(lambda - res.lambda));
//    CHECK(abs(eta - res.eta) < ERROR_LIMIT);
//    CHECK(abs(lambda - res.lambda) < ERROR_LIMIT);
//
//    auto r1_vec = as_vector<string>(data, "r1");
//    Complex r1 = Complex(boost::lexical_cast<Real>(r1_vec[0]), boost::lexical_cast<Real>(r1_vec[1]));
//    auto r2_vec = as_vector<string>(data, "r2");
//    Complex r2 = Complex(boost::lexical_cast<Real>(r2_vec[0]), boost::lexical_cast<Real>(r2_vec[1]));
//    auto r3_vec = as_vector<string>(data, "r3");
//    Complex r3 = Complex(boost::lexical_cast<Real>(r3_vec[0]), boost::lexical_cast<Real>(r3_vec[1]));
//    auto r4_vec = as_vector<string>(data, "r4");
//    Complex r4 = Complex(boost::lexical_cast<Real>(r4_vec[0]), boost::lexical_cast<Real>(r4_vec[1]));
//
//    Real theta_f = boost::lexical_cast<Real>(data.get<string>("theta_f"));
//    Real phi_f = boost::lexical_cast<Real>(data.get<string>("phi_f"));
//
//    CHECK(abs(res.r1_c - r1) < ERROR_LIMIT);
//    CHECK(abs(res.r2_c - r2) < ERROR_LIMIT);
//    CHECK(abs(res.r3_c - r3) < ERROR_LIMIT);
//    CHECK(abs(res.r4_c - r4) < ERROR_LIMIT);
//    CHECK(abs(res.theta_f - theta_f) < ERROR_LIMIT);
//    CHECK(abs(res.phi_f - phi_f) < ERROR_LIMIT);
//  }
//}