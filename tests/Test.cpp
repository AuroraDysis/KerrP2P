#include "ForwardRayTracing.h"

#include <tuple>
#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestData.h"

using std::string;
using Float64 = std::tuple<double, std::complex<double>>;
constexpr double FLOAT64_ERROR_LIMIT = 1e-12;

#ifdef FLOAT128
using Float128 = std::tuple<boost::multiprecision::float128, boost::multiprecision::complex128>;
const boost::multiprecision::float128 FLOAT128_ERROR_LIMIT{"1e-30"};
#endif

#ifdef BIGFLOAT
using BigFloat = std::tuple<boost::multiprecision::mpfr_float_50, boost::multiprecision::mpc_complex_50>;
const boost::multiprecision::mpfr_float_50 BIGFLOAT_ERROR_LIMIT{"1e-45"};
#endif

#if defined(FLOAT128) && defined(BIGFLOAT)
#define TEST_TYPES Float64, Float128, BigFloat
#elif defined(FLOAT128) && !defined(BIGFLOAT)
#define TEST_TYPES Float64, Float128
#elif !defined(FLOAT128) && defined(BIGFLOAT)
#define TEST_TYPES Float64, BigFloat
#else
#define TEST_TYPES Float64
#endif

template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const& pt, boost::property_tree::ptree::key_type const& key) {
  std::vector<T> r;
  for (auto &item: pt.get_child(key))
    r.push_back(item.second.get_value<T>());
  return r;
}

TEMPLATE_TEST_CASE("Forward Function", "[forward]", TEST_TYPES) {
  using Real = std::tuple_element_t<0u, TestType>;
  using Complex = std::tuple_element_t<1u, TestType>;

  for (const auto &data : TEST_DATA) {
    Real a = boost::lexical_cast<Real>(data.get<string>("a"));
    Real r_s = boost::lexical_cast<Real>(data.get<string>("r_s"));
    Real theta_s = boost::lexical_cast<Real>(data.get<string>("theta_s"));
    Real r_o = boost::lexical_cast<Real>(data.get<string>("r_o"));
    Sign nu_r = data.get<int>("nu_r") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
    Sign nu_theta = data.get<int>("nu_theta") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;

    Real lambda = boost::lexical_cast<Real>(data.get<string>("lambda"));
    Real eta = boost::lexical_cast<Real>(data.get<string>("eta"));

    ForwardRayTracing<Real, Complex> forward;
    forward.calc_t_f = true;
    auto status = forward.calc_ray_by_lambda_q(a, r_s, theta_s, r_o, nu_r, nu_theta, lambda, sqrt(eta));
    REQUIRE(status == RayStatus::NORMAL);

    auto r1_vec = as_vector<string>(data, "r1");
    Complex r1 = Complex(boost::lexical_cast<Real>(r1_vec[0]), boost::lexical_cast<Real>(r1_vec[1]));
    auto r2_vec = as_vector<string>(data, "r2");
    Complex r2 = Complex(boost::lexical_cast<Real>(r2_vec[0]), boost::lexical_cast<Real>(r2_vec[1]));
    auto r3_vec = as_vector<string>(data, "r3");
    Complex r3 = Complex(boost::lexical_cast<Real>(r3_vec[0]), boost::lexical_cast<Real>(r3_vec[1]));
    auto r4_vec = as_vector<string>(data, "r4");
    Complex r4 = Complex(boost::lexical_cast<Real>(r4_vec[0]), boost::lexical_cast<Real>(r4_vec[1]));

    std::vector<Real> radial_integrals;
    radial_integrals.reserve(3);
    for (const auto& str : as_vector<string>(data, "radial_integrals")) {
      radial_integrals.push_back(boost::lexical_cast<Real>(str));
    }

    std::vector<Real> angular_integrals;
    angular_integrals.reserve(3);
    for (const auto& str : as_vector<string>(data, "angular_integrals")) {
      angular_integrals.push_back(boost::lexical_cast<Real>(str));
    }

    Real theta_f = boost::lexical_cast<Real>(data.get<string>("theta_f"));
    Real phi_f = boost::lexical_cast<Real>(data.get<string>("phi_f"));

    // is big float
    Real ERROR_LIMIT;
    if constexpr (std::is_same_v<TestType, Float64>) {
      ERROR_LIMIT = FLOAT64_ERROR_LIMIT;
    }
#ifdef FLOAT128
    if constexpr (std::is_same_v<TestType, Float128>) {
      ERROR_LIMIT = FLOAT128_ERROR_LIMIT;
    }
#endif
#ifdef BIGFLOAT
    if constexpr (std::is_same_v<TestType, BigFloat>) {
      ERROR_LIMIT = BIGFLOAT_ERROR_LIMIT;
    }
#endif
    CHECK(abs(forward.r1_c - r1) < ERROR_LIMIT);
    CHECK(abs(forward.r2_c - r2) < ERROR_LIMIT);
    CHECK(abs(forward.r3_c - r3) < ERROR_LIMIT);
    CHECK(abs(forward.r4_c - r4) < ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[0] - radial_integrals[0]) < ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[1] - radial_integrals[1]) < ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[2] - radial_integrals[2]) < ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[0] - angular_integrals[0]) < ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[1] - angular_integrals[1]) < ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[2] - angular_integrals[2]) < ERROR_LIMIT);
    CHECK(abs(forward.theta_f - theta_f) < ERROR_LIMIT);
    CHECK(abs(forward.phi_f - phi_f) < ERROR_LIMIT);
  }
}
