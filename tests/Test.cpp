#include <tuple>
#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestData.h"
#include "ForwardRayTracing.h"
#include "Utils.h"

// ErrorLimit
template <typename T>
struct ErrorLimit {
  static const T Value;
};

using std::string;
using Float64 = std::tuple<double, std::complex<double>>;

template <>
const double ErrorLimit<double>::Value = 1e-10;

#ifdef FLOAT128
using Float128 = std::tuple<boost::multiprecision::float128, boost::multiprecision::complex128>;

template <>
const boost::multiprecision::float128 ErrorLimit<boost::multiprecision::float128>::Value{"1e-29"};
#endif

using BigFloat = std::tuple<BigFloatReal , BigFloatComplex>;
template <>
const BigFloatReal ErrorLimit<BigFloatReal>::Value{"1e-45"};

#if defined(FLOAT128)
#define TEST_TYPES Float64, Float128, BigFloat
#else
#define TEST_TYPES Float64, BigFloat
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
    ForwardRayTracingParams<Real> params;
    params.a = boost::lexical_cast<Real>(data.get<string>("a"));
    params.r_s = boost::lexical_cast<Real>(data.get<string>("r_s"));
    params.theta_s = boost::lexical_cast<Real>(data.get<string>("theta_s"));
    params.r_o = boost::lexical_cast<Real>(data.get<string>("r_o"));
    params.nu_r = data.get<int>("nu_r") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
    params.nu_theta = data.get<int>("nu_theta") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;

    params.lambda = boost::lexical_cast<Real>(data.get<string>("lambda"));
    params.q = sqrt(boost::lexical_cast<Real>(data.get<string>("eta")));
    params.calc_t_f = true;

    ForwardRayTracing<Real, Complex> forward;
    forward.calc_ray(params);
    CHECK(forward.ray_status == RayStatus::NORMAL);

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

    Real ERROR_LIMIT = ErrorLimit<Real>::Value;

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

TEMPLATE_TEST_CASE("Find Root Function", "[root]", Float64) {
  using Real = std::tuple_element_t<0u, TestType>;
  using Complex = std::tuple_element_t<1u, TestType>;

  for (const auto &data : TEST_DATA) {
    ForwardRayTracingParams<Real> params;
    params.a = boost::lexical_cast<Real>(data.get<string>("a"));
    params.r_s = boost::lexical_cast<Real>(data.get<string>("r_s"));
    params.theta_s = boost::lexical_cast<Real>(data.get<string>("theta_s"));
    params.r_o = boost::lexical_cast<Real>(data.get<string>("r_o"));
    params.nu_r = data.get<int>("nu_r") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
    params.nu_theta = data.get<int>("nu_theta") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;

    params.rc = boost::lexical_cast<Real>(data.get<string>("rc"));
    Real d = boost::lexical_cast<Real>(data.get<string>("d"));
    params.lgd_sign = d > 0 ? Sign::POSITIVE : Sign::NEGATIVE;
    params.lgd = log10(abs(d));
    params.calc_t_f = false;

    auto res = ForwardRayTracingUtils<Real, Complex>::find_result(params);
    REQUIRE(res.ray_status == RayStatus::NORMAL);

    auto r1_vec = as_vector<string>(data, "r1");
    Complex r1 = Complex(boost::lexical_cast<Real>(r1_vec[0]), boost::lexical_cast<Real>(r1_vec[1]));
    auto r2_vec = as_vector<string>(data, "r2");
    Complex r2 = Complex(boost::lexical_cast<Real>(r2_vec[0]), boost::lexical_cast<Real>(r2_vec[1]));
    auto r3_vec = as_vector<string>(data, "r3");
    Complex r3 = Complex(boost::lexical_cast<Real>(r3_vec[0]), boost::lexical_cast<Real>(r3_vec[1]));
    auto r4_vec = as_vector<string>(data, "r4");
    Complex r4 = Complex(boost::lexical_cast<Real>(r4_vec[0]), boost::lexical_cast<Real>(r4_vec[1]));

    Real theta_f = boost::lexical_cast<Real>(data.get<string>("theta_f"));
    Real phi_f = boost::lexical_cast<Real>(data.get<string>("phi_f"));

    Real ERROR_LIMIT = ErrorLimit<Real>::Value;

    CHECK(abs(res.r1_c - r1) < ERROR_LIMIT);
    CHECK(abs(res.r2_c - r2) < ERROR_LIMIT);
    CHECK(abs(res.r3_c - r3) < ERROR_LIMIT);
    CHECK(abs(res.r4_c - r4) < ERROR_LIMIT);
    CHECK(abs(res.theta_f - theta_f) < ERROR_LIMIT);
    CHECK(abs(res.phi_f - phi_f) < ERROR_LIMIT);
  }
}