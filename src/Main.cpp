#include <string>
#include <iostream>

#include <boost/filesystem.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace boost::property_tree;

#include "ForwardRayTracing.h"
#include "Utils.h"

using std::string;
using Float64 = std::tuple<double, std::complex<double>>;

#ifdef FLOAT128
using Float128 = std::tuple<boost::multiprecision::float128, boost::multiprecision::complex128>;
#endif

using BigFloat = std::tuple<BigFloatReal , BigFloatComplex>;

#if defined(FLOAT128)
#define TEST_TYPES Float64, Float128, BigFloat
#else
#define TEST_TYPES Float64, BigFloat
#endif

inline std::vector<boost::property_tree::ptree> TEST_DATA;
void get_test_data(std::string &path);
template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const& pt, boost::property_tree::ptree::key_type const& key) {
  std::vector<T> r;
  for (auto &item: pt.get_child(key))
    r.push_back(item.second.get_value<T>());
  return r;
}

#define CHECK(x) if (!(x)) { \
  std::cout << "Test failed: " << #x << std::endl; \
  std::cout << "File: " << __FILE__ << std::endl; \
  std::cout << "Line: " << __LINE__ << std::endl; \
  std::cout << "Function: " << __FUNCTION__ << std::endl; \
  std::cout << "Test data: " << std::endl; \
  return; \
}

template <typename Real, typename Complex>
void test() {
  for (const auto &data : TEST_DATA) {
    ForwardRayTracingParams<Real> params;
    params.a = boost::lexical_cast<Real>(data.get<string>("a"));
    params.r_s = boost::lexical_cast<Real>(data.get<string>("r_s"));
    params.theta_s = boost::lexical_cast<Real>(data.get<string>("theta_s"));
    params.r_o = boost::lexical_cast<Real>(data.get<string>("r_o"));
    params.nu_r = data.get<int>("nu_r") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;
    params.nu_theta = data.get<int>("nu_theta") == 1 ? Sign::POSITIVE : Sign::NEGATIVE;

    Real lambda = boost::lexical_cast<Real>(data.get<string>("lambda"));
    Real eta = boost::lexical_cast<Real>(data.get<string>("eta"));

    Real rc = boost::lexical_cast<Real>(data.get<string>("rc"));
    Real d = boost::lexical_cast<Real>(data.get<string>("d"));
    params.lgd_sign = d > 0 ? Sign::POSITIVE : Sign::NEGATIVE;
    Real lgd = log10(abs(d));
    params.calc_t_f = false;

    Real theta_o = boost::lexical_cast<Real>(data.get<string>("theta_f"));
    Real phi_o = boost::lexical_cast<Real>(data.get<string>("phi_f"));

    params.rc = rc + 0.001;
    params.lgd = lgd - 0.001;
    using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
        boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;
    int period = RealToInt::convert(theta_o / boost::math::constants::two_pi<Real>());
    auto res = ForwardRayTracingUtils<Real, Complex>::find_result(params, period, theta_o, phi_o);
    CHECK(res.ray_status == RayStatus::NORMAL);

    Real ERROR_LIMIT = ErrorLimit<Real>::Value * 1000;

    fmt::println("eta: {}, lambad: {}", eta, lambda);
    fmt::println("eta: {}, lambad: {}", res.eta, res.lambda);
    fmt::println("delta eta: {}, delta lambad: {}", abs(eta - res.eta), abs(lambda - res.lambda));
    CHECK(abs(eta - res.eta) < ERROR_LIMIT);
    CHECK(abs(lambda - res.lambda) < ERROR_LIMIT);

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

    CHECK(abs(res.r1_c - r1) < ERROR_LIMIT);
    CHECK(abs(res.r2_c - r2) < ERROR_LIMIT);
    CHECK(abs(res.r3_c - r3) < ERROR_LIMIT);
    CHECK(abs(res.r4_c - r4) < ERROR_LIMIT);
    CHECK(abs(res.theta_f - theta_f) < ERROR_LIMIT);
    CHECK(abs(res.phi_f - phi_f) < ERROR_LIMIT);
  }
}

void test2() {
  ForwardRayTracingParams<double> params;
  params.a = 0.8;
  params.r_s = 10;
  params.theta_s = 85 * boost::math::constants::pi<double>() / 180;
  params.r_o = 1000;
  params.nu_r = Sign::NEGATIVE;
  params.nu_theta = Sign::NEGATIVE;
  double theta_o = 17 * boost::math::constants::pi<double>() / 180;
  double phi_o = boost::math::constants::pi<double>() / 4;
  params.lgd_sign = Sign::POSITIVE;
  std::vector<double> rc_list;
  std::vector<double> lgd_list;
  rc_list.reserve(100);
  lgd_list.reserve(100);
  for (int i = 0; i < 1000; ++i) {
    rc_list.push_back(2 + i * 0.0015);
  }
  for (int i = 0; i < 1000; ++i) {
    lgd_list.push_back(-10 + i * 0.01);
  }
  auto sweep_result = ForwardRayTracingUtils<double, std::complex<double>>::sweep_rc_d(params, theta_o, phi_o, rc_list,
                                                                                       lgd_list, 50);
  fmt::println("count: {}", sweep_result.results.size());
  for (const auto &res: sweep_result.results) {
    fmt::println("rc: {}, lgd: {}, eta: {}, lambda: {}", res.rc, res.lgd, res.eta, res.lambda);
  }
}

int main(int argc, char **argv) {
//  if (argc < 2) {
//    std::cout << "Usage: " << argv[0] << " <path to test data>" << std::endl;
//    return 1;
//  }
//  std::string path = argv[1];
//  get_test_data(path);
//
//  test<double, std::complex<double>>();
//
  test2();
  return 0;
}

void get_test_data(std::string &path) {
  using namespace boost::filesystem;
  std::vector<ptree> test_data;

  if (!is_directory(path) || !exists(path)) {
    return;
  }

  for(auto& entry : boost::make_iterator_range(directory_iterator(path), {})) {
    // find .json files and parse them
    if (entry.path().extension() == ".json") {
      std::ifstream file(entry.path().string());
      ptree jsontree;
      read_json(file, jsontree);
      TEST_DATA.push_back(std::move(jsontree));
    }
  }
}
