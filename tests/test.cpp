#include "Forward.h"

#include <tuple>
#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <boost/math/constants/constants.hpp>

using Float64 = std::tuple<double, std::complex<double>>;

#ifdef FLOAT128
using Float128 = std::tuple<boost::multiprecision::float128, boost::multiprecision::complex128>;
#endif

#ifdef BIGFLOAT
using BigFloat = std::tuple<boost::multiprecision::mpfr_float_50, boost::multiprecision::mpc_complex_50>;
const boost::multiprecision::mpfr_float_50 BIGFLOAT_ERROR_LIMIT{"1e-45"};
#endif

TEMPLATE_TEST_CASE("Forward Function", "[forward]", BigFloat) {
  using Real = std::tuple_element_t<0u, TestType>;
  using Complex = std::tuple_element_t<1u, TestType>;
  Real a = boost::lexical_cast<Real>("0.8");
  Real r_s = boost::lexical_cast<Real>("10");
  Real theta_s = boost::math::constants::half_pi<Real>();
  Real r_o = boost::lexical_cast<Real>("1000");
  Sign nu_r = Sign::NEGATIVE;
  Sign nu_theta = Sign::NEGATIVE;
  Real lambda = boost::lexical_cast<Real>("-0.751131614119698035182184635438749173900590136995321237367912291363607606735312");
  Real eta = boost::lexical_cast<Real>("26.57242896970946927251987627934469966678305414295689726675509814044032238307038");

  ForwardRayTracing<Real, Complex> forward(a, r_s, theta_s, r_o);
  auto status = forward.calc_ray_by_lambda_q(lambda, sqrt(eta), nu_r, nu_theta);
  REQUIRE(status == RayStatus::NORMAL);

  Real theta_o = 17 * boost::math::constants::pi<Real>() / 180;
  Real phi_o = boost::math::constants::pi<Real>() / 4;

  std::array<Real, 3> radial_integral;
  radial_integral[0] = boost::lexical_cast<Real>("1.449918724448707783003361604770409888333642555992048976602177112826476986925901542332371494798805");
  radial_integral[1] = boost::lexical_cast<Real>("1.288507079276017729837985995399653209220514835403767425334993622951859438005362259829818057986322");
  radial_integral[2] = boost::lexical_cast<Real>("1050.230729556838727739610680482910331549495652276193749581856090200250574410573699314114422714724239");

  // {1.449918724448707783003361604770409888333642555992048976602177112826476986925901542332371494800308`74.93029591402417, 9.034760480706797702786289961054367423398533894470186848502176405935371833747610240436842225095523`72.26327202415743, 0.684953649593118508702107588106324458927421602634607435896061608959165805326668612330938804832442`72.54800818746092}
  std::array<Real, 3> angular_integral;
  angular_integral[0] = boost::lexical_cast<Real>("1.449918724448707783003361604770409888333642555992048976602177112826476986925901542332371494800308");
  angular_integral[1] = boost::lexical_cast<Real>("9.034760480706797702786289961054367423398533894470186848502176405935371833747610240436842225095523");
  angular_integral[2] = boost::lexical_cast<Real>("0.684953649593118508702107588106324458927421602634607435896061608959165805326668612330938804832442");

  // is big float
#ifdef FLOAT128
  if constexpr (std::is_same_v<TestType, FLOAT128>) {

  }
#endif
#ifdef BIGFLOAT
  if constexpr (std::is_same_v<TestType, BigFloat>) {
    CHECK(abs(forward.radial_integrals[0] - radial_integral[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[1] - radial_integral[1]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[2] - radial_integral[2]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[0] - angular_integral[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[1] - angular_integral[1]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[2] - angular_integral[2]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.theta_f - theta_o) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.phi_f - phi_o) < BIGFLOAT_ERROR_LIMIT);
  }
#endif
}
