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
  forward.calc_ray_by_lambda_q(lambda, sqrt(eta), nu_r, nu_theta);

  Real theta_o = 17 * boost::math::constants::pi<Real>() / 180;
  Real phi_o = boost::math::constants::pi<Real>() / 4;
  REQUIRE_THAT((forward.theta_f - theta_o).template convert_to<double>(), Catch::Matchers::WithinAbs(0, 1e-15));

//    ForwardRayTracing forward();
//    SECTION("lamqv2") {
//        Real a = RCONST(0.9);
//        Real rc = RCONST(3);
//        Real d = RCONST(-1);
//        auto res = lamqv2(a, rc, d);
//        std::array<Real, 2> expected = {-RCONST(9.0) / RCONST(5.0),
//                                        RCONST(-1.0) + RCONST(3.0) * std::sqrt(RCONST(3.0))};
//        for (int i = 0; i < 2; ++i) {
//            REQUIRE_THAT(res[i],
//                         Catch::Matchers::WithinRel(expected[i], 1e-15) ||
//                         Catch::Matchers::WithinAbs(expected[i], 1e-15));
//        }
//    }
//
//    SECTION("lametav2") {
//        Real a = RCONST(0.9);
//        Real rc = RCONST(3);
//        Real d = RCONST(-1);
//        auto res = lametav2(a, rc, d);
//        std::array<Real, 2> expected = {-RCONST(9.0) / RCONST(5.0),
//                                        RCONST(28.0) - RCONST(6.0) * std::sqrt(RCONST(3.0))};
//        for (int i = 0; i < 2; ++i) {
//            REQUIRE_THAT(res[i],
//                         Catch::Matchers::WithinRel(expected[i], 1e-15) ||
//                         Catch::Matchers::WithinAbs(expected[i], 1e-15));
//        }
//    }
}

