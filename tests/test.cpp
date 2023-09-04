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
  Real theta_s = 85 * boost::math::constants::pi<Real>() / 180;
  Real r_o = boost::lexical_cast<Real>("1000");
  Sign nu_r = Sign::NEGATIVE;
  Sign nu_theta = Sign::NEGATIVE;

  Real lambda = boost::lexical_cast<Real>("-0.731502407319534083613773987326004666314914414040701439754512574606612766652829");
  Real eta = boost::lexical_cast<Real>("26.56713111403703065164311251315681994126311681197727971968634598181153818059136");

  ForwardRayTracing<Real, Complex> forward(a, r_s, theta_s, r_o);
  forward.calc_t_f = true;
  auto status = forward.calc_ray_by_lambda_q(lambda, sqrt(eta), nu_r, nu_theta);
  REQUIRE(status == RayStatus::NORMAL);

  Real theta_o = 17 * boost::math::constants::pi<Real>() / 180;
  Real phi_o = boost::math::constants::pi<Real>() / 4;

  Real phi_f = boost::lexical_cast<Real>("-5.49778714378213816730962592073913004734504644890643518670615303653867871100086575");

  Real theta_p = boost::lexical_cast<Real>("3.0022293217454433556239461277650734910842674677850392541644705513025725256775");
  Real theta_m = boost::lexical_cast<Real>("0.1393633318443498828386972555144293931129019315900665668104740410052438806087");

  // {-6.04147718092362125054643971198534121943543777647776939740166176118788875572961,0.34976541948735184842467414106274102580862314222399429077625428266783415920229,2.6168580554264530185672062303512116038763894932178737111007094277275553075681,3.0748537060098163835545593405713885897504251410359013955246980507924992889593}
  std::array<Real, 4> roots;
  roots[0] = boost::lexical_cast<Real>("-6.04147718092362125054643971198534121943543777647776939740166176118788875572961");
  roots[1] = boost::lexical_cast<Real>("0.34976541948735184842467414106274102580862314222399429077625428266783415920229");
  roots[2] = boost::lexical_cast<Real>("2.6168580554264530185672062303512116038763894932178737111007094277275553075681");
  roots[3] = boost::lexical_cast<Real>("3.0748537060098163835545593405713885897504251410359013955246980507924992889593");

  std::array<Real, 3> radial_integral;
  radial_integral[0] = boost::lexical_cast<Real>("1.433501833183245888566100056540013190533965940293184660020030286950625692472");
  radial_integral[1] = boost::lexical_cast<Real>("1.261919912675584169592042845980479096921937800252921639946564428836950");
  radial_integral[2] = boost::lexical_cast<Real>("1049.780426179475588631603631554404912236398484925016942523912140358002099");

  std::array<Real, 3> angular_integral;
  angular_integral[0] = boost::lexical_cast<Real>("1.433501833183245888566100056540013190533965940293184660020030286950625692472215318033250206865704");
  angular_integral[1] = boost::lexical_cast<Real>("9.240854150060171263733998416165055249550570577097752258451983153800429441636281690694659074595944");
  angular_integral[2] = boost::lexical_cast<Real>("0.685701771865994128202517526289971980646599537219338301723052463245343065236549586433762191681527");

  // is big float
#ifdef FLOAT128
  if constexpr (std::is_same_v<TestType, FLOAT128>) {

  }
#endif
#ifdef BIGFLOAT
  if constexpr (std::is_same_v<TestType, BigFloat>) {
    CHECK(abs(forward.r1 - roots[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.r2 - roots[1]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.r3 - roots[2]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.r4 - roots[3]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.theta_p - theta_p) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.theta_m - theta_m) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[0] - radial_integral[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[0] - radial_integral[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[1] - radial_integral[1]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.radial_integrals[2] - radial_integral[2]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[0] - angular_integral[0]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[1] - angular_integral[1]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.angular_integrals[2] - angular_integral[2]) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.theta_f - theta_o) < BIGFLOAT_ERROR_LIMIT);
    CHECK(abs(forward.phi_f - phi_f) < BIGFLOAT_ERROR_LIMIT);
  }
#endif
}
