#include "forward.h"

#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

//bool compare(Real expected) {
//  return Catch::Matchers::WithinRel(expected, 1e-15) || Catch::Matchers::WithinAbs(expected, 1e-15)
//}

auto match_float(Real expected) {
  return ;
}

TEST_CASE("Forward Functions", "[forward]") {
  SECTION("lamqv2") {
    Real a = RCONST(0.9);
    Real rc = RCONST(3);
    Real d = RCONST(-1);
    auto res = lamqv2(a, rc, d);
    std::array<Real, 2> expected = {-RCONST(9.0) / RCONST(5.0), RCONST(-1.0) + RCONST(3.0) * std::sqrt(RCONST(3.0))};
    for (int i = 0; i < 2; ++i) {
      REQUIRE_THAT(res[i],
                   Catch::Matchers::WithinRel(expected[i], 1e-15) || Catch::Matchers::WithinAbs(expected[i], 1e-15));
    }
  }

  SECTION("lametav2") {
    Real a = RCONST(0.9);
    Real rc = RCONST(3);
    Real d = RCONST(-1);
    auto res = lametav2(a, rc, d);
    std::array<Real, 2> expected = {-RCONST(9.0) / RCONST(5.0), RCONST(28.0) - RCONST(6.0) * std::sqrt(RCONST(3.0))};
    for (int i = 0; i < 2; ++i) {
      REQUIRE_THAT(res[i],
                   Catch::Matchers::WithinRel(expected[i], 1e-15) || Catch::Matchers::WithinAbs(expected[i], 1e-15));
    }
  }
}

