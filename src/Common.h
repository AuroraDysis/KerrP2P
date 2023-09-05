#pragma once

#include <cmath>
#include <array>
#include <complex>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>


using boost::math::ellint_1;
using boost::math::ellint_2;
using boost::math::ellint_3;
using boost::math::constants::half;
using boost::math::constants::third;
using boost::math::constants::sixth;

enum class RayStatus {
  NORMAL,
  FALLS_IN,
  CONFINED,
  ETA_OUT_OF_RANGE, // eta should be positive
  THETA_OUT_OF_RANGE, // theta should be in [theta_m, theta_p]
  UNKOWN_ERROR,
};

enum class Sign : int {
  POSITIVE = 1,
  NEGATIVE = -1,
};

#define square(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))

template<typename E>
constexpr auto to_integral(E e) -> typename std::underlying_type<E>::type {
  return static_cast<typename std::underlying_type<E>::type>(e);
}

template<typename Real, typename Complex>
class ForwardRayTracing;

using std::real;
using std::isinf;

#ifdef FLOAT128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>

template <>
struct fmt::formatter<boost::multiprecision::float128> : fmt::ostream_formatter {};
#endif

#ifdef BIGFLOAT
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#endif

#if defined(FLOAT128) || defined(BIGFLOAT)
using boost::multiprecision::isinf;
using boost::multiprecision::real;

template <typename T>
struct fmt::formatter<boost::multiprecision::number<T>> : fmt::ostream_formatter {};
#endif
