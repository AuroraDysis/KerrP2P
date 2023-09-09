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
  ARGUMENT_ERROR,
  INTERNAL_ERROR, // may be caused by not enough precision
  UNKOWN_ERROR,
};

constexpr const char *ray_status_to_str(RayStatus status) {
  switch (status) {
    case RayStatus::NORMAL:
      return "NORMAL";
    case RayStatus::FALLS_IN:
      return "FALLS_IN";
    case RayStatus::CONFINED:
      return "CONFINED";
    case RayStatus::ETA_OUT_OF_RANGE:
      return "ETA_OUT_OF_RANGE";
    case RayStatus::THETA_OUT_OF_RANGE:
      return "THETA_OUT_OF_RANGE";
    case RayStatus::ARGUMENT_ERROR:
      return "ARGUMENT_ERROR";
    case RayStatus::INTERNAL_ERROR:
      return "INTERNAL_ERROR";
    case RayStatus::UNKOWN_ERROR:
      return "UNKOWN_ERROR";
  }
  return "UNKOWN_ERROR";
}

enum class Sign : int {
  POSITIVE = 1,
  NEGATIVE = -1,
};

#define MY_SQUARE(x) ((x) * (x))
#define MY_CUBE(x) ((x) * (x) * (x))

template<typename E>
constexpr auto to_integral(E e) -> typename std::underlying_type<E>::type {
  return static_cast<typename std::underlying_type<E>::type>(e);
}

template<typename Real, typename Complex>
class ForwardRayTracing;

using std::real;
using std::isinf;
using std::isnan;

#ifdef FLOAT128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>

template <>
struct fmt::formatter<boost::multiprecision::float128> : fmt::ostream_formatter {};
#endif

#ifdef BIGFLOAT_MPFR
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

using BigFloatReal = boost::multiprecision::mpfr_float_50;
using BigFloatComplex = boost::multiprecision::mpc_complex_50;
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using BigFloatReal = boost::multiprecision::cpp_bin_float_50;
using BigFloatComplex = boost::multiprecision::cpp_complex_50;;
#endif

using boost::multiprecision::real;
using boost::multiprecision::isinf;
using boost::multiprecision::isnan;

template <typename T>
struct fmt::formatter<boost::multiprecision::number<T>> : fmt::ostream_formatter {};

// helper return higher precision type
template <typename T>
struct HigherPrecision {
  using Type = T;
};

#if defined(FLOAT128)
template <>
struct HigherPrecision<double> {
  using Type = boost::multiprecision::float128;
};

template <>
struct HigherPrecision<std::complex<double>> {
  using Type = boost::multiprecision::complex128;
};

template <>
struct HigherPrecision<boost::multiprecision::float128> {
  using Type = boost::multiprecision::float128;
};

template <>
struct HigherPrecision<std::complex<boost::multiprecision::float128>> {
  using Type = BigFloatComplex;
};
#else
template <>
struct HigherPrecision<double> {
  using Type = BigFloatReal;
};

template <>
struct HigherPrecision<std::complex<double>> {
  using Type = BigFloatComplex;
};
#endif

// ErrorLimit
template <typename T>
struct ErrorLimit {
  inline static const T Value = std::numeric_limits<T>::epsilon() * 1000000;
};
