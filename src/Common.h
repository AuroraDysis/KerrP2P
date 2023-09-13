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

template<typename Real>
using MY_FLOOR = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>, boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;

template<typename Real, typename Complex>
class ForwardRayTracing;

using std::real;
using std::isinf;
using std::isnan;

template <typename T>
struct TypeName
{
    static std::string Get()
    {
        return typeid(T).name();
    }
};

#ifdef FLOAT128_NATIVE
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>

using Float128 = boost::multiprecision::float128;
using Complex128 = boost::multiprecision::complex128;

template <>
struct TypeName<boost::multiprecision::float128>
{
    static std::string Get()
    {
        return "float128";
    }
};

template <>
struct TypeName<boost::multiprecision::complex128>
{
    static std::string Get()
    {
		return "complex128";
	}
};
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using Float128 = boost::multiprecision::cpp_bin_float_quad;
using Complex128 = boost::multiprecision::cpp_complex_quad;

template <>
struct TypeName<boost::multiprecision::cpp_bin_float_quad>
{
	static std::string Get()
	{
		return "cpp_bin_float_quad";
	}
};

template <>
struct TypeName<boost::multiprecision::cpp_complex_quad>
{
	static std::string Get()
	{
		return "cpp_complex_quad";
	}
};
#endif

#ifdef ENABLE_MPFR
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

using mpfr_float_oct = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<71>>;
using mpc_complex_oct = boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<71>>;

using Float256 = mpfr_float_oct;
using Complex256 = mpc_complex_oct;

template <>
struct TypeName<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<71>>>
{
    static std::string Get()
    {
		return "mpfr_float_oct";
	}
};

template <>
struct TypeName<mpc_complex_oct>
{
    static std::string Get()
    {
		return "mpc_complex_oct";
	}
};
#else
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using Float256 = boost::multiprecision::cpp_bin_float_oct;
using Complex256 = boost::multiprecision::cpp_complex_oct;

template <>
struct TypeName<boost::multiprecision::cpp_bin_float_oct>
{
    static std::string Get()
    {
		return "cpp_bin_float_oct";
	}
};

template <>
struct TypeName<boost::multiprecision::cpp_complex_oct>
{
    static std::string Get()
    {
		return "cpp_complex_oct";
	}
};
#endif

using boost::multiprecision::real;
using boost::multiprecision::isinf;
using boost::multiprecision::isnan;

template<typename T>
struct fmt::formatter<boost::multiprecision::number<T>> : fmt::ostream_formatter {
};

// helper return higher precision type

template<typename T>
struct HigherPrecision {
    using Type = T;
};

template <>
struct HigherPrecision<double> {
    using Type = Float128;
};

template <>
struct HigherPrecision<std::complex<double>> {
    using Type = Complex128;
};

template <>
struct HigherPrecision<Float128> {
    using Type = Float256;
};

template <>
struct HigherPrecision<Complex128> {
    using Type = Complex256;
};

#ifdef ENABLE_MPFR
template <unsigned digits10>
struct HigherPrecision<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10>>> {
    using Type = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10 + 20>>;
};

template <unsigned digits10>
struct HigherPrecision<boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<digits10>>> {
    using Type = boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<digits10 + 20>>;
};

template <unsigned digits10>
struct TypeName<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<digits10>>>
{
    static std::string Get()
    {
        return fmt::format("mpfr_float_{}", digits10).c_str();
    }
};

template <unsigned digits10>
struct TypeName<boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<digits10>>>
{
    static std::string Get()
    {
        return fmt::format("mpc_complex_{}", digits10).c_str();
    }
};
#else
template <>
struct HigherPrecision<boost::multiprecision::cpp_bin_float_oct> {
    using Type = boost::multiprecision::cpp_bin_float_100;
};

template <>
struct HigherPrecision<boost::multiprecision::cpp_complex_oct> {
    using Type = boost::multiprecision::cpp_complex_100;
};

template <unsigned digits10>
struct HigherPrecision<boost::multiprecision::cpp_bin_float<digits10>> {
    using Type = boost::multiprecision::cpp_bin_float<digits10 + 20>;
};

template <unsigned digits10>
struct HigherPrecision<boost::multiprecision::cpp_complex<digits10>> {
    using Type = boost::multiprecision::cpp_complex<digits10 + 20>;
};

template <unsigned digits10>
struct TypeName<boost::multiprecision::cpp_bin_float<digits10>>
{
    static std::string Get()
    {
        return fmt::format("cpp_bin_float_{}", digits10);
    }
};

template <unsigned digits10>
struct TypeName<boost::multiprecision::cpp_complex<digits10>>
{
    static std::string Get()
    {
        return fmt::format("cpp_complex_{}", digits10);
    }
};
#endif

// ErrorLimit
template<typename T>
struct ErrorLimit {
    inline static const T Value = std::numeric_limits<T>::epsilon() * 1000000;
};
