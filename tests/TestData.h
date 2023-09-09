#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


#include "Utils.h"

inline std::vector<boost::property_tree::ptree> TEST_DATA;
void get_test_data(std::string &path);

using Float64 = std::tuple<double, std::complex<double>>;
#ifdef FLOAT128
using Float128 = std::tuple<boost::multiprecision::float128, boost::multiprecision::complex128>;
#endif
using BigFloat = std::tuple<BigFloatReal, BigFloatComplex>;

#if defined(FLOAT128)
#define TEST_TYPES Float64, Float128, BigFloat
#else
#define TEST_TYPES Float64, BigFloat
#endif