#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


#include "Utils.h"

inline std::vector<boost::property_tree::ptree> TEST_DATA;
void get_test_data(std::string &path);

using Test64 = std::tuple<double, std::complex<double>>;
using Test128 = std::tuple<Float128, Complex128>;
using Test256 = std::tuple<Float256, Complex256>;

#define TEST_TYPES Test64, Test128, Test256
