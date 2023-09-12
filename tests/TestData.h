#pragma once

#include <boost/lexical_cast.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Utils.h"

#include <vector>
#include <array>
#include <string>

inline std::vector<std::array<std::string, 9>> TEST_DATA_PP;
inline std::vector<std::array<std::string, 9>> TEST_DATA_PM;
inline std::vector<std::array<std::string, 9>> TEST_DATA_MP;
inline std::vector<std::array<std::string, 9>> TEST_DATA_MM;

void get_test_data(std::string &path);

using Test64 = std::tuple<double, std::complex<double>>;
using Test128 = std::tuple<Float128, Complex128>;
using Test256 = std::tuple<Float256, Complex256>;

#define TEST_TYPES Test64, Test128, Test256
