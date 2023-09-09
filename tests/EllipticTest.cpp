
#include <catch2/catch_template_test_macros.hpp>

#include "TestData.h"

TEMPLATE_TEST_CASE("Elliptic Integral Function", "[elliptic]", TEST_TYPES) {
	using Real = std::tuple_element_t<0u, TestType>;
	using Complex = std::tuple_element_t<1u, TestType>;
	const Real half_pi = boost::math::constants::half_pi<Real>();
	//Eigen::Vector<Real, 100> phi_list = Eigen::Vector<Real, 100>::LinSpaced(100, -half_pi, half_pi);
	//for (auto &phi : phi_list)
	//{
	//	Vector m_list = Vector::LinSpaced(100, 0.0, 1.0);
	//}
}
