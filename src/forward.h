#pragma once

#include <array>

// using Real = double;
// #define RCONST(x) x

using Real = long double;
#define RCONST(x) x##L

enum class RayStatus {
  NORMAL,
  FALLS_IN,
  CONFINED,
  ETA_OUT_OF_RANGE, // eta should be positive
  THETA_OUT_OF_RANGE, // theta should be in [theta_m, theta_p]
};

enum class Sign {
  POSITIVE,
  NEGATIVE,
};

// convert rc, d to lambda, q. Here the vertical distance d is defined on the (lambda, q) plane, q = sqrt(eta)
std::array<Real, 2> lamqv2(Real a, Real r_c, Real d);

// convert rc, d to lambda, eta
std::array<Real, 2> lametav2(Real a, Real r_c, Real d);

// find roots
std::array<Real, 4> roots(Real a, Real lambda, Real eta);

std::pair<RayStatus, std::array<Real, 3>> calc_I(Real a, Real lambda, Real eta, Real r_s, Real r_inf, Real nu_r);

// calculate angular path integral: G_theta, G_phi, G_t
std::array<Real, 3> calc_G(Real a, Real lambda, Real eta, Real theta_s, Real theta_inf, int m, Real nu_theta);
