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
