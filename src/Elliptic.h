#pragma once

#include "Common.h"

template <typename Real>
class EllipticIntegral {
private:
    Real x, y, z, rho;
public:
    Real ellint_k(const Real &phi, const Real &cos_phi, const Real sin_phi, const Real &m) {
    }

    Real ellint_f(const Real &phi, const Real &cos_phi, const Real sin_phi, const Real &m) {
    }

    Real ellint_pi(const Real &phi, const Real &cos_phi, const Real sin_phi, const Real &m) {
        // https://github.com/achael/kgeo/blob/main/kgeo/gsl_ellip_binding.py
        x = MY_SQUARE(cos_phi);
        y = 1 - m * MY_SQUARE(sin_phi);
        z = 1;
        rho = 1 - MY_SQUARE(sin_phi);

        if (abs(phi) > boost::math::constants::pi<Real>()) {
        }
    }
};