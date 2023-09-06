
#pragma once

// version

#ifndef BMO_VERSION_MAJOR
    #define BMO_VERSION_MAJOR 0
#endif

#ifndef BMO_VERSION_MINOR
    #define BMO_VERSION_MINOR 2
#endif

#ifndef BMO_VERSION_PATCH
    #define BMO_VERSION_PATCH 1
#endif

//

#ifdef BMO_ENABLE_EIGEN_WRAPPERS

    #if EIGEN_VERSION_AT_LEAST(3,4,50)
        #define BMO_EIGEN_INDEX_ALL Eigen::indexing::all
    #elif EIGEN_VERSION_AT_LEAST(3,4,0)
        #define BMO_EIGEN_INDEX_ALL Eigen::all
    #else
        #error Eigen must be version 3.4.0 or above
    #endif

    using scalar_t = bmo::ColVec_t::Scalar;
#endif

#ifdef BMO_ENABLE_ARMA_WRAPPERS
    using scalar_t = bmo::ColVec_t::elem_type;
#endif
