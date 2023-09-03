cmake_minimum_required(VERSION 3.17)
project(forward_ray_tracing)

set(CMAKE_CXX_STANDARD 17)

option(ENABLE_TESTING "Enable Test Builds" ON)
option(PRECISION "Precision: [FLOAT64, FLOAT128, BIGFLOAT]" BIGFLOAT)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost 1.82 REQUIRED)
find_package(Catch2 3 REQUIRED)
set(LIBRARIES Boost::boost)

if (PRECISION STREQUAL "FLOAT64")
    add_definitions(-DFLOAT64)
elseif (PRECISION STREQUAL "FLOAT128")
    set(LIBRARIES ${LIBRARIES} quadmath)
    add_definitions(-DFLOAT128)
elseif (PRECISION STREQUAL "BIGFLOAT")
    find_package(GMP REQUIRED)
    find_package(MPFR REQUIRED)
    find_package(MPC REQUIRED)
    include_directories(${GMP_INCLUDES} ${MPFR_INCLUDES} ${MPC_INCLUDES})
    set(LIBRARIES ${LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})
    add_definitions(-DBIGFLOAT)
endif ()

include_directories(${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES src/forward.h src/forward.cpp)

add_executable(forward_ray_tracing src/main.cpp ${SOURCE_FILES})
target_link_libraries(forward_ray_tracing PRIVATE ${LIBRARIES})

if (ENABLE_TESTING)
add_executable(tests tests/test.cpp ${SOURCE_FILES})
target_link_libraries(tests PRIVATE ${LIBRARIES} Catch2::Catch2WithMain)
endif()

#find_package(Python REQUIRED COMPONENTS Interpreter Development)
#find_package(pybind11 CONFIG REQUIRED)

#pybind11_add_module(py_forward_ray_tracing pybind.cpp ${SOURCE_FILES})
#target_link_libraries(py_forward_ray_tracing PUBLIC ${LIBRARIES})
