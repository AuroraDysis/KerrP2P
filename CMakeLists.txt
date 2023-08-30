cmake_minimum_required(VERSION 3.17)
project(forward_ray_tracing)

set(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
find_package(TBB CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})

set(LIBRARIES Eigen3::Eigen hdf5::hdf5 hdf5::hdf5_cpp hdf5::hdf5_hl_cpp)
set(SOURCE_FILES src/main.cpp)
add_executable(forward_ray_tracing ${SOURCE_FILES})
target_link_libraries(forward_ray_tracing PUBLIC ${LIBRARIES} TBB::tbb)

#find_package(Python REQUIRED COMPONENTS Interpreter Development)
#find_package(pybind11 CONFIG REQUIRED)

#pybind11_add_module(py_forward_ray_tracing pybind.cpp ${SOURCE_FILES})
#target_link_libraries(py_forward_ray_tracing PUBLIC ${LIBRARIES})
