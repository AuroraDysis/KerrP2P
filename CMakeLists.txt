cmake_minimum_required(VERSION 3.17)
project(forward_ray_tracing)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost 1.82 REQUIRED)
find_package(Catch2 3 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/src)

set(LIBRARIES Boost::boost)
set(SOURCE_FILES src/forward.h src/forward.cpp)

add_executable(forward_ray_tracing src/main.cpp ${SOURCE_FILES})
target_link_libraries(forward_ray_tracing PRIVATE ${LIBRARIES})

add_executable(tests tests/test.cpp ${SOURCE_FILES})
target_link_libraries(tests PRIVATE ${LIBRARIES} Catch2::Catch2WithMain)

#find_package(Python REQUIRED COMPONENTS Interpreter Development)
#find_package(pybind11 CONFIG REQUIRED)

#pybind11_add_module(py_forward_ray_tracing pybind.cpp ${SOURCE_FILES})
#target_link_libraries(py_forward_ray_tracing PUBLIC ${LIBRARIES})
