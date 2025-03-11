cmake_minimum_required(VERSION 3.14)

project(matrix)

enable_testing()
add_subdirectory(tests/unit_test)
add_subdirectory(tests/e2e_test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(matrix ./main.cpp)
target_include_directories(matrix PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
