cmake_minimum_required(VERSION 3.14)

project(matrix)

enable_testing()
add_subdirectory(unit_tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(matrix ./main.cpp)
target_include_directories(matrix PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(matrix PRIVATE -O2)
