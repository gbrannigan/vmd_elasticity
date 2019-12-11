cmake_minimum_required(VERSION 3.15)
project(bending_modulus)

set(CMAKE_CXX_STANDARD 14)

find_package(TCL)

include_directories(.
        ${TCL_INCLUDE_PATH})

add_executable(bending_modulus
        bending_modulus.cpp
        floatarray2d.h)

target_compile_options(bending_modulus PUBLIC -Wall -g -O3)
target_link_libraries(bending_modulus fftw3f ${TCL_LIBRARY})