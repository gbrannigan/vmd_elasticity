cmake_minimum_required(VERSION 3.15)
project(bending_modulus)

set(CMAKE_CXX_STANDARD 14)

find_package(TCL)

include_directories(.
        ${TCL_INCLUDE_PATH})

add_executable(bending_modulus
        commandline.cpp
        bending_modulus.cpp
        bending_modulus.h
        floatarray2d.h floatarray2d.cpp)

add_library(vmd_bending_modulus SHARED
        vmd_plugin_interface.cpp
        bending_modulus.cpp
        bending_modulus.h
        floatarray2d.h floatarray2d.cpp)

target_compile_options(bending_modulus PUBLIC -Wall -g -O3)
target_compile_options(vmd_bending_modulus PUBLIC -g)
target_link_libraries(bending_modulus fftw3f)
target_link_libraries(vmd_bending_modulus fftw3f ${TCL_LIBRARY})
