cmake_minimum_required(VERSION 3.25)
project(InterpolationFD C)

set(CMAKE_C_STANDARD 17)

add_executable(InterpolationFD main.c main_aux.c main_aux.h)

find_package(GSL REQUIRED)

target_link_libraries(InterpolationFD GSL::gsl GSL::gslcblas)
