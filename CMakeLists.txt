cmake_minimum_required(VERSION 3.13)

project(ldpc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable (ldpc main.cpp fonctions.hpp fonctions.cpp)