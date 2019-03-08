cmake_minimum_required(VERSION 3.13)

project(ldpc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS "-fopenmp -O0 -ftree-parallelize-loops=6 -floop-parallelize-all")

add_executable (ldpc main.cpp fonctions.hpp fonctions.cpp)
#target_compile_definitions(ldpc PRIVATE _GLIBCXX_PARALLEL)