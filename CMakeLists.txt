# Set minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set project name
project(MyProject)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Add executable target with test.cpp and benchmark.cpp
add_executable(my_executable sources/test.cpp sources/benchmark.cpp)

# Link OpenMP flags to target
target_link_libraries(my_executable PUBLIC OpenMP::OpenMP_CXX)

# Alternatively, set CMAKE_CXX_FLAGS to include -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Include headers folder
include_directories(headers)

# Alternatively, use target_include_directories for specific target
target_include_directories(my_executable PRIVATE sources headers)

# Add library target with benchmark.cpp
add_library(benchmark sources/benchmark.cpp)

# Link benchmark library to executable target
target_link_libraries(my_executable PUBLIC benchmark)
