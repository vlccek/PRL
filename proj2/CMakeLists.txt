cmake_minimum_required(VERSION 3.31)
project(PRJ2)

# Povolit MPI
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_executable(vuv vuv.cpp)

target_include_directories(vuv PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(vuv PRIVATE MPI::MPI_CXX)

