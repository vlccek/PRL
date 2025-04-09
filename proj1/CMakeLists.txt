cmake_minimum_required(VERSION 3.10)
project(PRL1)


# Povolit MPI
find_package(MPI REQUIRED)

# Nastavit kompilátor (jen pokud není detekován automaticky)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Přidat výchozí volby pro MPI programy
add_executable(oets oets.cpp)

# Přidat hlavičkové soubory a linkování s MPI
target_include_directories(oets PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(oets PRIVATE MPI::MPI_CXX)
