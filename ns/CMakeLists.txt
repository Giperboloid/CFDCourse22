cmake_minimum_required(VERSION 3.0)

Project(cfdapp)

set(CMAKE_CXX_STANDARD 17)

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${BOOST_ROOT})
include_directories(${AMGCL_ROOT})

set(CFDLIB_TARGET cfdlib)

add_subdirectory(lib)
add_subdirectory(progs)
