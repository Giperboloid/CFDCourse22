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

# default input directory to "input_data/"
add_definitions(-DCFDLIB_INPUT_DIR_PATH="${CMAKE_SOURCE_DIR}/input_data/")

# default output directory
if (MSVS)
	# build/output for visual studio
	set(OUTPUT_DIR_PATH ${CMAKE_BINARY_DIR}/output/)
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="${OUTPUT_DIR_PATH}")
	file(MAKE_DIRECTORY ${OUTPUT_DIR_PATH})
else(MSVS)
	# ./ for linux
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="")
endif(MSVS)

add_subdirectory(lib)
add_subdirectory(progs)
